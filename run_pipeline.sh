#!/usr/bin/env bash
#
# Oxford Nanopore (MinION) reference-mapping pipeline
#
# Folder-name based matching: for every experiment folder that exists
# under BOTH references_dir and data_dir (same folder name, e.g.
# "20260610_pUC19_test"), the pipeline maps reads against the matching
# reference(s). Output files are named after the reference FASTA (the
# real vector/sample name), not the barcode ID, so no manual renaming
# is needed afterwards.
#
# Layout:
#   references/<date>_<experiment_name>/<reference_name>.fasta
#   data/<date>_<experiment_name>/...
#   results/<date>_<experiment_name>_<run_YYMMDD_HHMM>/<reference_name>/...
#
# Per-reference fastq matching, in order of preference:
#   1. data/<experiment>/**/<reference_name>.<ext>  (exact filename match)
#   2. data/<experiment>/<reference_name>/**/*.fastq(.gz)  (alias/name folder)
#   3. data/<experiment>/barcodeNN/  (reference filename starts with NN_...)
#   4. shared pool (single-reference-per-experiment fallback)
#
# Consensus outputs:
#   <reference_name>.samtools.consensus.fasta
#       Raw read-pileup consensus from samtools consensus. Useful for checking
#       ambiguous N sites, low-confidence regions, and ONT homopolymer behavior.
#   <reference_name>.consensus.fasta
#       Final SnapGene/review consensus. Built by applying high-confidence VCF
#       variants to the original reference using bcftools consensus. Sites not
#       called as variants remain as the reference base, reducing false N/indel
#       artifacts in homopolymer/low-complexity ONT regions.
#
# Usage:
#   ./run_pipeline.sh [config.yaml]
#
# Requires: minimap2, samtools, bcftools (or medaka if variant_caller=medaka)
#           NanoFilt (optional, only used if min_read_length/quality > 0)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-$SCRIPT_DIR/config.yaml}"

# ---- minimal flat YAML reader (key: value, no nesting) -------------------
read_cfg() {
    local key="$1" default="$2"
    local val
    val=$(grep -E "^${key}:" "$CONFIG" | head -n1 | sed -E "s/^${key}:[[:space:]]*//; s/[[:space:]]*(#.*)?$//")
    val=$(sed -E 's/^"(.*)"$/\1/; s/^'"'"'(.*)'"'"'$/\1/' <<< "$val")
    val=$(sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//' <<< "$val")
    if [[ -z "$val" ]]; then echo "$default"; else echo "$val"; fi
}

resolve_cfg_path() {
    local value="$1"
    if [[ "$value" == /* ]]; then
        echo "$value"
    else
        echo "$SCRIPT_DIR/$value"
    fi
}

# Relative paths remain relative to the repository for backwards
# compatibility. Absolute paths are accepted so the local UI can create an
# isolated job without modifying the normal references/data/results folders.
REF_ROOT="$(resolve_cfg_path "$(read_cfg references_dir references)")"
DATA_ROOT="$(resolve_cfg_path "$(read_cfg data_dir data)")"
RESULTS_ROOT="$(resolve_cfg_path "$(read_cfg results_dir results)")"
PRESET="$(read_cfg minimap2_preset map-ont)"
THREADS="$(read_cfg threads 4)"
MIN_LEN="$(read_cfg min_read_length 0)"
MIN_QUAL="$(read_cfg min_read_quality 0)"
VARIANT_CALLER="$(read_cfg variant_caller bcftools)"
PARALLEL_JOBS="$(read_cfg parallel_jobs 1)"
PILON_JAR="$(read_cfg pilon_jar "")"
if [[ -n "$PILON_JAR" && "$PILON_JAR" != /* ]]; then
    PILON_JAR="$SCRIPT_DIR/$PILON_JAR"
fi
PILON_MEM="$(read_cfg pilon_mem 16G)"

RUN_TIMESTAMP="$(date +%y%m%d_%H%M)"

# Prefer the modern ONT SUP profile when the installed bcftools provides it.
# The legacy `-X ont` profile contains `-I` (skip indels), so it must not be
# used for a plasmid SNP/indel analyser. Older versions use compatible options
# which keep indel calling enabled.
BCFTOOLS_PLATFORM_MODE="compatible"
BCFTOOLS_PROFILE="ONT compatible (-B -Q5)"
if command -v bcftools >/dev/null 2>&1; then
    BCFTOOLS_CONFIGS="$(bcftools mpileup -X list 2>&1 || true)"
    if grep -qE '(^|[[:space:],])ont-sup([[:space:],]|$)' <<< "$BCFTOOLS_CONFIGS"; then
        BCFTOOLS_PLATFORM_MODE="ont-sup"
        BCFTOOLS_PROFILE="ont-sup"
    fi
fi

echo "==================================================================="
echo " Nanopore reference-mapping pipeline"
echo "   references  : $REF_ROOT"
echo "   data        : $DATA_ROOT"
echo "   results     : $RESULTS_ROOT"
echo "   preset      : $PRESET   threads/sample: $THREADS   parallel samples: $PARALLEL_JOBS"
echo "   bcftools    : $BCFTOOLS_PROFILE"
echo "==================================================================="

mkdir -p "$RESULTS_ROOT"
shopt -s nullglob

# ---------------------------------------------------------------------------
# Interactive experiment selection
# ---------------------------------------------------------------------------
AVAILABLE_EXPS=()
for EXP_DIR in "$REF_ROOT"/*/; do
    EXP_NAME="$(basename "$EXP_DIR")"
    if [[ -d "$DATA_ROOT/$EXP_NAME" ]]; then
        AVAILABLE_EXPS+=("$EXP_NAME")
    fi
done

if [[ ${#AVAILABLE_EXPS[@]} -eq 0 ]]; then
    echo "No experiments found (need a folder with the same name under both" >&2
    echo "  $REF_ROOT/ and $DATA_ROOT/)." >&2
    exit 1
fi

SELECTED_EXPS=()
if [[ -t 0 ]]; then
    echo "Available experiments:"
    for i in "${!AVAILABLE_EXPS[@]}"; do
        printf '  %2d) %s\n' "$((i + 1))" "${AVAILABLE_EXPS[$i]}"
    done
    echo ""
    read -r -p "Run which experiment(s)? (number(s), comma-separated, or Enter for all): " SELECTION
    if [[ -z "$SELECTION" ]]; then
        SELECTED_EXPS=("${AVAILABLE_EXPS[@]}")
    else
        IFS=',' read -r -a CHOICES <<< "$SELECTION"
        for choice in "${CHOICES[@]}"; do
            choice="$(tr -d '[:space:]' <<< "$choice")"
            [[ -z "$choice" ]] && continue
            if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#AVAILABLE_EXPS[@]} )); then
                SELECTED_EXPS+=("${AVAILABLE_EXPS[$((choice - 1))]}")
            else
                for e in "${AVAILABLE_EXPS[@]}"; do
                    [[ "$e" == "$choice" ]] && SELECTED_EXPS+=("$e")
                done
            fi
        done
        if [[ ${#SELECTED_EXPS[@]} -eq 0 ]]; then
            echo "No valid experiment selected, exiting." >&2
            exit 1
        fi
    fi
else
    SELECTED_EXPS=("${AVAILABLE_EXPS[@]}")
fi

echo ""
echo "Selected experiment(s): ${SELECTED_EXPS[*]}"

# ---------------------------------------------------------------------------
# process_one EXP_DIR REF_PATH
# ---------------------------------------------------------------------------
process_one() {
    set -euo pipefail
    local EXP_DIR="$1" REF_PATH="$2"
    local EXP_NAME DATA_EXP_DIR EXP_RESULTS_DIR REF_FILE REF_NAME SAMPLE_DIR LOG

    EXP_NAME="$(basename "$EXP_DIR")"
    DATA_EXP_DIR="$DATA_ROOT/$EXP_NAME"
    EXP_RESULTS_DIR="$RESULTS_ROOT/${EXP_NAME}_${RUN_TIMESTAMP}"
    REF_FILE="$(basename "$REF_PATH")"
    REF_NAME="${REF_FILE%.*}"
    SAMPLE_DIR="$EXP_RESULTS_DIR/$REF_NAME"
    LOG="[$EXP_NAME/$REF_NAME]"
    mkdir -p "$SAMPLE_DIR"

    echo "$LOG --- Reference: $REF_FILE -> results/${EXP_NAME}_${RUN_TIMESTAMP}/$REF_NAME/ ---"

    local REF_FILES_ALL=("$EXP_DIR"*.fasta "$EXP_DIR"*.fa "$EXP_DIR"*.fna)
    local REF_NAMES=() rp rf
    for rp in "${REF_FILES_ALL[@]}"; do
        rf="$(basename "$rp")"
        REF_NAMES+=("${rf%.*}")
    done

    # Batch UI stages uploaded reads as symbolic links. Accept both regular
    # files and links so those reads are not skipped before mapping.
    local READ_FIND_EXPR=( \( -type f -o -type l \) \( -iname "*.fastq.gz" -o -iname "*.fastq" -o -iname "*.fq.gz" -o -iname "*.fq" -o -iname "*.gz" \) )

    strip_read_ext() {
        local n="$1"
        n="${n%.gz}"
        n="${n%.fastq}"
        n="${n%.fq}"
        echo "$n"
    }

    local SHARED_FASTQ=() f rel top is_ref_dir is_ref_file r
    while IFS= read -r -d '' f; do
        rel="${f#"$DATA_EXP_DIR"/}"
        top="${rel%%/*}"
        is_ref_dir=false
        is_ref_file=false
        for r in "${REF_NAMES[@]}"; do
            if [[ "$top" == "$r" ]]; then is_ref_dir=true; fi
            if [[ "$(strip_read_ext "$(basename "$f")")" == "$r" ]]; then is_ref_file=true; fi
        done
        if [[ "$is_ref_dir" == false && "$is_ref_file" == false && ! "$top" =~ ^[Bb]arcode0*[0-9]+$ ]]; then
            SHARED_FASTQ+=("$f")
        fi
    done < <(find "$DATA_EXP_DIR" "${READ_FIND_EXPR[@]}" -print0)

    local FASTQ_FILES=()
    while IFS= read -r -d '' f; do
        if [[ "$(strip_read_ext "$(basename "$f")")" == "$REF_NAME" ]]; then
            FASTQ_FILES+=("$f")
        fi
    done < <(find "$DATA_EXP_DIR" "${READ_FIND_EXPR[@]}" -print0)

    local REF_SUBDIR="$DATA_EXP_DIR/$REF_NAME"
    if [[ ${#FASTQ_FILES[@]} -gt 0 ]]; then
        echo "$LOG Reads: ${#FASTQ_FILES[@]} file(s) matched by filename '$REF_NAME.*'"
    elif [[ -d "$REF_SUBDIR" ]]; then
        while IFS= read -r -d '' f; do FASTQ_FILES+=("$f"); done \
            < <(find "$REF_SUBDIR" "${READ_FIND_EXPR[@]}" -print0)
        echo "$LOG Reads: ${#FASTQ_FILES[@]} file(s) from $(basename "$REF_SUBDIR")/ (matched by name)"
    elif [[ "$REF_NAME" =~ ^0*([0-9]+)[^0-9] ]]; then
        local REF_NUM="${BASH_REMATCH[1]}" BARCODE_SUBDIR="" d dname
        for d in "$DATA_EXP_DIR"/*/; do
            dname="$(basename "$d")"
            if [[ "$dname" =~ ^[Bb]arcode0*([0-9]+)$ && "${BASH_REMATCH[1]}" == "$REF_NUM" ]]; then
                BARCODE_SUBDIR="$d"
                break
            fi
        done
        if [[ -n "$BARCODE_SUBDIR" ]]; then
            while IFS= read -r -d '' f; do FASTQ_FILES+=("$f"); done \
                < <(find "$BARCODE_SUBDIR" "${READ_FIND_EXPR[@]}" -print0)
            echo "$LOG Reads: ${#FASTQ_FILES[@]} file(s) from $(basename "$BARCODE_SUBDIR")/ (matched by barcode number $REF_NUM)"
        else
            FASTQ_FILES=("${SHARED_FASTQ[@]}")
            echo "$LOG Reads: ${#FASTQ_FILES[@]} file(s) from $EXP_NAME/ (shared, single-reference layout)"
        fi
    else
        FASTQ_FILES=("${SHARED_FASTQ[@]}")
        echo "$LOG Reads: ${#FASTQ_FILES[@]} file(s) from $EXP_NAME/ (shared, single-reference layout)"
    fi

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "$LOG [SKIP] no fastq files found for reference $REF_NAME" >&2
        return 0
    fi

    # 1. Merge reads
    local MERGED_FASTQ="$SAMPLE_DIR/${REF_NAME}.merged.fastq.gz"
    if [[ ! -s "$MERGED_FASTQ" ]]; then
        echo "$LOG [1/6] Merging ${#FASTQ_FILES[@]} read file(s) -> $(basename "$MERGED_FASTQ")"
        : > "${MERGED_FASTQ%.gz}.tmp"
        for f in "${FASTQ_FILES[@]}"; do
            if [[ "$f" == *.gz ]]; then
                zcat "$f" >> "${MERGED_FASTQ%.gz}.tmp"
            else
                cat "$f" >> "${MERGED_FASTQ%.gz}.tmp"
            fi
        done
        gzip -c "${MERGED_FASTQ%.gz}.tmp" > "$MERGED_FASTQ"
        rm -f "${MERGED_FASTQ%.gz}.tmp"
    else
        echo "$LOG [1/6] Merged reads already exist, skipping"
    fi

    # 2. Optional QC filtering with NanoFilt
    local READS_FOR_MAPPING="$MERGED_FASTQ"
    if [[ "$MIN_LEN" -gt 0 || "$MIN_QUAL" -gt 0 ]]; then
        if command -v NanoFilt >/dev/null 2>&1; then
            local FILTERED_FASTQ="$SAMPLE_DIR/${REF_NAME}.filtered.fastq.gz"
            echo "$LOG [2/6] QC filtering (length>=$MIN_LEN, quality>=$MIN_QUAL) -> $(basename "$FILTERED_FASTQ")"
            zcat "$MERGED_FASTQ" | NanoFilt -l "$MIN_LEN" -q "$MIN_QUAL" | gzip > "$FILTERED_FASTQ"
            READS_FOR_MAPPING="$FILTERED_FASTQ"
        else
            echo "$LOG [2/6] WARNING: NanoFilt not found; requested QC filtering was skipped" >&2
        fi
    else
        echo "$LOG [2/6] QC filtering disabled"
    fi

    # 3. Map to reference with minimap2, sort/index with samtools
    local SORTED_BAM="$SAMPLE_DIR/${REF_NAME}.sorted.bam"
    echo "$LOG [3/6] Mapping -> $(basename "$SORTED_BAM")"
    minimap2 -ax "$PRESET" -t "$THREADS" "$REF_PATH" "$READS_FOR_MAPPING" \
        | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
    samtools index "$SORTED_BAM"
    samtools flagstat "$SORTED_BAM" > "$SAMPLE_DIR/${REF_NAME}.flagstat.txt"
    samtools depth -a "$SORTED_BAM" > "$SAMPLE_DIR/${REF_NAME}.depth.txt"

    # 4. Raw samtools consensus: useful for ambiguous N inspection
    local SAMTOOLS_CONSENSUS_FASTA="$SAMPLE_DIR/${REF_NAME}.samtools.consensus.fasta"
    echo "$LOG [4/6] Building raw samtools consensus -> $(basename "$SAMTOOLS_CONSENSUS_FASTA")"
    samtools consensus -a -f fasta "$SORTED_BAM" > "$SAMTOOLS_CONSENSUS_FASTA"
    sed -i "1s/.*/>${REF_NAME}_samtools_consensus/" "$SAMTOOLS_CONSENSUS_FASTA"

    # 5. Variant calling
    local VCF="$SAMPLE_DIR/${REF_NAME}.vcf.gz"
    echo "$LOG [5/6] Calling variants ($VARIANT_CALLER) -> $(basename "$VCF")"
    if [[ "$VARIANT_CALLER" == "medaka" ]] && command -v medaka_haploid_variant >/dev/null 2>&1; then
        medaka_haploid_variant -i "$READS_FOR_MAPPING" -r "$REF_PATH" -o "$SAMPLE_DIR/medaka"
        cp "$SAMPLE_DIR/medaka/medaka.annotated.vcf" "$SAMPLE_DIR/${REF_NAME}.vcf" 2>/dev/null || true
        bgzip -f "$SAMPLE_DIR/${REF_NAME}.vcf"
        bcftools index -f "$VCF"
    elif [[ "$VARIANT_CALLER" == "pilon" ]]; then
        if [[ -z "$PILON_JAR" || ! -f "$PILON_JAR" ]]; then
            echo "$LOG [SKIP] variant_caller=pilon but pilon_jar not found: '$PILON_JAR'" >&2
        else
            local PILON_DIR="$SAMPLE_DIR/pilon"
            mkdir -p "$PILON_DIR"
            java -Xmx"$PILON_MEM" -jar "$PILON_JAR" \
                --genome "$REF_PATH" --bam "$SORTED_BAM" \
                --output "$REF_NAME" --outdir "$PILON_DIR" \
                --fix all --mindepth 2.0 --changes --vcf --verbose --threads "$THREADS" \
                > "$PILON_DIR/${REF_NAME}.pilon.log" 2>&1

            if [[ -f "$PILON_DIR/${REF_NAME}.changes" ]]; then
                cp "$PILON_DIR/${REF_NAME}.changes" "$SAMPLE_DIR/${REF_NAME}.changes"
            fi

            if [[ -f "$PILON_DIR/${REF_NAME}.vcf" ]]; then
                bcftools view -e 'ALT="."' "$PILON_DIR/${REF_NAME}.vcf" -Oz -o "$VCF"
                bcftools index -f "$VCF"
            fi
        fi
    else
        # FORMAT/DP and FORMAT/AD let the UI calculate per-call depth and
        # alternate-allele fraction instead of relying on QUAL alone.
        local BCFTOOLS_PLATFORM_ARGS=(-B -Q 5)
        if [[ "$BCFTOOLS_PLATFORM_MODE" == "ont-sup" ]]; then
            BCFTOOLS_PLATFORM_ARGS=(-X ont-sup)
        fi
        bcftools mpileup "${BCFTOOLS_PLATFORM_ARGS[@]}" -Ou \
            -a FORMAT/DP,FORMAT/AD -f "$REF_PATH" "$SORTED_BAM" \
            | bcftools call --ploidy 1 -mv -Oz -o "$VCF"
        bcftools index -f "$VCF"
    fi

    # 6. Final reference-guided consensus for SnapGene/review
    #    This applies only called VCF variants to the reference. Non-variant
    #    ambiguous/homopolymer sites remain as the reference base.
    local FINAL_CONSENSUS_FASTA="$SAMPLE_DIR/${REF_NAME}.consensus.fasta"
    if [[ -s "$VCF" ]]; then
        echo "$LOG [6/6] Building final VCF-guided consensus -> $(basename "$FINAL_CONSENSUS_FASTA")"
        bcftools consensus -f "$REF_PATH" "$VCF" > "$FINAL_CONSENSUS_FASTA"
        sed -i "1s/.*/>${REF_NAME}/" "$FINAL_CONSENSUS_FASTA"
    else
        echo "$LOG [6/6] No VCF found; copying reference as final consensus -> $(basename "$FINAL_CONSENSUS_FASTA")"
        cp "$REF_PATH" "$FINAL_CONSENSUS_FASTA"
        sed -i "1s/.*/>${REF_NAME}/" "$FINAL_CONSENSUS_FASTA"
    fi

    # Per-sample summary report
    local REPORT="$SAMPLE_DIR/${REF_NAME}_report.md"
    local MAPPED MEAN_DEPTH N_COUNT
    MAPPED=$(grep "mapped (" "$SAMPLE_DIR/${REF_NAME}.flagstat.txt" | head -1)
    MEAN_DEPTH=$(awk '{sum+=$3; n++} END {if (n>0) printf "%.2f", sum/n; else print "0"}' "$SAMPLE_DIR/${REF_NAME}.depth.txt")
    N_COUNT=$(awk '!/^>/{line=$0; gsub(/[^Nn]/, "", line); n+=length(line)} END{print n+0}' "$SAMTOOLS_CONSENSUS_FASTA")
    {
        echo "# Report: $REF_NAME"
        echo ""
        echo "- Experiment: $EXP_NAME"
        echo "- Reference: $REF_FILE"
        echo "- Mapping: $MAPPED"
        echo "- Mean depth: ${MEAN_DEPTH}x"
        echo "- Final consensus for SnapGene: ${REF_NAME}.consensus.fasta"
        echo "- Raw samtools consensus: ${REF_NAME}.samtools.consensus.fasta"
        echo "- Raw samtools consensus N count: ${N_COUNT}"
        echo "- Variants: ${REF_NAME}.vcf.gz"
    } > "$REPORT"

    echo "$LOG Done: $SAMPLE_DIR"
}

# ---------------------------------------------------------------------------
# Build the job list
# ---------------------------------------------------------------------------
JOBLIST="$(mktemp)"
trap 'rm -f "$JOBLIST"' EXIT

JOB_COUNT=0
for EXP_NAME in "${SELECTED_EXPS[@]}"; do
    EXP_DIR="$REF_ROOT/$EXP_NAME/"
    DATA_EXP_DIR="$DATA_ROOT/$EXP_NAME"

    echo ""
    echo "==================================================================="
    echo "Experiment: $EXP_NAME"
    echo "==================================================================="

    if [[ ! -d "$DATA_EXP_DIR" ]]; then
        echo "  [SKIP] no matching data folder: $DATA_EXP_DIR" >&2
        continue
    fi

    REF_FILES=("$EXP_DIR"*.fasta "$EXP_DIR"*.fa "$EXP_DIR"*.fna)
    if [[ ${#REF_FILES[@]} -eq 0 ]]; then
        echo "  [SKIP] no reference fasta found in $EXP_DIR" >&2
        continue
    fi

    mkdir -p "$RESULTS_ROOT/${EXP_NAME}_${RUN_TIMESTAMP}"

    for REF_PATH in "${REF_FILES[@]}"; do
        printf '%s\0%s\0' "$EXP_DIR" "$REF_PATH" >> "$JOBLIST"
        JOB_COUNT=$((JOB_COUNT + 1))
    done
done

echo ""
echo "==================================================================="
echo " Running $JOB_COUNT sample(s), $PARALLEL_JOBS at a time"
echo "==================================================================="

export DATA_ROOT RESULTS_ROOT PRESET THREADS MIN_LEN MIN_QUAL VARIANT_CALLER BCFTOOLS_PLATFORM_MODE PILON_JAR PILON_MEM RUN_TIMESTAMP
export -f process_one

if [[ "$JOB_COUNT" -gt 0 ]]; then
    if ! xargs -0 -n2 -P "$PARALLEL_JOBS" bash -c 'process_one "$1" "$2"' _ < "$JOBLIST"; then
        echo "" >&2
        echo "WARNING: one or more samples failed - see [SKIP]/error messages above" >&2
    fi
fi

echo ""
echo "Pipeline finished. Results in: $RESULTS_ROOT"

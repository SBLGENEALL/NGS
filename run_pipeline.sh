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
#   results/<date>_<experiment_name>/<reference_name>/...
#
# Per-reference fastq matching, in order of preference:
#   1. data/<experiment>/**/<reference_name>.<ext>  (exact filename match,
#      ext can be .fastq(.gz), .fq(.gz) or plain .gz, e.g.
#      R34.141-DAR-revSPG23_TIR90.fa <-> R34.141-DAR-revSPG23_TIR90.gz)
#   2. data/<experiment>/<reference_name>/**/*.fastq(.gz)
#      (e.g. demultiplexed with a sample sheet whose alias = reference_name
#      -- see scripts/generate_samplesheet.py)
#   3. data/<experiment>/barcodeNN/  (reference filename starts with a
#      number NN that matches barcode number NN)
#   4. Any fastq(.gz) directly under data/<experiment>/ that is NOT
#      inside another reference-named subfolder/file or a barcodeNN/
#      folder (single-reference-per-experiment case, e.g.
#      data/<experiment>/barcode01/*.fastq.gz)
#
# Parallelism:
#   Each reference (sample) is processed independently, so on a
#   high-core workstation you can process many samples at once.
#   Set "threads" (per-sample minimap2/samtools threads) and
#   "parallel_jobs" (how many samples to process concurrently) in
#   config.yaml so that threads * parallel_jobs <= number of CPU cores.
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
    # Strip surrounding quotes, e.g. pilon_jar: "/path/to/pilon.jar"
    val=$(sed -E 's/^"(.*)"$/\1/; s/^'"'"'(.*)'"'"'$/\1/' <<< "$val")
    if [[ -z "$val" ]]; then echo "$default"; else echo "$val"; fi
}

REF_ROOT="$SCRIPT_DIR/$(read_cfg references_dir references)"
DATA_ROOT="$SCRIPT_DIR/$(read_cfg data_dir data)"
RESULTS_ROOT="$SCRIPT_DIR/$(read_cfg results_dir results)"
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

# bcftools mpileup supports a "-X ont" config preset that tunes indel/SNP
# calling parameters for ONT's error profile (much fewer false-positive
# indels in homopolymer runs than the default/Illumina-tuned settings).
BCFTOOLS_PLATFORM_OPT=""
if command -v bcftools >/dev/null 2>&1 && bcftools mpileup 2>&1 | grep -qE -- '-X|--config'; then
    BCFTOOLS_PLATFORM_OPT="-X ont"
fi

echo "==================================================================="
echo " Nanopore reference-mapping pipeline"
echo "   references  : $REF_ROOT"
echo "   data        : $DATA_ROOT"
echo "   results     : $RESULTS_ROOT"
echo "   preset      : $PRESET   threads/sample: $THREADS   parallel samples: $PARALLEL_JOBS"
echo "==================================================================="

mkdir -p "$RESULTS_ROOT"

shopt -s nullglob

# ---------------------------------------------------------------------------
# process_one EXP_DIR REF_PATH
#
# Runs the full per-sample pipeline (merge -> QC -> mapping -> consensus ->
# variant calling -> report) for a single reference fasta. Self-contained so
# it can be exported and run in parallel via xargs.
# ---------------------------------------------------------------------------
process_one() {
    set -euo pipefail
    local EXP_DIR="$1" REF_PATH="$2"
    local EXP_NAME DATA_EXP_DIR EXP_RESULTS_DIR REF_FILE REF_NAME SAMPLE_DIR LOG

    EXP_NAME="$(basename "$EXP_DIR")"
    DATA_EXP_DIR="$DATA_ROOT/$EXP_NAME"
    EXP_RESULTS_DIR="$RESULTS_ROOT/$EXP_NAME"
    REF_FILE="$(basename "$REF_PATH")"
    REF_NAME="${REF_FILE%.*}"
    SAMPLE_DIR="$EXP_RESULTS_DIR/$REF_NAME"
    LOG="[$EXP_NAME/$REF_NAME]"
    mkdir -p "$SAMPLE_DIR"

    echo "$LOG --- Reference: $REF_FILE -> results/$EXP_NAME/$REF_NAME/ ---"

    # All reference names in this experiment, used to detect per-reference
    # fastq subfolders/files vs. the shared "leftover" pool below.
    local REF_FILES_ALL=("$EXP_DIR"*.fasta "$EXP_DIR"*.fa "$EXP_DIR"*.fna)
    local REF_NAMES=() rp rf
    for rp in "${REF_FILES_ALL[@]}"; do
        rf="$(basename "$rp")"
        REF_NAMES+=("${rf%.*}")
    done

    # Read-file patterns: ONT exports are sometimes named *.fastq.gz / *.fastq,
    # but also just *.fq(.gz) or even plain *.gz with no "fastq" in the name.
    local READ_FIND_EXPR=( -type f \( -iname "*.fastq.gz" -o -iname "*.fastq" -o -iname "*.fq.gz" -o -iname "*.fq" -o -iname "*.gz" \) )

    # Strip read-file extensions to compare against a reference name
    # (handles "<name>.gz", "<name>.fastq.gz", "<name>.fq", ...).
    strip_read_ext() {
        local n="$1"
        n="${n%.gz}"
        n="${n%.fastq}"
        n="${n%.fq}"
        echo "$n"
    }

    # Fastq files directly under DATA_EXP_DIR that are NOT inside a
    # subfolder named after one of this experiment's references, not inside
    # a barcodeNN/ folder (reserved for number-based matching), and whose
    # name doesn't exactly match one of this experiment's reference names
    # (reserved for exact-file matching). Used as the fallback pool for
    # references that don't match anything else (e.g. a single reference
    # per experiment with data/<exp>/barcode01/).
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

    # Pick the fastq files for this reference, in order of preference:
    #   1. data/<exp>/**/<reference_name>.<ext>  (exact filename match,
    #      e.g. R34.141-DAR-revSPG23_TIR90.fa <-> R34.141-DAR-revSPG23_TIR90.gz)
    #   2. data/<exp>/<reference_name>/          (exact name folder match)
    #   3. data/<exp>/barcodeNN/                  (reference filename
    #      starts with a number "NN_..." that matches barcode number NN)
    #   4. shared pool (single-reference-per-experiment fallback)
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
        echo "$LOG [1/5] Merging ${#FASTQ_FILES[@]} read file(s) -> $(basename "$MERGED_FASTQ")"
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
        echo "$LOG [1/5] Merged reads already exist, skipping"
    fi

    # 2. Optional QC filtering with NanoFilt
    local READS_FOR_MAPPING="$MERGED_FASTQ"
    if [[ "$MIN_LEN" -gt 0 || "$MIN_QUAL" -gt 0 ]] && command -v NanoFilt >/dev/null 2>&1; then
        local FILTERED_FASTQ="$SAMPLE_DIR/${REF_NAME}.filtered.fastq.gz"
        echo "$LOG [2/5] QC filtering (length>=$MIN_LEN, quality>=$MIN_QUAL) -> $(basename "$FILTERED_FASTQ")"
        zcat "$MERGED_FASTQ" | NanoFilt -l "$MIN_LEN" -q "$MIN_QUAL" | gzip > "$FILTERED_FASTQ"
        READS_FOR_MAPPING="$FILTERED_FASTQ"
    else
        echo "$LOG [2/5] QC filtering skipped"
    fi

    # 3. Map to reference with minimap2, sort with samtools
    local SORTED_BAM="$SAMPLE_DIR/${REF_NAME}.sorted.bam"
    echo "$LOG [3/5] Mapping -> $(basename "$SORTED_BAM")"
    minimap2 -ax "$PRESET" -t "$THREADS" "$REF_PATH" "$READS_FOR_MAPPING" \
        | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
    samtools index "$SORTED_BAM"
    samtools flagstat "$SORTED_BAM" > "$SAMPLE_DIR/${REF_NAME}.flagstat.txt"
    samtools depth -a "$SORTED_BAM" > "$SAMPLE_DIR/${REF_NAME}.depth.txt"

    # 4. Consensus sequence
    local CONSENSUS_FASTA="$SAMPLE_DIR/${REF_NAME}.consensus.fasta"
    echo "$LOG [4/5] Building consensus -> $(basename "$CONSENSUS_FASTA")"
    samtools consensus -a -f fasta "$SORTED_BAM" > "$CONSENSUS_FASTA"
    sed -i "1s/.*/>${REF_NAME}/" "$CONSENSUS_FASTA"

    # 5. Variant calling
    local VCF="$SAMPLE_DIR/${REF_NAME}.vcf.gz"
    echo "$LOG [5/5] Calling variants ($VARIANT_CALLER) -> $(basename "$VCF")"
    if [[ "$VARIANT_CALLER" == "medaka" ]] && command -v medaka_haploid_variant >/dev/null 2>&1; then
        medaka_haploid_variant -i "$READS_FOR_MAPPING" -r "$REF_PATH" -o "$SAMPLE_DIR/medaka"
        cp "$SAMPLE_DIR/medaka/medaka.annotated.vcf" "$SAMPLE_DIR/${REF_NAME}.vcf" 2>/dev/null || true
        bgzip -f "$SAMPLE_DIR/${REF_NAME}.vcf"
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

            # .changes file: one line per correction (POS REF->ALT), copied
            # next to the other per-sample outputs for easy inspection.
            if [[ -f "$PILON_DIR/${REF_NAME}.changes" ]]; then
                cp "$PILON_DIR/${REF_NAME}.changes" "$SAMPLE_DIR/${REF_NAME}.changes"
            fi

            # Pilon's --vcf includes every reference position (ALT="." where
            # no change was made); keep only the actual variant records so
            # it matches the bcftools/medaka output format.
            if [[ -f "$PILON_DIR/${REF_NAME}.vcf" ]]; then
                bcftools view -e 'ALT="."' "$PILON_DIR/${REF_NAME}.vcf" -Oz -o "$VCF"
                bcftools index -f "$VCF"
            fi
        fi
    else
        bcftools mpileup $BCFTOOLS_PLATFORM_OPT -f "$REF_PATH" "$SORTED_BAM" 2>/dev/null \
            | bcftools call --ploidy 1 -mv -Oz -o "$VCF"
        bcftools index -f "$VCF"
    fi

    # Per-sample summary report
    local REPORT="$SAMPLE_DIR/${REF_NAME}_report.md"
    local MAPPED MEAN_DEPTH
    MAPPED=$(grep "mapped (" "$SAMPLE_DIR/${REF_NAME}.flagstat.txt" | head -1)
    MEAN_DEPTH=$(awk '{sum+=$3; n++} END {if (n>0) printf "%.2f", sum/n; else print "0"}' "$SAMPLE_DIR/${REF_NAME}.depth.txt")
    {
        echo "# Report: $REF_NAME"
        echo ""
        echo "- Experiment: $EXP_NAME"
        echo "- Reference: $REF_FILE"
        echo "- Mapping: $MAPPED"
        echo "- Mean depth: ${MEAN_DEPTH}x"
        echo "- Consensus: ${REF_NAME}.consensus.fasta"
        echo "- Variants: ${REF_NAME}.vcf.gz"
    } > "$REPORT"

    echo "$LOG Done: $SAMPLE_DIR"
}

# ---------------------------------------------------------------------------
# Build the job list: one (EXP_DIR, REF_PATH) pair per reference fasta of
# every experiment that has a matching data folder.
# ---------------------------------------------------------------------------
JOBLIST="$(mktemp)"
trap 'rm -f "$JOBLIST"' EXIT

JOB_COUNT=0
for EXP_DIR in "$REF_ROOT"/*/; do
    EXP_NAME="$(basename "$EXP_DIR")"
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

    mkdir -p "$RESULTS_ROOT/$EXP_NAME"

    for REF_PATH in "${REF_FILES[@]}"; do
        printf '%s\0%s\0' "$EXP_DIR" "$REF_PATH" >> "$JOBLIST"
        JOB_COUNT=$((JOB_COUNT + 1))
    done
done

echo ""
echo "==================================================================="
echo " Running $JOB_COUNT sample(s), $PARALLEL_JOBS at a time"
echo "==================================================================="

export DATA_ROOT RESULTS_ROOT PRESET THREADS MIN_LEN MIN_QUAL VARIANT_CALLER BCFTOOLS_PLATFORM_OPT PILON_JAR PILON_MEM
export -f process_one

if [[ "$JOB_COUNT" -gt 0 ]]; then
    if ! xargs -0 -n2 -P "$PARALLEL_JOBS" bash -c 'process_one "$1" "$2"' _ < "$JOBLIST"; then
        echo "" >&2
        echo "WARNING: one or more samples failed - see [SKIP]/error messages above" >&2
    fi
fi

echo ""
echo "Pipeline finished. Results in: $RESULTS_ROOT"

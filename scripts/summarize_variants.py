#!/usr/bin/env python3
"""
summarize_variants.py

Scan results/<experiment>/<reference_name>/ for each sample's
<reference_name>.vcf.gz and <reference_name>_report.md, and produce a
single summary table across ALL samples: mapping rate, mean depth, and a
short description of every variant found (point mutation / insertion /
deletion), so you don't have to open each VCF individually.

Usage:
    python scripts/summarize_variants.py --results results --output results/summary_report.md

    # CSV instead of markdown
    python scripts/summarize_variants.py --results results --output results/summary_report.csv

    # Filter out low-confidence calls (common with ONT homopolymer indel
    # noise from bcftools default settings)
    python scripts/summarize_variants.py --results results --output results/summary_report.md \
        --min-qual 20 --min-depth 10
"""
import argparse
import csv
import gzip
import os
import re


def read_vcf_variants(vcf_path, min_qual=0.0, min_depth=0):
    """Return a list of (pos, type, ref, alt) for each ALT allele in the VCF
    that passes the QUAL/DP thresholds."""
    variants = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            pos, ref, alt_field, qual_str, info = cols[1], cols[3], cols[4], cols[5], cols[7]

            try:
                qual = float(qual_str)
            except ValueError:
                qual = 0.0
            if qual < min_qual:
                continue

            dp_match = re.search(r"(?:^|;)DP=(\d+)", info)
            dp = int(dp_match.group(1)) if dp_match else 0
            if dp < min_depth:
                continue

            for alt in alt_field.split(","):
                if alt == "." or alt == "":
                    continue
                if len(ref) == 1 and len(alt) == 1:
                    vtype = "point mutation"
                elif len(alt) > len(ref):
                    vtype = "insertion"
                elif len(alt) < len(ref):
                    vtype = "deletion"
                else:
                    vtype = "complex"
                variants.append((pos, vtype, ref, alt))
    return variants


def read_changes_variants(changes_path):
    """Parse a Pilon *.changes file into the same (pos, type, ref, alt) tuples
    used for VCF variants, so format_variants() can be reused.

    Each line looks like:
        <old_locus> <new_locus> <old_seq> <new_seq>
    e.g.
        pUC19:1234-1234 pUC19:1234-1234 A C
    where old_seq/new_seq is "." for pure insertions/deletions.
    """
    variants = []
    with open(changes_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            if len(cols) < 4:
                continue
            old_locus, _new_locus, old_seq, new_seq = cols[0], cols[1], cols[2], cols[3]

            m = re.match(r"^[^:]+:(\d+)", old_locus)
            pos = m.group(1) if m else "?"

            ref = "" if old_seq == "." else old_seq
            alt = "" if new_seq == "." else new_seq

            if len(ref) == 1 and len(alt) == 1:
                vtype = "point mutation"
            elif len(alt) > len(ref):
                vtype = "insertion"
            elif len(alt) < len(ref):
                vtype = "deletion"
            else:
                vtype = "complex"
            variants.append((pos, vtype, ref if ref else ".", alt if alt else "."))
    return variants


def read_report_fields(report_path):
    """Pull 'Mapping' and 'Mean depth' lines out of a *_report.md file."""
    mapping, depth = "", ""
    if not os.path.isfile(report_path):
        return mapping, depth
    with open(report_path) as f:
        for line in f:
            m = re.match(r"-\s*Mapping:\s*(.*)", line)
            if m:
                mapping = m.group(1).strip()
            m = re.match(r"-\s*Mean depth:\s*(.*)", line)
            if m:
                depth = m.group(1).strip()
    return mapping, depth


def format_variants(variants):
    if not variants:
        return "None"
    parts = []
    for pos, vtype, ref, alt in variants:
        label = {"point mutation": "SNP", "insertion": "ins", "deletion": "del", "complex": "complex"}[vtype]
        parts.append(f"{label} @{pos} {ref}>{alt}")
    return "; ".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results", help="Results root directory (default: results)")
    parser.add_argument("--output", required=True, help="Output summary file (.md or .csv)")
    parser.add_argument("--min-qual", type=float, default=0.0,
                         help="Minimum VCF QUAL to count a variant (default: 0, no filter)")
    parser.add_argument("--min-depth", type=int, default=0,
                         help="Minimum read depth (INFO/DP) to count a variant (default: 0, no filter)")
    args = parser.parse_args()

    rows = []
    for exp_name in sorted(os.listdir(args.results)):
        exp_dir = os.path.join(args.results, exp_name)
        if not os.path.isdir(exp_dir):
            continue
        for sample_name in sorted(os.listdir(exp_dir)):
            sample_dir = os.path.join(exp_dir, sample_name)
            if not os.path.isdir(sample_dir):
                continue
            vcf_path = os.path.join(sample_dir, f"{sample_name}.vcf.gz")
            changes_path = os.path.join(sample_dir, f"{sample_name}.changes")
            report_path = os.path.join(sample_dir, f"{sample_name}_report.md")
            if not os.path.isfile(vcf_path) and not os.path.isfile(changes_path):
                continue

            mapping, depth = read_report_fields(report_path)

            if os.path.isfile(changes_path):
                # Pilon: .changes already lists only the corrections made,
                # so no QUAL/DP filtering applies.
                method = "pilon"
                variants = read_changes_variants(changes_path)
            else:
                method = "bcftools/medaka"
                variants = read_vcf_variants(vcf_path, min_qual=args.min_qual, min_depth=args.min_depth)

            n_snp = sum(1 for v in variants if v[1] == "point mutation")
            n_ins = sum(1 for v in variants if v[1] == "insertion")
            n_del = sum(1 for v in variants if v[1] == "deletion")

            rows.append({
                "experiment": exp_name,
                "sample": sample_name,
                "mapping": mapping,
                "mean_depth": depth,
                "method": method,
                "n_snp": n_snp,
                "n_ins": n_ins,
                "n_del": n_del,
                "variants": format_variants(variants),
            })

    if not rows:
        raise SystemExit(f"No *.vcf.gz found under {args.results}/<experiment>/<sample>/")

    if args.output.lower().endswith(".csv"):
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["experiment", "sample", "mapping", "mean_depth", "method",
                                                     "n_snp", "n_ins", "n_del", "variants"])
            writer.writeheader()
            writer.writerows(rows)
    else:
        with open(args.output, "w") as f:
            f.write("# Variant summary\n\n")
            f.write("| Experiment | Sample | Mapping | Mean depth | Method | SNP | Ins | Del | Variants |\n")
            f.write("|---|---|---|---|---|---|---|---|---|\n")
            for r in rows:
                f.write(f"| {r['experiment']} | {r['sample']} | {r['mapping']} | {r['mean_depth']} | {r['method']} "
                        f"| {r['n_snp']} | {r['n_ins']} | {r['n_del']} | {r['variants']} |\n")

    print(f"Wrote summary for {len(rows)} sample(s) to {args.output}")
    n_with_variants = sum(1 for r in rows if r["variants"] != "None")
    print(f"  {n_with_variants} sample(s) have at least one variant")


if __name__ == "__main__":
    main()

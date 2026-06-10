#!/usr/bin/env python3
"""
generate_samplesheet.py

Generate a MinKNOW sample sheet (CSV) BEFORE sequencing, so that each
barcode/well is pre-labelled with the alias of its matching reference
FASTA. When MinKNOW/Guppy/Dorado demultiplex the run using this sample
sheet, the output folders are already named after the reference
("alias") instead of "barcode01", "barcode02", ... — no renaming step
needed afterwards.

Each *.fasta/*.fa/*.fna file found in --references is treated as one
well/sample and assigned the next barcode number (01, 02, ...). The
"alias" column is set to the reference filename (without extension), so
it matches the reference name 1:1 — this is exactly the name that
run_pipeline.sh later uses for its results/<experiment>/<reference_name>/
output folder.

Usage:
    python scripts/generate_samplesheet.py \
        --references references/20260610_pUC19_test \
        --experiment-id 20260610_pUC19_test \
        --output references/20260610_pUC19_test/sample_sheet.csv

    # If your plate order doesn't match alphabetical filename order,
    # list the files explicitly in plate order with --order:
    python scripts/generate_samplesheet.py \
        --references references/20260610_pUC19_test \
        --experiment-id 20260610_pUC19_test \
        --order order.txt \
        --output references/20260610_pUC19_test/sample_sheet.csv

Workflow:
    1. Put one reference FASTA per well into references/<experiment>/.
    2. Run this script to create sample_sheet.csv.
    3. Load sample_sheet.csv into MinKNOW before starting the run
       (Start run -> Sample sheet -> Browse).
    4. Sequence + basecall as usual. Output folders are named by alias.
    5. Copy the basecalled output into data/<experiment>/ (same
       experiment folder name as references/<experiment>/) and run
       run_pipeline.sh — no manual renaming required.
"""
import argparse
import csv
import os

REF_EXTENSIONS = (".fasta", ".fa", ".fna")


def list_reference_names(references_dir):
    names = []
    for entry in sorted(os.listdir(references_dir)):
        if entry.lower().endswith(REF_EXTENSIONS):
            names.append(os.path.splitext(entry)[0])
    return names


def load_order(order_file, references_dir):
    names = []
    with open(order_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().endswith(REF_EXTENSIONS):
                line = os.path.splitext(line)[0]
            names.append(line)
    available = set(list_reference_names(references_dir))
    missing = [n for n in names if n not in available]
    if missing:
        raise SystemExit(f"Names in --order not found under {references_dir}: {missing}")
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--references", required=True,
                         help="Folder containing one reference FASTA per well (e.g. references/<experiment>/)")
    parser.add_argument("--experiment-id", required=True,
                         help="Experiment name, used as experiment_id (should match the references/<name> "
                              "and data/<name> folder name)")
    parser.add_argument("--output", required=True, help="Output sample sheet CSV path")
    parser.add_argument("--kit", default="SQK-NBD114-96",
                         help="Barcoding kit name, written as a comment in the output (default: SQK-NBD114-96)")
    parser.add_argument("--start-barcode", type=int, default=1,
                         help="First barcode number to assign (default: 1 -> barcode01)")
    parser.add_argument("--flow-cell-id", default="",
                         help="Optional flow cell ID to fill in every row")
    parser.add_argument("--type", default="test_sample",
                         help="MinKNOW 'type' column value (default: test_sample)")
    parser.add_argument("--order", default=None,
                         help="Optional text file listing reference file/sample names in plate order "
                              "(one per line). Default: alphabetical order of files in --references.")
    args = parser.parse_args()

    if args.order:
        names = load_order(args.order, args.references)
    else:
        names = list_reference_names(args.references)

    if not names:
        raise SystemExit(f"No {REF_EXTENSIONS} reference files found in {args.references}")

    n_barcodes = len(names) + args.start_barcode - 1
    if n_barcodes > 96:
        raise SystemExit(f"{n_barcodes} barcodes needed but a 96-barcode kit only supports up to 96.")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as f:
        f.write(f"# kit: {args.kit}\n")
        f.write(f"# Generated from references in: {args.references}\n")
        f.write("# Load this file into MinKNOW (Start run -> Sample sheet) BEFORE sequencing,\n")
        f.write("# so demultiplexed output folders are named by 'alias' instead of 'barcodeNN'.\n")
        writer = csv.writer(f)
        writer.writerow(["flow_cell_id", "position_id", "experiment_id", "sample_id", "alias", "type", "barcode"])
        for i, name in enumerate(names):
            barcode_num = args.start_barcode + i
            barcode = f"barcode{barcode_num:02d}"
            writer.writerow([args.flow_cell_id, "", args.experiment_id, name, name, args.type, barcode])

    print(f"Wrote {len(names)} rows to {args.output}")
    for i, name in enumerate(names):
        print(f"  barcode{args.start_barcode + i:02d} -> {name}")


if __name__ == "__main__":
    main()

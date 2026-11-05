#!/usr/bin/env python3
"""Parse Supplementary Tables 1 (clinical) and 2 (scDNA QC) of Lin et al., Nature 2024
from the supplementary PDF (41586_2024_8129_MOESM1_ESM.pdf) into TSV.

Requires `pdftotext` (poppler).

Usage:
  python3 scripts/00_metadata/parse_suppl_tables.py \
      --pdf data/metadata/41586_2024_8129_MOESM1_ESM.pdf --out-dir data/metadata
Outputs:
  paper_clinical.tsv  Suppl. Table 1
  paper_scdna_qc.tsv  Suppl. Table 2, plus aneuploid_frac = aneuploid / cells_pass_qc
"""
import argparse
import csv
import os
import re
import subprocess

CLIN_COLS = ["donor", "tissue_label", "age", "ethnicity", "menopause", "bmi",
             "parity", "hyperplasia", "metaplasia"]
QC_COLS = ["donor", "mean_total_reads", "mean_reads_retained", "coverage_x",
           "duplication_rate", "paired_end", "cells_sequenced", "cells_pass_qc",
           "cells_pass_qc_pct", "diploid_cells", "aneuploid_cells"]

CLIN_RE = re.compile(
    r"^\s*(P\d{2})\s+(P\d{2}[LR])\s+(\d+)\s+(.+?)\s{2,}(\S+)\s+([\d.]+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")
QC_RE = re.compile(
    r"^\s*(P\d{2})\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)%\s+([TF])\s+(\d+)\s+(\d+)"
    r"\s+([\d.]+)%\s+(\d+)\s+(\d+)\s*$")


def write(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    text = subprocess.run(["pdftotext", "-layout", args.pdf, "-"],
                          check=True, capture_output=True, text=True).stdout
    clin, qc = {}, {}
    for line in text.splitlines():
        if m := QC_RE.match(line):
            qc[m.group(1)] = list(m.groups())
        elif m := CLIN_RE.match(line):
            clin[m.group(1)] = list(m.groups())

    for d, row in qc.items():
        row.append(f"{int(row[10]) / int(row[7]):.4f}")   # aneuploid / cells_pass_qc

    os.makedirs(args.out_dir, exist_ok=True)
    write(f"{args.out_dir}/paper_clinical.tsv", CLIN_COLS, [clin[k] for k in sorted(clin)])
    write(f"{args.out_dir}/paper_scdna_qc.tsv", QC_COLS + ["aneuploid_frac"],
          [qc[k] for k in sorted(qc)])
    print(f"clinical rows: {len(clin)}  qc rows: {len(qc)}")


if __name__ == "__main__":
    main()

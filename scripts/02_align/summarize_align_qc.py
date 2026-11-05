#!/usr/bin/env python3
"""Per-cell alignment QC table from the align workflow outputs.

Reads qc/<cell>.cutadapt.txt, qc/<cell>.flagstat.txt and qc/<cell>.markdup.txt, and
optionally counts reads QDNAseq will use from bam/<cell>.bam (primary, mapped, not
duplicate, MAPQ >= 37: binReadCounts defaults). Expected input reads come from ENA
(cells.tsv, sum of read_count over the cell's runs; R1 only, so spots = reads).

Columns:
  reads_ena        sum of ENA read_count over the cell's runs
  reads_in         reads entering cutadapt (should equal reads_ena)
  reads_short      reads dropped by cutadapt -m
  mapped_pct       primary mapped / reads after trimming
  dup_pct          duplicate / examined reads (samtools markdup, single-end)
  usable_reads     primary, mapped, non-duplicate, MAPQ >= 37 (needs --count-usable)
  pass_<bin>kb     usable_reads >= 15 * 2 * usable_bins / 0.91 (diploid threshold)

Usage:
  python3 scripts/02_align/summarize_align_qc.py --run-dir <pubdir> \
      --cells data/metadata/cells.tsv --usable-bins data/metadata/usable_bins.tsv \
      --count-usable --out <table.tsv>
--count-usable runs `samtools view -c` (samtools on PATH, or --samtools "docker run ...").
"""
import argparse
import csv
import os
import re
import shlex
import subprocess
from collections import defaultdict


def grab(path, pattern, cast=int):
    with open(path) as f:
        m = re.search(pattern, f.read(), re.M)
    return cast(m.group(1).replace(",", "")) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="directory with bam/ and qc/")
    ap.add_argument("--cells", required=True)
    ap.add_argument("--usable-bins", required=True)
    ap.add_argument("--count-usable", action="store_true")
    ap.add_argument("--samtools", default="samtools", help="samtools command prefix")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ena = defaultdict(int)
    meta = {}
    with open(args.cells) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ena[r["cell_id"]] += int(r["read_count"])
            meta.setdefault(r["cell_id"], r)
    with open(args.usable_bins) as f:
        bins = {int(r["bin_kb"]): int(r["usable_bins"]) for r in csv.DictReader(f, delimiter="\t")}
    thresholds = {b: 15 * 2 * n / 0.91 for b, n in sorted(bins.items())}

    qc = os.path.join(args.run_dir, "qc")
    cells = sorted(f[:-len(".flagstat.txt")] for f in os.listdir(qc) if f.endswith(".flagstat.txt"))
    rows = []
    for c in cells:
        p = lambda ext: os.path.join(qc, f"{c}.{ext}")
        reads_in = grab(p("cutadapt.txt"), r"Total reads processed:\s+([\d,]+)")
        short = grab(p("cutadapt.txt"), r"Reads that were too short:\s+([\d,]+)")
        total = grab(p("flagstat.txt"), r"^(\d+) \+ \d+ in total")
        mapped = grab(p("flagstat.txt"), r"^(\d+) \+ \d+ primary mapped")
        examined = grab(p("markdup.txt"), r"^EXAMINED: (\d+)")
        dups = grab(p("markdup.txt"), r"^DUPLICATE TOTAL: (\d+)")
        row = {"cell_id": c, "donor": meta[c]["donor"], "coassay": meta[c]["coassay"],
               "runs": meta[c]["runs_per_cell"], "reads_ena": ena[c], "reads_in": reads_in,
               "reads_short": short, "mapped_pct": round(100 * mapped / total, 2),
               "dup_pct": round(100 * dups / examined, 2)}
        if args.count_usable:
            bam = os.path.join(args.run_dir, "bam", f"{c}.bam")
            cmd = shlex.split(args.samtools) + ["view", "-c", "-q", "37", "-F", "0xD04", bam]
            n = int(subprocess.run(cmd, check=True, capture_output=True, text=True).stdout)
            row["usable_reads"] = n
            for b, t in thresholds.items():
                row[f"pass_{b}kb"] = int(n >= t)
        rows.append(row)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} cells -> {args.out}")


if __name__ == "__main__":
    main()

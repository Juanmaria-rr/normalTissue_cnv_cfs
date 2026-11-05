#!/usr/bin/env python3
"""Build an alignment samplesheet (cell_id <TAB> runs) from the per-donor pilot tables.

Optionally draws a reproducible random subset of cells per donor (for test runs).
Cells keep all their runs, which the align workflow concatenates.

Usage:
  python3 scripts/02_align/make_samplesheet.py --pilot-dir data/pilot \
      --donors P08 P07 P36 P20 P18 --n-per-donor 4 --seed 1 --out data/pilot/samplesheet_test20.tsv
  (omit --n-per-donor to keep every cell)
"""
import argparse
import csv
import os
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot-dir", required=True)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--n-per-donor", type=int)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    rows = []
    for donor in args.donors:
        with open(os.path.join(args.pilot_dir, f"cells_{donor}.tsv")) as f:
            cells = list(csv.DictReader(f, delimiter="\t"))
        if args.n_per_donor:
            cells = rng.sample(cells, min(args.n_per_donor, len(cells)))
        rows += cells

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_id", "runs"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} cells -> {args.out}")


if __name__ == "__main__":
    main()

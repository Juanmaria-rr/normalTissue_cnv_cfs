#!/usr/bin/env python3
"""Select pilot cells from data/metadata/cells.tsv.

Inclusion: donor in --donors, imaging qc_flag == "Good", cell_flag == "yes",
not a control well. No depth filter here: the 15-nrpbpcc filter is applied after
alignment/dedup by fitting_absCN (--filters powered), since raw ENA counts overestimate
usable reads. The raw-depth pass flags (per-cell totals) are carried along for reference.

A cell may have been sequenced in more than one run (top-up sequencing of the same
library); its runs are listed together and must be concatenated before alignment.

Outputs (in --out-dir):
  pilot_runs.tsv           selected runs (all cells.tsv columns)
  cells_<donor>.tsv        cell_id <TAB> comma-separated SRR runs, one row per cell
  pilot_summary.tsv        per-donor counts (cells, not runs)

Usage:
  python3 scripts/00_metadata/select_pilot_cells.py \
      --cells data/metadata/cells.tsv --donors P08 P07 P36 P20 P18 --out-dir data/pilot
"""
import argparse
import csv
import os
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with open(args.cells) as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames
        rows = [r for r in reader if r["donor"] in args.donors]

    sel = [r for r in rows if r["qc_flag"] == "Good" and r["cell_flag"] == "yes"
           and r["control_well"] == "0"]

    os.makedirs(args.out_dir, exist_ok=True)
    with open(f"{args.out_dir}/pilot_runs.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sel)

    runs_by_cell = defaultdict(list)
    cell_row = {}
    for r in sel:
        runs_by_cell[r["cell_id"]].append(r["run_accession"])
        cell_row[r["cell_id"]] = r
    for d in args.donors:
        with open(f"{args.out_dir}/cells_{d}.tsv", "w") as f:
            f.write("cell_id\truns\n")
            for cid in sorted(c for c in runs_by_cell if cell_row[c]["donor"] == d):
                f.write(f"{cid}\t{','.join(sorted(runs_by_cell[cid]))}\n")

    depth_cols = [c for c in cols if c.startswith("pass_depth_")]
    summ = []
    for d in args.donors:
        cs = [cell_row[c] for c in runs_by_cell if cell_row[c]["donor"] == d]
        rs = [r for r in sel if r["donor"] == d]
        row = {"donor": d,
               "n_cells": len(cs),
               "n_runs": len(rs),
               "n_multi_run_cells": sum(int(r["runs_per_cell"]) > 1 for r in cs),
               "n_coassay": sum(r["coassay"] == "1" for r in cs),
               "n_atac_barcode": len({r["cell_id"] for r in rs if r["atac_barcode"]}),
               # sequenced bases, both mates (download size scales with this)
               "gbases": round(sum(int(r["read_count"]) * int(r["read_length"]) * 2
                                   for r in rs) / 1e9, 1)}
        for c in depth_cols:
            row[c.replace("pass_depth_", "raw_pass_")] = sum(r[c] == "1" for r in cs)
        summ.append(row)
    with open(f"{args.out_dir}/pilot_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summ)
    print(f"selected {len(runs_by_cell)} cells ({len(sel)} runs) across {len(args.donors)} donors")


if __name__ == "__main__":
    main()

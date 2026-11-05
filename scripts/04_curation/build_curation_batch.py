#!/usr/bin/env python3
"""Build a curation-app batch from the long relative-CN export (export_relcn.R).

Writes, into <out-dir>/<batch>/:
  <batch>_copyNumbersSegmented_<bin>kb.h5   one group per sample: bins/{chromosome,start,
                                            end,use,copynumber,segmented}, samples/total_reads
  <batch>_gridSearch.csv                    sample, ploidy, purity, clonality, powered, LOHok
  <batch>_cellinfo.tsv                      per-sample QC shown in the app (optional --qc)
  samples.txt                               checklist, in the (shuffled) curation order
and the blinding key to <key-dir>/<batch>_key.tsv (sample code -> cell_id), outside the
folder the curator opens. No downsampled file is written: single cells are not downsampled.

The grid reproduces pureSamples_absCN.R (fitting_absCN @56db0d7), so the app's heatmap
shows the same surface the pipeline searched:
  ploidy 1.4-8 (0.1), purity 0.90-1.00 (0.01)
  rel_ploidy = mean(copynumber over used bins); seqdepth = rel_ploidy / (p*q + 2(1-q))
  abs = (seg/seqdepth - 2(1-q)) / q; clonality = mean((abs - round(abs))^2)
  powered: (((2(1-q)+qp)/(pq))/q) * 15 * (2(1-q)+qp) * n_used_bins / 0.91 < total_reads
  LOHok:   sum(abs < 0.8) < 1e7 / (bin_kb * 100)   (upstream formula, kept as is)
Unlike upstream, every purity is kept (upstream keeps the best purity per ploidy).

Usage:
  python3 scripts/04_curation/build_curation_batch.py --export curation/staging/test20_200kb \
      --bin 200 --batch test20_200kb --out-dir curation/data --key-dir curation/keys \
      [--qc cn_qc_200kb.tsv] [--no-blind] [--seed 1]
"""
import argparse
import csv
import gzip
import os
import random
import string

import h5py
import numpy as np
import pandas as pd

PLOIDIES = np.round(np.arange(1.4, 8.0001, 0.1), 1)
PURITIES = np.round(np.arange(0.90, 1.0001, 0.01), 2)
QC_COLS = ["reads", "nrpbpcc", "noise_excess", "breakpoints", "frac_nondiploid", "rt_rho"]


def grid_for_cell(cn, seg, total_reads, n_used, bin_kb, downs=15):
    rel_ploidy = np.nanmean(cn)
    s = seg[np.isfinite(seg)]
    rows = []
    for p in PLOIDIES:
        for q in PURITIES:
            cellploidy = p * q + 2 * (1 - q)
            seqdepth = rel_ploidy / cellploidy
            ab = (s / seqdepth - 2 * (1 - q)) / q
            clon = float(np.mean((ab - np.round(ab)) ** 2))
            need = (((2 * (1 - q) + q * p) / (p * q)) / q) * downs * (2 * (1 - q) + q * p) * n_used / 0.91
            loh = int(np.sum(ab < 0.8) < 1e7 / (bin_kb * 100))
            rows.append((p, q, clon, int(need < total_reads), loh))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", required=True, help="prefix written by export_relcn.R")
    ap.add_argument("--bin", type=int, required=True)
    ap.add_argument("--batch", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--key-dir", required=True)
    ap.add_argument("--qc", help="cn_qc_<bin>kb.tsv from cn_qc_metrics.R")
    ap.add_argument("--no-blind", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    bins = pd.read_csv(f"{a.export}_bins.tsv.gz", sep="\t", dtype={"chromosome": str})
    cells = pd.read_csv(f"{a.export}_cells.tsv", sep="\t").set_index("cell_id")
    ids = list(cells.index)

    rng = random.Random(a.seed)
    if a.no_blind:
        code = {c: c for c in ids}
    else:
        alphabet = string.ascii_uppercase + string.digits
        codes = set()
        while len(codes) < len(ids):
            codes.add("".join(rng.choice(alphabet) for _ in range(6)))
        code = dict(zip(ids, sorted(codes)))
    order = ids[:]
    rng.shuffle(order)

    bdir = os.path.join(a.out_dir, a.batch)
    os.makedirs(bdir, exist_ok=True)
    os.makedirs(a.key_dir, exist_ok=True)

    h5_path = os.path.join(bdir, f"{a.batch}_copyNumbersSegmented_{a.bin}kb.h5")
    grid_rows = []
    with h5py.File(h5_path, "w") as h5:
        for c in order:
            b = bins[bins.cell_id == c]
            g = h5.create_group(code[c])
            bg = g.create_group("bins")
            bg["chromosome"] = b.chromosome.astype("S2").values
            bg["start"] = b.start.values.astype(np.int64)
            bg["end"] = b.end.values.astype(np.int64)
            bg["use"] = b["use"].astype(bool).values
            bg["copynumber"] = b.copynumber.values.astype(np.float64)
            bg["segmented"] = b.segmented.values.astype(np.float64)
            sg = g.create_group("samples")
            sg["total_reads"] = np.array([int(cells.loc[c, "total_reads"])], dtype=np.int64)
            used = b["use"].astype(bool).values
            for p, q, clon, pw, loh in grid_for_cell(
                    b.copynumber.values[used], b.segmented.values[used],
                    cells.loc[c, "total_reads"], int(cells.loc[c, "used_bins"]), a.bin):
                grid_rows.append((code[c], p, q, clon, pw, loh))

    grid = pd.DataFrame(grid_rows, columns=["sample", "ploidy", "purity", "clonality", "powered", "LOHok"])
    grid.to_csv(os.path.join(bdir, f"{a.batch}_gridSearch.csv"), index=False)

    with open(os.path.join(bdir, "samples.txt"), "w") as f:
        f.writelines(f"[ ] {code[c]}\n" for c in order)
    with open(os.path.join(a.key_dir, f"{a.batch}_key.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "cell_id", "bin_kb"])
        w.writerows((code[c], c, a.bin) for c in order)

    if a.qc:
        qc = pd.read_csv(a.qc, sep="\t").set_index("cell_id")
        info = qc.loc[order, [x for x in QC_COLS if x in qc.columns]].copy()
        info.index = [code[c] for c in order]
        info.index.name = "sample"
        info.to_csv(os.path.join(bdir, f"{a.batch}_cellinfo.tsv"), sep="\t")

    # run_curation.sh lists batches from <out-dir>/batches.tsv (label <TAB> flag)
    manifest = os.path.join(a.out_dir, "batches.tsv")
    listed = open(manifest).read().split() if os.path.exists(manifest) else []
    if a.batch not in listed:
        with open(manifest, "a") as f:
            f.write(f"{a.batch}\t\n")

    print(f"{len(ids)} cells -> {bdir} ({'not ' if a.no_blind else ''}blinded); key -> {a.key_dir}")


if __name__ == "__main__":
    main()

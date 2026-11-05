#!/usr/bin/env python3
"""Apply the recalibrated depth_ok rule of cn_qc_metrics.R to existing cn_qc_<bin>kb.tsv tables.

depth_ok       mean corrected reads per used bin / 2 >= 15, i.e. nrpbpcc x ploidy / 2 >= 15
               (nrpbpcc = mean corrected reads per used bin / fitted ploidy, so the product is
               the per-bin mean and the threshold is the diploid one, independent of the fit)
depth_ok_reads former rule, total reads >= 15 x 2 x used_bins / 0.91, kept for reference

Equivalent to re-running cn_qc_metrics.R (same inputs) without re-reading the RDS objects.
Idempotent: the former rule is stored once in depth_ok_reads and never overwritten.

Usage: python3 scripts/03_cn_qc/recalibrate_depth.py --qc results/pilot_cn/qc
"""
import argparse

import pandas as pd

BINS = [200, 500, 1000]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    a = ap.parse_args()
    for b in BINS:
        f = f"{a.qc}/cn_qc_{b}kb.tsv"
        q = pd.read_csv(f, sep="\t")
        if "depth_ok_reads" not in q.columns:
            q.insert(q.columns.get_loc("depth_ok") + 1, "depth_ok_reads", q.depth_ok)
        q["depth_ok"] = (q.nrpbpcc * q.ploidy / 2 >= 15).astype(int)
        q.to_csv(f, sep="\t", index=False)
        print(f"{b} kb: depth_ok {q.depth_ok.sum()} (former rule {q.depth_ok_reads.sum()}) of {len(q)} cells")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Genome coverage of each CFS catalogue and rough power of the breakend overlap test.

For each catalogue x cell type (plus Le Tallec "all" and "epithelial" =
HCT116/LoVo/LS174T/MCF10A/CAL51), on autosomes:
  cov_<pad>  % of the hg19 autosomal genome (2,881 Mb) covered by sites padded by <pad> on
             each side, after merging overlaps. Padding approximates the breakend interval:
             half the median interval width (0.5 Mb at 500 kb bins, 1 Mb at 1 Mb bins).
  pw<f>_<bin> power to detect fold f (1.5 or 2) with a one-sided binomial test at alpha 0.05,
             for n independent events: n = 380 at 500 kb and 230 at 1 Mb (pilot, collapsed).
Limitations: assumes independent events placed uniformly on the genome (no RT
stratification, no usable-bin mask), a single test (no multiple-testing correction) and the
pilot's uncalibrated event counts. It is an order of magnitude, not a design calculation.

Usage: python3 scripts/05_features/cfs_power.py --cfs data/features/cfs/cfs_catalogs_hg19.tsv \
    --out data/features/cfs/cfs_power.tsv
"""
import argparse

import pandas as pd
from scipy.stats import binom

GENOME = 2881e6
EPITHELIAL = ["HCT116", "LoVo", "LS174T", "MCF10A", "CAL51"]
SCENARIOS = [("500kb", 0.5e6, 380), ("1Mb", 1e6, 230)]


def coverage(g, pad):
    tot = 0
    for _, x in g.groupby("chrom"):
        iv = sorted(zip(x.start - pad, x.end + pad))
        s, e = iv[0]
        for a, b in iv[1:]:
            if a <= e:
                e = max(e, b)
            else:
                tot += e - s
                s, e = a, b
        tot += e - s
    return tot / GENOME


def power(n, p, fold):
    k = binom.ppf(0.95, n, p)
    return 1 - binom.cdf(k, n, min(1.0, p * fold))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = pd.read_csv(a.cfs, sep="\t")
    d = d[d.chrom != "chrX"]
    lt = d[d.catalog == "letallec2013"]
    groups = list(d.groupby(["catalog", "cell_type"])) + [
        (("letallec2013", "all"), lt), (("letallec2013", "epithelial"), lt[lt.cell_type.isin(EPITHELIAL)])]
    rows = []
    for (cat, ct), g in groups:
        r = dict(catalog=cat, cell_type=ct, n_sites=len(g), cov_0=round(100 * coverage(g, 0), 1))
        for lab, pad, n in SCENARIOS:
            p = coverage(g, pad)
            r[f"cov_{lab}"] = round(100 * p, 1)
            for f in (1.5, 2):
                r[f"pw{f}_{lab}"] = round(power(n, p, f), 2)
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

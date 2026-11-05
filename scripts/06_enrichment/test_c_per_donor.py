#!/usr/bin/env python3
"""Test C: breakend-level observed / expected overlap per donor, and heterogeneity between donors.

Unit: each interstitial breakend of each cell (not collapsed: a clone contributes one
breakend per cell, so this test weights events by the number of cells carrying them).
For each breakend, o = its [lo, hi] overlaps the feature set; p = its exact null overlap
probability (same chromosome, RT quintile and width; grid-aligned placements; outside
pericentromeric blocks), as in test A.
Per donor: O_d = sum o, E_d = sum p, ratio = O_d / E_d, 95% CI from 2,000 bootstrap resamples
of cells within the donor. Pooled ratio = sum O / sum E with a cell bootstrap stratified by
donor. Heterogeneity: X2 = sum_d (O_d - R E_d)^2 / (R E_d), R = pooled ratio, df = donors - 1
(Poisson approximation; donors with E_d < 1 are left out of X2).
When O = 0 the bootstrap interval is degenerate [0, 0]; the upper limit is then the exact
one-sided 95% Poisson bound -ln(0.05) / E = 3.0 / E and ci = "one-sided".
Caveat: breakend-level units weight clones by their number of cells (one large clone, such as
P08's gain clone, dominates its donor), and X2 ignores that dependence.
Change from the plan (mixed logistic model with donor random effect): with 5 donors the
random-effect variance is not estimable, so per-donor ratios and a heterogeneity test are
reported instead.

Usage: python3 scripts/06_enrichment/test_c_per_donor.py --events results/enrichment/pilot/main \
    --cfs data/features/cfs/cfs_catalogs_hg19.tsv --genes data/features/genes --out results/enrichment/pilot/main
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(__file__))
from test_a_permutation import Null, feature_sets, overlaps  # noqa: E402

BINS = [200, 500, 1000]
SETS = ["macheret2020:all", "macheret2020:U2OS", "macheret2020:HeLa", "macheret2020:HS68", "ji2020:all",
        "ji2020:HBE", "humcfs:all", "letallec2013:all", "letallec2013:epithelial", "genes_1mb",
        "brison2019:all", "sarni2020:all", "wilson2015:all", "simpson2021:all", "consensusA:core", "consensusA:reproducible", "consensusA:all", "breast_predicted:all", "A_plus_C:all", "union_all:all"]
N_BOOT = 2000
UPPER0 = -np.log(0.05)


def boot_ratio(df, rng):
    # resample cells within each donor; ratio = sum o / sum p over the resampled cells
    acc_o, acc_p = np.zeros(N_BOOT), np.zeros(N_BOOT)
    for _, g in df.groupby("donor"):
        cells = g.groupby("cell_id").agg(o=("o", "sum"), p=("p", "sum"))
        idx = rng.integers(0, len(cells), (N_BOOT, len(cells)))
        acc_o += cells.o.values[idx].sum(1)
        acc_p += cells.p.values[idx].sum(1)
    r = acc_o / np.maximum(acc_p, 1e-12)
    return np.quantile(r, .025), np.quantile(r, .975)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    fs = feature_sets(a.cfs, a.genes)
    rows = []
    for b in BINS:
        bins = pd.read_csv(f"{a.events}/usable_bins_rt_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        null = Null(bins)
        edges = null.edges
        be = pd.read_csv(f"{a.events}/breakends_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        be = be[(be.type == "interstitial") & be.rt.notna()].copy()
        be["rtq"] = np.searchsorted(edges, be.rt) + 1
        be["main_def"] = be.stratum.isin(["E2", "E3", "E4"])
        for name in SETS:
            feat = fs[name]
            o = np.zeros(len(be), bool); p = np.full(len(be), np.nan)
            cache = {}
            for i, r in enumerate(be.itertuples()):
                o[i] = overlaps(feat, r.chromosome, r.lo, r.hi)[0]
                key = (r.chromosome, r.rtq, int(r.hi - r.lo))
                if key not in cache:
                    x, y = null.candidates(r.chromosome, int(r.rtq), int(r.hi - r.lo))
                    cache[key] = overlaps(feat, r.chromosome, x, y).mean() if len(x) else np.nan
                p[i] = cache[key]
            b2 = be.assign(o=o, p=p).dropna(subset=["p"])
            for d in ("loss", "gain", "all"):
                x = b2[b2.main_def] if d == "all" else b2[b2.main_def & (b2.direction == d)]
                if x.empty:
                    continue
                O, E = x.o.sum(), x.p.sum()
                R = O / E if E > 0 else np.nan
                lo, hi = boot_ratio(x, rng)
                ci = "bootstrap"
                if O == 0 and E > 0:
                    lo, hi, ci = 0.0, UPPER0 / E, "one-sided"
                dd = x.groupby("donor").agg(O=("o", "sum"), E=("p", "sum"), cells=("cell_id", "nunique"), n=("o", "size"))
                k = dd.E >= 1
                X2 = (((dd.O - R * dd.E) ** 2) / (R * dd.E))[k].sum() if R > 0 else np.nan
                p_het = chi2.sf(X2, k.sum() - 1) if k.sum() > 1 else np.nan
                rows.append(dict(bin_kb=b, feature_set=name, direction=d, donor="pooled", breakends=len(x),
                                 cells=x.cell_id.nunique(), observed=int(O), expected=round(E, 2), ratio=round(R, 3),
                                 ratio_lo=round(lo, 3), ratio_hi=round(hi, 3), ci=ci, het_X2=round(X2, 3), het_p=p_het))
                for dn, r in dd.iterrows():
                    l2, h2 = boot_ratio(x[x.donor == dn], rng)
                    c2 = "bootstrap"
                    if r.O == 0 and r.E > 0:
                        l2, h2, c2 = 0.0, UPPER0 / r.E, "one-sided"
                    rows.append(dict(bin_kb=b, feature_set=name, direction=d, donor=dn, breakends=int(r.n),
                                     cells=int(r.cells), observed=int(r.O), expected=round(r.E, 2),
                                     ratio=round(r.O / r.E, 3) if r.E > 0 else np.nan,
                                     ratio_lo=round(l2, 3), ratio_hi=round(h2, 3), ci=c2))
        print(f"{b} kb done")
    t = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    t.to_csv(f"{a.out}/test_c.tsv", sep="\t", index=False)
    s = t[(t.direction == "loss") & t.feature_set.isin(["macheret2020:all", "macheret2020:HS68"])]
    print(s[["bin_kb", "feature_set", "donor", "breakends", "observed", "expected", "ratio", "ratio_lo", "ratio_hi", "het_p"]].to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Projected power of test A when the analysis is scaled from the 5-donor pilot to all 49 donors.

For each set x bin x direction (main stratum, interstitial):
  E_full = E_pilot x k, with k = cells in the full dataset / cells in the pilot (91,374 /
           10,102 = 9.05): assumes the same events per cell and the same mix of donors
  power(f, alpha) = P(O >= c), O ~ Poisson(f x E_full), c = smallest count with
           P(Poisson(E_full) >= c) <= alpha (one-sided exact Poisson test, a close approximation
           of the Poisson-binomial when every p_i is small)
  min_fold_80 = smallest fold with power >= 0.8, at alpha = 0.05 (one test) and at
           alpha = 0.05 / m (m = number of main-stratum CFS rows in the pilot, a Bonferroni-type
           stand-in for the global correction)
Limitations: ignores the extra variance from event clustering (the site-shuffle p-values of
the pilot were ~1-6x larger than the exact ones), assumes the pilot's event rate per cell, and
treats the fold as constant across donors.

Usage: python3 scripts/06_enrichment/power_scaleup.py --res results/enrichment/pilot --out results/enrichment/pilot/power_scaleup.tsv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import poisson

CELLS_FULL, CELLS_PILOT = 91374, 10102
CONTROLS = ["genes_300kb", "genes_1mb"]


def crit(E, alpha):
    c = int(poisson.ppf(1 - alpha, E)) + 1
    while c > 0 and poisson.sf(c - 2, E) <= alpha:
        c -= 1
    while poisson.sf(c - 1, E) > alpha:
        c += 1
    return c


def power(f, E, alpha):
    return poisson.sf(crit(E, alpha) - 1, f * E)


def min_fold(E, alpha, target=0.8):
    for f in np.arange(1.0, 10.01, 0.01):
        if power(f, E, alpha) >= target:
            return round(float(f), 2)
    return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    k = CELLS_FULL / CELLS_PILOT
    t = pd.read_csv(f"{a.res}/main/test_a.tsv", sep="\t")
    t = t[(t.stratum == "main") & (t.expected > 0)]
    m = int((~t.feature_set.isin(CONTROLS)).sum())
    rows = []
    for r in t.itertuples():
        E = r.expected * k
        rows.append(dict(bin_kb=r.bin_kb, feature_set=r.feature_set, direction=r.direction,
                         expected_pilot=r.expected, expected_full=round(E, 1), fold_pilot=r.fold,
                         power_fold1_5=round(power(1.5, E, .05), 3), power_fold2=round(power(2.0, E, .05), 3),
                         min_fold_80=min_fold(E, .05), min_fold_80_global=min_fold(E, .05 / m)))
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"scale factor {k:.2f}, global m = {m}")
    print(out[(out.direction == "loss") & (out.bin_kb == 500)].to_string(index=False))


if __name__ == "__main__":
    main()

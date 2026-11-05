#!/usr/bin/env python3
"""Robustness tables for test A (main stratum = E2-E4, interstitial events).

From <res>/main/test_a_event_probs.tsv (o_i, p_i per event and set, written by
test_a_permutation.py), without re-running the null:
  robustness_loco.tsv      leave-one-chromosome-out: O, E, fold and exact Poisson-binomial
                           p_enrich after removing each autosome in turn
  robustness_terminal.tsv  events split by `terminal` (the altered flank runs to the chromosome
                           end) versus internal (flanked on both sides by retained segments)
From the test_a*.tsv of every run (main stratum rows):
  robustness_variants.tsv  one row per variant x bin x set x direction: main, nocuts,
                           samecells, cells500 (500 kb main cells called at every bin),
                           refined (500 kb breakends refined to 200 kb), linkage 1 / 5 Mb,
                           collapse on E2-E4 breakends only, and leave-one-donor-out

Usage: python3 scripts/06_enrichment/robustness.py --res results/enrichment/pilot
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from test_a_permutation import pb_tails  # noqa: E402

VARIANTS = [("main", "main/test_a.tsv"), ("nocuts", "nocuts/test_a.tsv"), ("samecells", "samecells/test_a.tsv"),
            ("cells500", "cells500/test_a.tsv"), ("refined", "refined/test_a.tsv"),
            ("link 1 Mb", "main_link1mb/test_a.tsv"), ("link 5 Mb", "main_link5mb/test_a.tsv"),
            ("collapse E2-E4", "main_collapse_main/test_a.tsv")] + \
           [(f"without {d}", f"main/test_a_without_{d}.tsv") for d in ("P07", "P08", "P18", "P20", "P36")]
DIRS = ("loss", "gain", "all")


def summarise(g):
    g = g[np.isfinite(g.p)]
    O, E = int(g.o.sum()), float(g.p.sum())
    return dict(n_events=len(g), observed=O, expected=round(E, 3), fold=round(O / E, 3) if E > 0 else np.nan,
                p_enrich=pb_tails(g.p.values, O)[0] if len(g) else np.nan)


def by_direction(ev):
    for d in DIRS:
        yield d, (ev if d == "all" else ev[ev.direction == d])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    a = ap.parse_args()
    pr = pd.read_csv(f"{a.res}/main/test_a_event_probs.tsv", sep="\t", dtype={"chromosome": str})
    pr = pr[(pr.type == "interstitial") & pr.main_def]

    loco, term = [], []
    for (b, name), g in pr.groupby(["bin_kb", "feature_set"]):
        for d, x in by_direction(g):
            for c in sorted(x.chromosome.unique(), key=int):
                loco.append(dict(bin_kb=b, feature_set=name, direction=d, chromosome_out=c, **summarise(x[x.chromosome != c])))
            for lab, y in (("terminal", x[x.terminal]), ("internal", x[~x.terminal])):
                term.append(dict(bin_kb=b, feature_set=name, direction=d, events=lab, **summarise(y)))
    pd.DataFrame(loco).to_csv(f"{a.res}/robustness_loco.tsv", sep="\t", index=False)
    pd.DataFrame(term).to_csv(f"{a.res}/robustness_terminal.tsv", sep="\t", index=False)

    rows = []
    for lab, f in VARIANTS:
        p = f"{a.res}/{f}"
        if not os.path.exists(p):
            print("missing:", p)
            continue
        t = pd.read_csv(p, sep="\t")
        t = t[t.stratum == "main"]
        keep = [c for c in ("bin_kb", "feature_set", "direction", "n_events", "observed", "expected", "fold",
                            "fold_lo", "fold_hi", "p_enrich", "p_perm") if c in t]
        rows.append(t[keep].assign(variant=lab))
    v = pd.concat(rows, ignore_index=True)
    v.to_csv(f"{a.res}/robustness_variants.tsv", sep="\t", index=False)
    x = v[(v.feature_set == "macheret2020:HS68") & (v.direction == "loss")]
    print(x.pivot_table(index="variant", columns="bin_kb", values="fold").round(2).to_string())
    lc = pd.DataFrame(loco)
    y = lc[(lc.feature_set == "macheret2020:HS68") & (lc.direction == "loss") & (lc.bin_kb == 500)]
    print("HS68 500 kb loss, leave-one-chromosome-out fold range:", y.fold.min(), "-", y.fold.max(),
          "| worst:", y.loc[y.fold.idxmin(), "chromosome_out"])


if __name__ == "__main__":
    main()

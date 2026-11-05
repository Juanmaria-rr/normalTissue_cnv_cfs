#!/usr/bin/env python3
"""Test A: positional enrichment of interstitial breakend events in genomic feature sets.

Null: each event is re-placed on its own chromosome, centred on a usable bin of the same
replication-timing quintile, keeping its interval width, never overlapping a pericentromeric
block (+/- 1 Mb). Because placements are uniform over that candidate set, each event's null
overlap probability p_i is computed exactly by enumerating all its candidate placements.
  observed O = number of events whose [lo, hi] overlaps a feature
  expected E = sum p_i;  fold = O / E
  p (enrichment / depletion) from 20,000 simulated sums of independent Bernoulli(p_i),
  equivalent to 20,000 joint permutations of all events; fold 95% CI from 2,000 bootstrap
  resamples of events (o_i, p_i)
  q = BH within each feature set (family) across bins x directions x strata
Subsets: interstitial events; direction loss / gain / all; stratum E1, E2, E3, E4, main
(E2-E4 = the paper's 1.3 / 0.7).
Feature sets: each CFS catalogue x cell type and unions (macheret2020, ji2020, letallec2013
all and epithelial, humcfs, brison2019, sarni2020, wilson2015, simpson2021), and the derived
catalogues of cfs_derived_hg19.tsv when present (consensusA core / reproducible / all,
breast_predicted, A_plus_C, union_all); controls: genes > 300 kb and > 1 Mb (positive), and
RT-matched shuffles of selected sets (negative: sites moved to random usable bins of the
same chromosome and RT quintile, width kept; fold should be ~1).

Inputs: <events>/events_<bin>kb.tsv and usable_bins_rt_<bin>kb.tsv (build_event_table.py),
data/features/cfs/cfs_catalogs_hg19.tsv, data/features/genes/large_genes_{300kb,1mb}.tsv.
Options: --shuffle unions|all|none (which sets get shuffle controls), --exclude-donor
(leave-one-donor-out), --tag (output suffix).
Outputs: <out>/test_a<tag>.tsv, <out>/test_a_shuffle_controls<tag>.tsv.

Usage: python3 scripts/06_enrichment/test_a_permutation.py --events results/enrichment/pilot/main \
    --cfs data/features/cfs/cfs_catalogs_hg19.tsv --genes data/features/genes --out results/enrichment/pilot/main
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from genome_annot import PERICEN  # noqa: E402

BINS = [200, 500, 1000]
EPITHELIAL = ["HCT116", "LoVo", "LS174T", "MCF10A", "CAL51"]
N_SIM, N_BOOT, N_SHUF = 20000, 2000, 20
SHUFFLE_SETS = ["macheret2020:all", "ji2020:all", "humcfs:all", "letallec2013:all", "genes_1mb",
                "consensusA:all", "breast_predicted:all", "A_plus_C:all", "union_all:all", "simpson2021:all"]


def merge(df):
    out = {}
    for c, g in df.groupby("chromosome"):
        iv = sorted(zip(g.start.astype(np.int64), g.end.astype(np.int64)))
        S, E = [], []
        for s, e in iv:
            if S and s <= E[-1]:
                E[-1] = max(E[-1], e)
            else:
                S.append(s); E.append(e)
        out[str(c)] = (np.array(S), np.array(E))
    return out


def overlaps(feat, chrom, a, b):
    a, b = np.atleast_1d(a), np.atleast_1d(b)
    if chrom not in feat:
        return np.zeros(len(a), bool)
    S, E = feat[chrom]
    i = np.searchsorted(E, a, side="left")
    ok = i < len(S)
    res = np.zeros(len(a), bool)
    res[ok] = S[i[ok]] <= b[ok]
    return res


def feature_sets(cfs_file, genes_dir):
    c = pd.read_csv(cfs_file, sep="\t", dtype={"chrom": str})
    derived = os.path.join(os.path.dirname(cfs_file), "cfs_derived_hg19.tsv")
    if os.path.exists(derived):  # consensusA, breast_predicted, A_plus_C, union_all (build_cfs_consensus.py)
        c = pd.concat([c, pd.read_csv(derived, sep="\t", dtype={"chrom": str})], ignore_index=True)
    c = c.assign(chromosome=c.chrom.str.replace("chr", ""))
    c = c[c.chromosome.isin([str(i) for i in range(1, 23)])]
    fs = {}
    for (cat, ct), g in c.groupby(["catalog", "cell_type"]):
        fs[f"{cat}:{ct}"] = g
    for cat, g in c.groupby("catalog"):
        fs[f"{cat}:all"] = g
    lt = c[c.catalog == "letallec2013"]
    fs["letallec2013:epithelial"] = lt[lt.cell_type.isin(EPITHELIAL)]
    for lab in ("300kb", "1mb"):
        g = pd.read_csv(f"{genes_dir}/large_genes_{lab}.tsv", sep="\t", dtype={"chromosome": str})
        fs[f"genes_{lab}"] = g
    return {k: merge(v[["chromosome", "start", "end"]]) for k, v in fs.items()}


class Null:
    """Candidate placements per (chromosome, RT quintile): usable-bin midpoints."""

    def __init__(self, bins):
        q = np.nanquantile(bins.rt, [.2, .4, .6, .8])
        bins = bins.assign(rtq=np.searchsorted(q, bins.rt) + 1, mid=(bins.start + bins.end) // 2)
        self.mids = {k: g.mid.values.astype(np.int64) for k, g in bins.groupby(["chromosome", "rtq"])}
        self.bins = bins

    def candidates(self, chrom, rtq, width):
        m = self.mids.get((chrom, rtq), np.array([], np.int64))
        a, b = m - width // 2, m + width // 2
        bl, bh = PERICEN[int(chrom)]
        keep = ~((a / 1e6 <= bh + 1) & (b / 1e6 >= bl - 1))
        return a[keep], b[keep]


def event_probs(ev, null, feat):
    o, p = np.zeros(len(ev), bool), np.full(len(ev), np.nan)
    for i, r in enumerate(ev.itertuples()):
        o[i] = overlaps(feat, r.chromosome, r.lo, r.hi)[0]
        a, b = null.candidates(r.chromosome, int(r.rt_quintile), int(r.hi - r.lo))
        if len(a):
            p[i] = overlaps(feat, r.chromosome, a, b).mean()
    return o, p


def stats(o, p, rng):
    k = np.isfinite(p)
    o, p = o[k], p[k]
    n = len(o)
    if n == 0:
        return dict(n_events=0)
    O, E = int(o.sum()), float(p.sum())
    sims = (rng.random((N_SIM, n)) < p).sum(1)
    bi = rng.integers(0, n, (N_BOOT, n))
    bf = o[bi].sum(1) / np.maximum(p[bi].sum(1), 1e-12)
    return dict(n_events=n, observed=O, expected=round(E, 2), fold=round(O / E, 3) if E > 0 else np.nan,
                fold_lo=round(np.quantile(bf, .025), 3), fold_hi=round(np.quantile(bf, .975), 3),
                p_enrich=(1 + (sims >= O).sum()) / (N_SIM + 1), p_deplete=(1 + (sims <= O).sum()) / (N_SIM + 1))


def subsets(ev):
    ev = ev[ev.type == "interstitial"]
    for d in ("loss", "gain", "all"):
        x = ev if d == "all" else ev[ev.direction == d]
        for s in ("E1", "E2", "E3", "E4", "main"):
            y = x[x.main_def] if s == "main" else x[x.stratum == s]
            yield d, s, y


def shuffle_feature(feat, null, rng):
    rows = []
    for c, (S, E) in feat.items():
        for s, e in zip(S, E):
            b = null.bins[(null.bins.chromosome == c) & (null.bins.end >= s) & (null.bins.start <= e)]
            rtq = int(b.rtq.mode().iloc[0]) if len(b) else int(rng.integers(1, 6))
            m = null.mids.get((c, rtq))
            if m is None or not len(m):
                continue
            mid = rng.choice(m)
            rows.append((c, mid - (e - s) // 2, mid + (e - s) // 2))
    return merge(pd.DataFrame(rows, columns=["chromosome", "start", "end"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--shuffle", default="unions", choices=["unions", "all", "none"],
                    help="feature sets with RT-matched shuffle controls")
    ap.add_argument("--exclude-donor", default=None, help="leave-one-donor-out sensitivity")
    ap.add_argument("--tag", default="", help="suffix of the output files")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    fs = feature_sets(a.cfs, a.genes)
    shuffle_sets = {"unions": SHUFFLE_SETS, "all": list(fs), "none": []}[a.shuffle]

    rows, ctrl = [], []
    for b in BINS:
        ev = pd.read_csv(f"{a.events}/events_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        ev = ev[ev.rt_quintile.notna()]
        if a.exclude_donor:
            ev = ev[ev.donor != a.exclude_donor]
        null = Null(pd.read_csv(f"{a.events}/usable_bins_rt_{b}kb.tsv", sep="\t", dtype={"chromosome": str}))
        for name, feat in fs.items():
            o, p = event_probs(ev, null, feat)
            ev2 = ev.assign(o=o, p=p)
            for d, s, y in subsets(ev2):
                r = stats(y.o.values, y.p.values, rng)
                rows.append(dict(bin_kb=b, feature_set=name, direction=d, stratum=s, **r))
        # negative controls: RT-matched shuffles, main stratum, each direction (loss, gain, all)
        for name in shuffle_sets:
            for k in range(N_SHUF):
                sf = shuffle_feature(fs[name], null, rng)
                o, p = event_probs(ev, null, sf)
                ev2 = ev.assign(o=o, p=p)
                for d, s, y in subsets(ev2):
                    if s != "main":
                        continue
                    r = stats(y.o.values, y.p.values, rng)
                    ctrl.append(dict(bin_kb=b, feature_set=name, shuffle=k, direction=d, **{x: r.get(x) for x in ("n_events", "observed", "expected", "fold")}))
        print(f"{b} kb done")

    t = pd.DataFrame(rows)
    t["q_enrich"] = np.nan
    for name, g in t.groupby("feature_set"):
        k = g.p_enrich.notna()
        pv = g.loc[k, "p_enrich"].values
        order = np.argsort(pv); m = len(pv)
        qv = np.minimum.accumulate((pv[order] * m / np.arange(1, m + 1))[::-1])[::-1]
        qq = np.empty(m); qq[order] = np.minimum(qv, 1)
        t.loc[g.index[k], "q_enrich"] = qq
    os.makedirs(a.out, exist_ok=True)
    t.to_csv(f"{a.out}/test_a{a.tag}.tsv", sep="\t", index=False)
    c = pd.DataFrame(ctrl)
    if c.empty:
        return
    c.to_csv(f"{a.out}/test_a_shuffle_controls{a.tag}.tsv", sep="\t", index=False)
    print(c.groupby(["bin_kb", "feature_set", "direction"]).fold.describe(percentiles=[.025, .5, .975])[["count", "2.5%", "50%", "97.5%"]].round(2).to_string())


if __name__ == "__main__":
    main()

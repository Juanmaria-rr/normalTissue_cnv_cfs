#!/usr/bin/env python3
"""Test A: positional enrichment of interstitial breakend events in genomic feature sets.

Null (grid-aligned, revised 2026-11-05 after the independent review): each event is re-placed
on its own chromosome with its interval starting at a usable bin start, keeping its width,
never overlapping a pericentromeric block (+/- 1 Mb), and only at placements whose
replication timing falls in the event's RT quintile. The RT of a placement is computed
exactly as for events: mean RT of the usable bins overlapping the interval, quintile edges
from the genome-wide usable-bin RT of that bin size. Placements are uniform over that
candidate set, so each event's null overlap probability p_i is exact (all candidates are
enumerated).
  observed O = number of events whose [lo, hi] overlaps a feature
  expected E = sum p_i;  fold = O / E
  p_enrich / p_deplete  exact Poisson-binomial tails P(S >= O) / P(S <= O), S = sum of
                        independent Bernoulli(p_i) (no Monte-Carlo error)
  fold_lo / fold_hi     2.5 / 97.5% of 2,000 bootstrap resamples of events (o_i, p_i)
  fold_lo_donor / _hi   same, resampling donors (cluster bootstrap; with 5 donors a rough
                        guide to between-donor variability, not a calibrated interval)
  p_perm                main stratum only: site-shuffle permutation p. The feature set is
                        shuffled --perm times (default 5,000; each site moved to a random usable bin of the
                        same chromosome and RT quintile, width kept) and the whole test is
                        recomputed; p_perm = (1 + #{fold_shuffled >= fold_observed}) / (N + 1).
                        Unlike the exact p, it keeps the spatial clustering of events and of
                        sites, so it is the reference p-value.
  q_family              BH within each feature set across bins x directions x strata
  q_global / q_perm_global  BH over all main-stratum rows of every CFS set (bins x directions
                        x sets; positive controls genes_* excluded), on p_enrich / p_perm
Subsets: interstitial events; direction loss / gain / all; stratum E1, E2, E3, E4, main
(E2-E4 = the paper's 1.3 / 0.7).
Feature sets (deduplicated): <catalogue>:all for every catalogue; per cell type only where
cell types differ and are of interest (macheret2020 HS68 / HeLa / U2OS, ji2020 HBE / HeLa /
U2OS, letallec2013 MCF10A and epithelial = HCT116 / LoVo / LS174T / MCF10A / CAL51,
consensusA core / reproducible); controls genes > 300 kb and > 1 Mb (positive).

Inputs: <events>/events_<bin>kb.tsv and usable_bins_rt_<bin>kb.tsv (build_event_table.py),
data/features/cfs/cfs_catalogs_hg19.tsv (+ cfs_derived_hg19.tsv), data/features/genes/.
Options: --perm N (0 = no permutation), --bins, --exclude-donor (leave-one-donor-out),
--tag (output suffix), --workers.
Outputs: <out>/test_a<tag>.tsv, <out>/test_a_shuffle_controls<tag>.tsv (fold of every
shuffle, main stratum), <out>/test_a_event_probs<tag>.tsv (o_i, p_i per event and set, for
robustness.py).

Usage: python3 scripts/06_enrichment/test_a_permutation.py --events results/enrichment/pilot/main \
    --cfs data/features/cfs/cfs_catalogs_hg19.tsv --genes data/features/genes --out results/enrichment/pilot/main
"""
import argparse
import os
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from genome_annot import PERICEN  # noqa: E402

BINS = [200, 500, 1000]
EPITHELIAL = ["HCT116", "LoVo", "LS174T", "MCF10A", "CAL51"]
PER_CELL_TYPE = {"macheret2020": None, "ji2020": None, "consensusA": None, "letallec2013": ["MCF10A"]}
CONTROLS = ["genes_300kb", "genes_1mb"]
N_BOOT = 2000
DIRS = ("loss", "gain", "all")
STRATA = ("E1", "E2", "E3", "E4", "main")


def merge(df):
    out = {}
    for c, g in df.groupby("chromosome"):
        out[str(c)] = merge_arrays(g.start.values.astype(np.int64), g.end.values.astype(np.int64))
    return out


def merge_arrays(s, e):
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    ce = np.maximum.accumulate(e)
    new = np.r_[True, s[1:] > ce[:-1]]
    g = np.cumsum(new) - 1
    S = s[new]
    E = np.zeros(len(S), np.int64)
    np.maximum.at(E, g, e)
    return S, E


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
    for cat, g in c.groupby("catalog"):
        fs[f"{cat}:all"] = g
        if cat in PER_CELL_TYPE:
            for ct, h in g.groupby("cell_type"):
                if PER_CELL_TYPE[cat] is None or ct in PER_CELL_TYPE[cat]:
                    fs[f"{cat}:{ct}"] = h
    lt = c[c.catalog == "letallec2013"]
    fs["letallec2013:epithelial"] = lt[lt.cell_type.isin(EPITHELIAL)]
    for lab in ("300kb", "1mb"):
        g = pd.read_csv(f"{genes_dir}/large_genes_{lab}.tsv", sep="\t", dtype={"chromosome": str})
        fs[f"genes_{lab}"] = g
    return {k: merge(v[["chromosome", "start", "end"]]) for k, v in fs.items()}


class Null:
    """Grid-aligned candidate placements with RT computed as for events."""

    def __init__(self, bins):
        self.edges = np.nanquantile(bins.rt, [.2, .4, .6, .8])
        bins = bins.assign(rtq=np.searchsorted(self.edges, bins.rt) + 1, mid=(bins.start + bins.end) // 2)
        bins.loc[bins.rt.isna(), "rtq"] = 0
        self.bins = bins
        self.chrom = {}
        for c, g in bins.groupby("chromosome"):
            g = g.sort_values("start")
            r = g.rt.values
            fin = np.isfinite(r)
            self.chrom[c] = dict(start=g.start.values.astype(np.int64), end=g.end.values.astype(np.int64),
                                 cs=np.r_[0, np.cumsum(np.where(fin, r, 0))], cn=np.r_[0, np.cumsum(fin)],
                                 rtq=g.rtq.values.astype(int))
        # shuffle targets: usable-bin midpoints per (chromosome, bin RT quintile)
        self.mids = {k: g.mid.values.astype(np.int64) for k, g in bins[bins.rtq > 0].groupby(["chromosome", "rtq"])}
        self.cache = {}

    def interval_rtq(self, chrom, lo, hi):
        d = self.chrom[chrom]
        i = np.searchsorted(d["end"], lo, side="left")
        j = np.searchsorted(d["start"], hi, side="right")
        n = d["cn"][j] - d["cn"][i]
        with np.errstate(invalid="ignore", divide="ignore"):
            rt = (d["cs"][j] - d["cs"][i]) / n
        q = np.searchsorted(self.edges, rt) + 1
        q[~np.isfinite(rt)] = 0
        return q

    def candidates(self, chrom, rtq, width):
        key = (chrom, width)
        if key not in self.cache:
            if chrom not in self.chrom:
                self.cache[key] = (np.array([], np.int64),) * 2 + (np.array([], int),)
            else:
                a = self.chrom[chrom]["start"]
                b = a + width
                bl, bh = PERICEN[int(chrom)]
                keep = ~((a / 1e6 <= bh + 1) & (b / 1e6 >= bl - 1))
                a, b = a[keep], b[keep]
                self.cache[key] = (a, b, self.interval_rtq(chrom, a, b))
        a, b, q = self.cache[key]
        k = q == rtq
        return a[k], b[k]


def event_probs(ev, null, feat):
    o = np.zeros(len(ev), bool)
    p = np.full(len(ev), np.nan)
    chrom, lo, hi = ev.chromosome.values, ev.lo.values.astype(np.int64), ev.hi.values.astype(np.int64)
    rtq, width = ev.rt_quintile.values.astype(int), (ev.hi - ev.lo).values.astype(np.int64)
    for c in np.unique(chrom):
        k = np.where(chrom == c)[0]
        o[k] = overlaps(feat, c, lo[k], hi[k])
    pc = {}
    for i in range(len(ev)):
        key = (chrom[i], rtq[i], width[i])
        if key not in pc:
            a, b = null.candidates(*key)
            pc[key] = overlaps(feat, chrom[i], a, b).mean() if len(a) else np.nan
        p[i] = pc[key]
    return o, p


def pb_pmf(p):
    """Exact Poisson-binomial probability mass of S = sum of independent Bernoulli(p_i)."""
    pmf = np.zeros(len(p) + 1)
    pmf[0] = 1.0
    for k, x in enumerate(p, start=1):
        pmf[1:k + 1] = pmf[1:k + 1] * (1 - x) + pmf[:k] * x
        pmf[0] *= 1 - x
    return pmf


def pb_tails(p, O):
    """Exact Poisson-binomial P(S >= O), P(S <= O)."""
    pmf = pb_pmf(p)
    return float(min(1.0, pmf[O:].sum())), float(min(1.0, pmf[:O + 1].sum()))


def rng_for(seed, *key):
    return np.random.default_rng([seed, zlib.crc32("|".join(map(str, key)).encode())])


def stats(o, p, donors, rng):
    k = np.isfinite(p)
    o, p, donors = o[k], p[k], donors[k]
    n = len(o)
    if n == 0:
        return dict(n_events=0)
    O, E = int(o.sum()), float(p.sum())
    pe, pd_ = pb_tails(p, O)
    bi = rng.integers(0, n, (N_BOOT, n))
    bf = o[bi].sum(1) / np.maximum(p[bi].sum(1), 1e-12)
    ud = np.unique(donors)
    do = np.array([o[donors == d].sum() for d in ud])
    dp = np.array([p[donors == d].sum() for d in ud])
    di = rng.integers(0, len(ud), (N_BOOT, len(ud)))
    df = do[di].sum(1) / np.maximum(dp[di].sum(1), 1e-12)
    return dict(n_events=n, observed=O, expected=round(E, 3), fold=round(O / E, 3) if E > 0 else np.nan,
                fold_lo=round(np.quantile(bf, .025), 3), fold_hi=round(np.quantile(bf, .975), 3),
                fold_lo_donor=round(np.quantile(df, .025), 3), fold_hi_donor=round(np.quantile(df, .975), 3),
                p_enrich=pe, p_deplete=pd_)


def subsets(ev):
    ev = ev[ev.type == "interstitial"]
    for d in DIRS:
        x = ev if d == "all" else ev[ev.direction == d]
        for s in STRATA:
            y = x[x.main_def] if s == "main" else x[x.stratum == s]
            yield d, s, y


class Shuffler:
    """Moves each merged site to a random usable-bin midpoint of its chromosome and RT quintile."""

    def __init__(self, feat, null):
        rows = []
        for c, (S, E) in feat.items():
            if c not in null.chrom:
                continue
            d = null.chrom[c]
            i = np.searchsorted(d["end"], S, side="left")
            j = np.searchsorted(d["start"], E, side="right")
            for s, e, a, b in zip(S, E, i, j):
                q = d["rtq"][a:b]
                q = q[q > 0]
                rows.append((c, int(e - s), int(np.bincount(q).argmax()) if len(q) else 0))
        self.sites = pd.DataFrame(rows, columns=["chromosome", "width", "rtq"])
        self.groups = []
        for (c, q), g in self.sites.groupby(["chromosome", "rtq"]):
            m = null.mids.get((c, q))
            if q == 0:  # no RT information: any usable bin of the chromosome
                m = null.bins.mid.values[null.bins.chromosome.values == c].astype(np.int64)
            if m is not None and len(m):
                self.groups.append((c, g.width.values.astype(np.int64), m))

    def draw(self, rng):
        per = {}
        for c, w, m in self.groups:
            mid = m[rng.integers(0, len(m), len(w))]
            per.setdefault(c, []).append((mid - w // 2, mid + w // 2))
        return {c: merge_arrays(np.concatenate([x[0] for x in v]), np.concatenate([x[1] for x in v])) for c, v in per.items()}


# worker state (set once per process)
_W = {}


def _init(events_dir, bins, cfs, genes, exclude):
    _W["fs"] = feature_sets(cfs, genes)
    _W["ev"], _W["null"] = {}, {}
    for b in bins:
        ev = pd.read_csv(f"{events_dir}/events_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        ev = ev[ev.rt_quintile.notna()].reset_index(drop=True)
        if exclude:
            ev = ev[ev.donor != exclude].reset_index(drop=True)
        _W["ev"][b] = ev
        _W["null"][b] = Null(pd.read_csv(f"{events_dir}/usable_bins_rt_{b}kb.tsv", sep="\t", dtype={"chromosome": str}))


def _task(args):
    b, name, n_perm, seed = args
    ev, null, feat = _W["ev"][b], _W["null"][b], _W["fs"][name]
    o, p = event_probs(ev, null, feat)
    ev2 = ev.assign(o=o, p=p)
    rows, ctrl = [], []
    main_fold = {}
    for d, s, y in subsets(ev2):
        r = stats(y.o.values, y.p.values, y.donor.values, rng_for(seed, b, name, d, s))
        rows.append(dict(bin_kb=b, feature_set=name, direction=d, stratum=s, **r))
        if s == "main":
            main_fold[d] = r.get("fold", np.nan)
    if n_perm:
        sh = Shuffler(feat, null)
        rng = rng_for(seed, "perm", b, name)
        inter = ev[ev.type == "interstitial"]
        idx = {d: (inter.main_def if d == "all" else inter.main_def & (inter.direction == d)).values for d in DIRS}
        for k in range(n_perm):
            sf = sh.draw(rng)
            o2, p2 = event_probs(inter, null, sf)
            for d in DIRS:
                m = idx[d] & np.isfinite(p2)
                O, E = o2[m].sum(), p2[m].sum()
                ctrl.append((b, name, k, d, int(m.sum()), int(O), round(float(E), 3), O / E if E > 0 else np.nan))
        c = pd.DataFrame(ctrl, columns=["bin_kb", "feature_set", "shuffle", "direction", "n_events", "observed", "expected", "fold"])
        for r in rows:
            if r["stratum"] == "main" and r.get("n_events", 0):
                f = c.fold[c.direction == r["direction"]].values
                f = f[np.isfinite(f)]
                r["p_perm"] = (1 + (f >= r["fold"] - 1e-12).sum()) / (len(f) + 1)
                r["n_perm"] = len(f)
    probs = ev2[["event_id", "donor", "chromosome", "direction", "type", "stratum", "main_def", "terminal", "o", "p"]].assign(bin_kb=b, feature_set=name)
    return rows, ctrl, probs


def bh(pv):
    pv = np.asarray(pv, float)
    m = len(pv)
    if m == 0:
        return pv
    o = np.argsort(pv)
    q = np.minimum.accumulate((pv[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(q, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--perm", type=int, default=5000, help="site shuffles per set and bin (0 = none)")
    ap.add_argument("--bins", default="200,500,1000")
    ap.add_argument("--exclude-donor", default=None, help="leave-one-donor-out sensitivity")
    ap.add_argument("--tag", default="", help="suffix of the output files")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    a = ap.parse_args()
    bins = [int(x) for x in a.bins.split(",")]
    names = list(feature_sets(a.cfs, a.genes))
    tasks = [(b, n, a.perm, a.seed) for b in bins for n in names]
    rows, ctrl, probs = [], [], []
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(a.events, bins, a.cfs, a.genes, a.exclude_donor)) as ex:
        for i, (r, c, p) in enumerate(ex.map(_task, tasks)):
            rows += r; ctrl += c; probs.append(p)
            if (i + 1) % 20 == 0:
                print(f"{i + 1}/{len(tasks)} tasks", flush=True)

    t = pd.DataFrame(rows)
    t["q_family"] = np.nan
    for _, g in t.groupby("feature_set"):
        k = g.p_enrich.notna()
        t.loc[g.index[k], "q_family"] = bh(g.loc[k, "p_enrich"])
    glob = (t.stratum == "main") & ~t.feature_set.isin(CONTROLS) & t.p_enrich.notna()
    t.loc[glob, "q_global"] = bh(t.loc[glob, "p_enrich"])
    if "p_perm" in t:
        gp = glob & t.p_perm.notna()
        t.loc[gp, "q_perm_global"] = bh(t.loc[gp, "p_perm"])
    os.makedirs(a.out, exist_ok=True)
    t.to_csv(f"{a.out}/test_a{a.tag}.tsv", sep="\t", index=False)
    pd.concat(probs).to_csv(f"{a.out}/test_a_event_probs{a.tag}.tsv", sep="\t", index=False)
    if ctrl:
        c = pd.DataFrame(ctrl, columns=["bin_kb", "feature_set", "shuffle", "direction", "n_events", "observed", "expected", "fold"])
        c.to_csv(f"{a.out}/test_a_shuffle_controls{a.tag}.tsv", sep="\t", index=False)
    show = t[(t.stratum == "main") & (t.direction == "loss")]
    cols = [x for x in ("bin_kb", "feature_set", "observed", "expected", "fold", "p_enrich", "p_perm", "q_family", "q_global") if x in show]
    print(show.sort_values("p_enrich")[cols].head(12).to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Test B: per-window regression of interstitial event counts on feature overlap + covariates.

Windows: usable 1 Mb bins (autosomes, applyFilters use == TRUE), outside pericentromeric
blocks +/- 1 Mb; the same windows for events of every bin size.
Rows: donor x window. Response: number of interstitial events (subset by direction and
stratum) assigned to the window. Each event is assigned to the eligible window with the
largest overlap with its [lo, hi] interval (revised 2026-11-05: assignment by midpoint
dropped ~26% of events whose midpoint window was not eligible); events overlapping no
eligible window are excluded and counted in n_events_excluded.
Covariates per window (from the 50 kb annotation of sWGS_relativeCN, averaged):
  feature    1 if the window overlaps the feature set (binary)
  rt         mean replication timing (z-scored)
  gc         mean GC % (z-scored)
  map        mean mappability (z-scored)
  genes300   fraction of the window inside genes > 300 kb (refGene, merged)
  offset     log(fraction of the window's 50 kb windows with use == TRUE)
Models (statsmodels GLM, NB2 with alpha from the Cameron-Trivedi auxiliary regression on the
Poisson fit, Poisson when alpha <= 0, flagged in `family`). Standard error of `feature` = the
larger of the model-based and the cluster-robust (by chromosome, 22 clusters) one: clustering
stops windows of a chromosome being treated as independent, but with few events and few
clusters the sandwich estimator can be too small, so it is never allowed to shrink the SE
(column se_source says which one was used):
  crude                events ~ feature + donor + offset
  adjusted             events ~ feature + rt + gc + map + chromosome + donor + offset
  adjusted_genes300    adjusted + genes300. Large genes are a likely mechanism of fragility
                       (a mediator rather than a confounder), so this model is conservative
When no event (or every event) falls in feature windows the likelihood has no finite
maximum (separation): irr / p are left empty and `note` says so (previously p = 0).
q_family = BH within feature set x model across bins x directions x strata; q_global = BH
per model over the main-stratum rows of all CFS sets (genes_1mb excluded).
The IRR of `feature` is the event rate in windows overlapping the feature relative to the
rest. Limitations: midpoints only (interval width ignored); a binary overlap ignores how
much of the window is covered; donor and chromosome as fixed effects.

Usage: python3 scripts/06_enrichment/test_b_regression.py --events results/enrichment/pilot/main \
    --rt external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv \
    --cfs data/features/cfs/cfs_catalogs_hg19.tsv --genes data/features/genes --out results/enrichment/pilot/main
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(__file__))
from genome_annot import PERICEN  # noqa: E402
from test_a_permutation import bh, feature_sets, merge, overlaps  # noqa: E402

BINS = [200, 500, 1000]
SETS = ["macheret2020:all", "macheret2020:U2OS", "macheret2020:HeLa", "macheret2020:HS68", "ji2020:all",
        "ji2020:HBE", "humcfs:all", "letallec2013:all", "letallec2013:epithelial", "genes_1mb",
        "brison2019:all", "sarni2020:all", "wilson2015:all", "simpson2021:all", "consensusA:core", "consensusA:reproducible", "consensusA:all", "breast_predicted:all", "A_plus_C:all", "union_all:all"]
warnings.filterwarnings("ignore")


def windows(events_dir, rt_file, genes_dir):
    w = pd.read_csv(f"{events_dir}/usable_bins_rt_1000kb.tsv", sep="\t", dtype={"chromosome": str})
    keep = [not ((s / 1e6 <= PERICEN[int(c)][1] + 1) and (e / 1e6 >= PERICEN[int(c)][0] - 1))
            for c, s, e in zip(w.chromosome, w.start, w.end)]
    w = w[keep].reset_index(drop=True)
    a = pd.read_csv(rt_file, sep="\t", usecols=["seqnames", "start", "gc", "mappability", "use", "replicationTime"])
    a["seqnames"] = a.seqnames.astype(str)
    a["wstart"] = ((a.start - 1) // 1_000_000) * 1_000_000 + 1
    g = a.groupby(["seqnames", "wstart"]).agg(gc=("gc", "mean"), map=("mappability", "mean"),
                                              rt=("replicationTime", "mean"), usable=("use", "mean"))
    g = g.reindex(list(zip(w.chromosome, w.start)))
    for c in ("gc", "map", "rt", "usable"):
        w[c] = g[c].values
    w = w[(w.usable > 0) & w.rt.notna()].reset_index(drop=True)
    genes = merge(pd.read_csv(f"{genes_dir}/large_genes_300kb.tsv", sep="\t", dtype={"chromosome": str}))
    w["genes300"] = [covered(genes, c, s, e) for c, s, e in zip(w.chromosome, w.start, w.end)]
    for c in ("rt", "gc", "map"):
        w[c] = (w[c] - w[c].mean()) / w[c].std()
    w["log_usable"] = np.log(w.usable)
    return w


def covered(feat, chrom, s, e):
    if chrom not in feat:
        return 0.0
    S, E = feat[chrom]
    ov = np.clip(np.minimum(E, e) - np.maximum(S, s), 0, None)
    return float(ov.sum()) / (e - s + 1)


MODELS = (("crude", "events ~ feature + C(donor)"),
          ("adjusted", "events ~ feature + rt + gc + map + C(chromosome) + C(donor)"),
          ("adjusted_genes300", "events ~ feature + rt + gc + map + genes300 + C(chromosome) + C(donor)"))
CONTROLS = ["genes_1mb"]


def assign_windows(ev, w):
    """Index of the eligible window with the largest overlap with each event's [lo, hi] (-1 if none)."""
    out = np.full(len(ev), -1)
    by = {c: (g.index.values, g.start.values, g.end.values) for c, g in w.groupby("chromosome")}
    for i, (c, lo, hi) in enumerate(zip(ev.chromosome, ev.lo, ev.hi)):
        if c not in by:
            continue
        idx, S, E = by[c]
        ov = np.minimum(E, hi) - np.maximum(S, lo)
        if (ov >= 0).any():
            out[i] = idx[np.argmax(ov)]
    return out


def feature_effect(df, formula):
    """IRR, 95% CI and p of `feature` with the conservative SE (max of model-based and clustered)."""
    from scipy.stats import norm
    m, fam, alpha = fit(df, formula)
    b = m.params["feature"]
    se_cl = m.bse["feature"]
    fam_obj = sm.families.Poisson() if fam == "poisson" else sm.families.NegativeBinomial(alpha=alpha)
    se_mb = smf.glm(formula, df, offset=df.log_usable, family=fam_obj).fit().bse["feature"]
    se, src = (se_cl, "cluster") if se_cl >= se_mb else (se_mb, "model")
    return dict(family=fam, alpha=alpha, irr=np.exp(b), irr_lo=np.exp(b - 1.96 * se), irr_hi=np.exp(b + 1.96 * se),
                p=2 * norm.sf(abs(b / se)), se_source=src,
                irr_rt=np.exp(m.params["rt"]) if "rt" in m.params else np.nan,
                irr_genes300=np.exp(m.params["genes300"]) if "genes300" in m.params else np.nan)


def fit(df, formula):
    # alpha by the Cameron-Trivedi auxiliary regression on Poisson residuals:
    # ((y - mu)^2 - y) / mu = alpha * mu + e  (OLS through the origin); NB2 GLM with that alpha.
    # The discrete NB maximum likelihood overflowed with the fixed effects and offset.
    cl = dict(cov_type="cluster", cov_kwds={"groups": pd.factorize(df.chromosome)[0]})
    po = smf.glm(formula, df, offset=df.log_usable, family=sm.families.Poisson()).fit()
    mu, y = po.mu, df["events"].values
    aux = ((y - mu) ** 2 - y) / mu
    alpha = float((aux * mu).sum() / (mu ** 2).sum())
    if not np.isfinite(alpha) or alpha <= 1e-6:
        return smf.glm(formula, df, offset=df.log_usable, family=sm.families.Poisson()).fit(**cl), "poisson", alpha
    m = smf.glm(formula, df, offset=df.log_usable, family=sm.families.NegativeBinomial(alpha=alpha)).fit(**cl)
    return m, "negbin", alpha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--rt", required=True)
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    w = windows(a.events, a.rt, a.genes)
    fs = feature_sets(a.cfs, a.genes)
    for name in SETS:
        w[name] = [int(overlaps(fs[name], c, s, e)[0]) for c, s, e in zip(w.chromosome, w.start, w.end)]
    rows = []
    for b in BINS:
        ev = pd.read_csv(f"{a.events}/events_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        ev = ev[ev.type == "interstitial"]
        donors = sorted(ev.donor.unique())
        for d in ("loss", "gain", "all"):
            for s in ("main", "E3"):
                x = ev if d == "all" else ev[ev.direction == d]
                x = x[x.main_def] if s == "main" else x[x.stratum == s]
                wi = assign_windows(x, w)
                x = x.assign(widx=wi)[wi >= 0]
                cnt = x.groupby(["donor", "widx"]).size()
                df = pd.concat([w.assign(donor=dn, widx=w.index) for dn in donors], ignore_index=True)
                df["events"] = cnt.reindex(list(zip(df.donor, df.widx))).fillna(0).values
                for name in SETS:
                    df["feature"] = df[name]
                    n_in = int(df.events[df.feature == 1].sum())
                    n_out = int(df.events[df.feature == 0].sum())
                    for model, f in MODELS:
                        base = dict(bin_kb=b, direction=d, stratum=s, feature_set=name, model=model,
                                    n_events_input=len(wi), n_events_excluded=int((wi < 0).sum()),
                                    n_events=int(df.events.sum()), n_events_feature=n_in, n_windows=len(w),
                                    n_windows_feature=int(w[name].sum()))
                        if n_in == 0 or n_out == 0:
                            rows.append(dict(base, note="separation: no events " + ("in" if n_in == 0 else "outside") + " feature windows"))
                            continue
                        rows.append(dict(base, **feature_effect(df, f)))
        print(f"{b} kb done")
    t = pd.DataFrame(rows)
    t["q_family"] = np.nan
    t["q_global"] = np.nan
    for _, g in t.groupby(["feature_set", "model"]):
        g = g[g.p.notna()]
        t.loc[g.index, "q_family"] = bh(g.p)
    for _, g in t[(t.stratum == "main") & ~t.feature_set.isin(CONTROLS) & t.p.notna()].groupby("model"):
        t.loc[g.index, "q_global"] = bh(g.p)
    t.round(5).to_csv(f"{a.out}/test_b.tsv", sep="\t", index=False)
    show = t[(t.stratum == "main") & (t.direction != "gain")]
    print(show.pivot_table(index=["feature_set", "direction"], columns=["bin_kb", "model"], values="irr").round(2).to_string())
    print("families:", t.family.value_counts().to_dict(), "| SE used:", t.se_source.value_counts().to_dict(), "| separation rows:", t.p.isna().sum())
    print("events excluded (no eligible window):", t.drop_duplicates(["bin_kb", "direction", "stratum"])[["bin_kb", "direction", "stratum", "n_events_input", "n_events_excluded"]].to_string(index=False))


if __name__ == "__main__":
    main()

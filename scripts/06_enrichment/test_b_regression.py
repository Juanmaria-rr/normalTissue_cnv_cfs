#!/usr/bin/env python3
"""Test B: per-window regression of interstitial event counts on feature overlap + covariates.

Windows: usable 1 Mb bins (autosomes, applyFilters use == TRUE), outside pericentromeric
blocks +/- 1 Mb; the same windows for events of every bin size.
Rows: donor x window. Response: number of interstitial events (subset by direction and
stratum) whose midpoint falls in the window.
Covariates per window (from the 50 kb annotation of sWGS_relativeCN, averaged):
  feature    1 if the window overlaps the feature set (binary)
  rt         mean replication timing (z-scored)
  gc         mean GC % (z-scored)
  map        mean mappability (z-scored)
  genes300   fraction of the window inside genes > 300 kb (refGene, merged)
  offset     log(fraction of the window's 50 kb windows with use == TRUE)
Models (statsmodels GLM, NB2 with alpha from the Cameron-Trivedi auxiliary regression on the
Poisson fit; Poisson with HC1 robust SE when alpha <= 0, flagged in `family`):
  unadjusted  events ~ feature + donor + offset
  adjusted    events ~ feature + rt + gc + map + genes300 + chromosome + donor + offset
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
from test_a_permutation import feature_sets, merge, overlaps  # noqa: E402

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


def fit(df, formula):
    # alpha by the Cameron-Trivedi auxiliary regression on Poisson residuals:
    # ((y - mu)^2 - y) / mu = alpha * mu + e  (OLS through the origin); NB2 GLM with that alpha.
    # The discrete NB maximum likelihood overflowed with the fixed effects and offset.
    po = smf.glm(formula, df, offset=df.log_usable, family=sm.families.Poisson()).fit()
    mu, y = po.mu, df["events"].values
    aux = ((y - mu) ** 2 - y) / mu
    alpha = float((aux * mu).sum() / (mu ** 2).sum())
    if not np.isfinite(alpha) or alpha <= 1e-6:
        return smf.glm(formula, df, offset=df.log_usable, family=sm.families.Poisson()).fit(cov_type="HC1"), "poisson_hc1", alpha
    m = smf.glm(formula, df, offset=df.log_usable, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
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
                x = x.assign(wstart=((x.pos - 1) // 1_000_000) * 1_000_000 + 1)
                cnt = x.groupby(["donor", "chromosome", "wstart"]).size()
                df = pd.concat([w.assign(donor=dn) for dn in donors], ignore_index=True)
                df["events"] = cnt.reindex(list(zip(df.donor, df.chromosome, df.start))).fillna(0).values
                for name in SETS:
                    df["feature"] = df[name]
                    for model, f in (("unadjusted", "events ~ feature + C(donor)"),
                                     ("adjusted", "events ~ feature + rt + gc + map + genes300 + C(chromosome) + C(donor)")):
                        m, fam, alpha = fit(df, f)
                        ci = m.conf_int().loc["feature"]
                        rows.append(dict(bin_kb=b, direction=d, stratum=s, feature_set=name, model=model, family=fam,
                                         alpha=alpha, n_events=int(df.events.sum()), n_windows=len(w),
                                         irr=np.exp(m.params["feature"]), irr_lo=np.exp(ci[0]), irr_hi=np.exp(ci[1]),
                                         p=m.pvalues["feature"],
                                         irr_rt=np.exp(m.params["rt"]) if "rt" in m.params else np.nan,
                                         irr_genes300=np.exp(m.params["genes300"]) if "genes300" in m.params else np.nan))
        print(f"{b} kb done")
    t = pd.DataFrame(rows)
    t["q"] = np.nan
    for (name, model), g in t.groupby(["feature_set", "model"]):
        pv = g.p.values; o = np.argsort(pv); m = len(pv)
        qv = np.minimum.accumulate((pv[o] * m / np.arange(1, m + 1))[::-1])[::-1]
        qq = np.empty(m); qq[o] = np.minimum(qv, 1)
        t.loc[g.index, "q"] = qq
    t.round(5).to_csv(f"{a.out}/test_b.tsv", sep="\t", index=False)
    show = t[(t.stratum == "main") & (t.direction != "gain")]
    print(show.pivot_table(index=["feature_set", "direction"], columns=["bin_kb", "model"], values="irr").round(2).to_string())
    print("families:", t.family.value_counts().to_dict())


if __name__ == "__main__":
    main()

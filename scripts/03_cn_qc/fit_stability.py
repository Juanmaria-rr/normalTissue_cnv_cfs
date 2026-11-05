#!/usr/bin/env python3
"""Stability of the absolute fit across bins in altered cells, and its effect on breakends.

Cells: depth_ok at 500 kb and >= 1 retained non-neutral segment at 500 kb (call_breakends.R,
1.3/0.7, >= 1 Mb). For each cell and bin:
  ploidy_best      unfiltered best-clonality ploidy (cn_qc_metrics.R)
  clon_gain_vs2    clonality at ploidy 2.0 (best purity) / best clonality, from
                   <abscn>/<donor>_<bin>kb/<cell>_resfitting.txt; 1 = the fit is no better
                   than diploid, large = the chosen ploidy fits much better
  uncertain_frac   fraction of segment length whose CN = P x rel is > 0.3 from an integer
  n_rel / n_abs / n_shared  breakends by the ratio rule and by integer CN = round(P x rel)
                   (same merging and >= 1 Mb rule; shared = same chromosome within 2 bins)
Stability class per cell (bins with depth_ok only):
  stable_diploid      every ploidy within 0.15 of 2.0
  stable_nondiploid   every ploidy within 0.15 of each other, not 2
  unstable            ploidies differ by > 0.15
  single_bin          depth_ok in only one bin
Approximation: absolute CN is taken as P x rel (purity 1, rel normalised by the median
instead of the mean used by fitting_absCN), so breakends near half-integers may differ
slightly from the exact conversion of rel_vs_abs_review.py.

Usage: python3 scripts/03_cn_qc/fit_stability.py --qc results/pilot_cn/qc --abscn results/pilot_cn/absCN \
    --out results/pilot_cn/fit_stability
"""
import argparse
import os

import numpy as np
import pandas as pd

BINS = [200, 500, 1000]


def clon_gain(abscn, cell, b):
    f = f"{abscn}/{cell.split('_')[0]}_{b}kb/{cell}_resfitting.txt"
    if not os.path.exists(f):
        return np.nan
    r = pd.read_csv(f, sep="\t")
    c2 = r[np.isclose(r.ploidy, 2.0)].clonality.min()
    return c2 / max(r.clonality.min(), 1e-9)


def merge_breakends(s, state, bin_bp):
    # s: one cell's segments sorted; merge equal consecutive states per chromosome, drop < 1 Mb
    k = s.assign(state=state)
    brk = np.r_[True, (k.chromosome.values[1:] != k.chromosome.values[:-1]) | (k.state.values[1:] != k.state.values[:-1])]
    k["g"] = np.cumsum(brk)
    m = k.groupby("g").agg(chromosome=("chromosome", "first"), start=("start", "min"), end=("end", "max"),
                           state=("state", "first"))
    m = m[(m.end - m.start + 1) >= 1e6].reset_index(drop=True)
    idx = np.flatnonzero((m.chromosome.values[1:] == m.chromosome.values[:-1]) & (m.state.values[1:] != m.state.values[:-1]))
    return pd.DataFrame({"chromosome": m.chromosome.values[idx], "pos": m.start.values[idx + 1]})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--abscn", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    q = {b: pd.read_csv(f"{a.qc}/cn_qc_{b}kb.tsv", sep="\t").set_index("cell_id") for b in BINS}
    seg = {b: pd.read_csv(f"{a.qc}/breakends_{b}kb_segments.tsv", sep="\t", dtype={"chromosome": str}) for b in BINS}
    c5 = pd.read_csv(f"{a.qc}/breakends_500kb_cells.tsv", sep="\t").set_index("cell_id")
    alt = [c for c in q[500].index[q[500].depth_ok == 1] if (c5.loc[c, "n_gain_seg"] + c5.loc[c, "n_loss_seg"]) > 0]
    print(len(alt), "altered cells at 500 kb")

    rows = []
    for b in BINS:
        sb = seg[b][seg[b].cell_id.isin(alt)]
        for c, s in sb.groupby("cell_id", sort=False):
            if c not in q[b].index:
                continue
            P = q[b].loc[c, "ploidy_best"]
            cn = P * s.rel.values
            ln = (s.end - s.start + 1).values
            unc = ln[np.abs(cn - np.round(cn)) > 0.3].sum() / ln.sum()
            rstate = np.where(s.rel > 1.3, 1, np.where(s.rel < 0.7, -1, 0))
            rb = merge_breakends(s, rstate, b * 1000)
            ab = merge_breakends(s, np.round(cn).astype(int), b * 1000)
            tol = 2 * b * 1000
            sh = sum(((ab.chromosome == r.chromosome) & (np.abs(ab.pos - r.pos) <= tol)).any() for r in rb.itertuples())
            sh_a = sum(((rb.chromosome == r.chromosome) & (np.abs(rb.pos - r.pos) <= tol)).any() for r in ab.itertuples())
            rows.append(dict(cell_id=c, donor=c.split("_")[0], bin_kb=b, depth_ok=int(q[b].loc[c, "depth_ok"]),
                             ploidy_best=P, clon_gain_vs2=clon_gain(a.abscn, c, b), uncertain_frac=unc,
                             n_rel=len(rb), n_abs=len(ab), rel_shared=sh, abs_shared=sh_a))
    t = pd.DataFrame(rows)
    t.round(4).to_csv(f"{a.out}/fit_stability_per_bin.tsv", sep="\t", index=False)

    cls = []
    for c, g in t[t.depth_ok == 1].groupby("cell_id"):
        p = g.ploidy_best.values
        if len(p) < 2:
            k = "single_bin"
        elif np.ptp(p) <= 0.15:
            k = "stable_diploid" if np.all(np.abs(p - 2) <= 0.15) else "stable_nondiploid"
        else:
            k = "unstable"
        cls.append(dict(cell_id=c, donor=g.donor.iloc[0], n_bins=len(p), ploidies="/".join(f"{x:.1f}" for x in g.sort_values("bin_kb").ploidy_best),
                        stability=k))
    cl = pd.DataFrame(cls)
    cl.to_csv(f"{a.out}/fit_stability_per_cell.tsv", sep="\t", index=False)
    print(pd.crosstab(cl.donor, cl.stability, margins=True))

    tt = t[t.depth_ok == 1].merge(cl[["cell_id", "stability"]], on="cell_id", how="left")
    summ = tt.groupby(["bin_kb", "stability"]).agg(cells=("cell_id", "size"), rel=("n_rel", "sum"), abs=("n_abs", "sum"),
                                                     rel_shared=("rel_shared", "sum"), abs_shared=("abs_shared", "sum"),
                                                     median_uncertain=("uncertain_frac", "median"),
                                                     median_clon_gain=("clon_gain_vs2", "median")).reset_index()
    summ["abs_only"] = summ["abs"] - summ.abs_shared
    summ["rel_only"] = summ.rel - summ.rel_shared
    summ.round(3).to_csv(f"{a.out}/fit_stability_summary.tsv", sep="\t", index=False)
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

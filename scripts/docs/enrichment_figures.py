#!/usr/bin/env python3
"""Figures of the CFS enrichment analysis (pilot), written to --out as PNG.

  fig_a_forest      test A fold (95% bootstrap CI) per feature set, interstitial events,
                    main stratum (E2-E4), losses and gains, three bins; grey bar = 2.5-97.5%
                    range of the fold over RT-matched shuffles of that set (losses / all)
  fig_a_strata      test A fold (95% bootstrap CI) per deviation stratum (E1-E4), losses and all
  fig_a_sensitivity test A fold (95% bootstrap CI) for losses (main): main set, no cuts, same
                    cells, and leaving out each donor
  fig_b_irr         test B IRR (95% CI), unadjusted vs adjusted, losses (main)
  fig_c_donors      test C ratio (95% cell bootstrap CI) per donor, losses (main)
  fig_regions       loss breakend intervals (cells of the main set, 500 kb) around FHIT,
                    LRP1B and IMMP2L, with MiDAS sites, HumCFS bands, gene spans and RT

Usage: python3 scripts/docs/enrichment_figures.py --res results/enrichment/pilot --out docs/enrichment
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SETS = ["macheret2020:all", "macheret2020:U2OS", "macheret2020:HeLa", "macheret2020:HS68", "ji2020:all",
        "ji2020:U2OS", "ji2020:HeLa", "ji2020:HBE", "humcfs:all", "letallec2013:all", "letallec2013:epithelial",
        "letallec2013:MCF10A", "brison2019:all", "sarni2020:all", "wilson2015:all", "simpson2021:all",
        "consensusA:core", "consensusA:reproducible", "consensusA:all", "breast_predicted:all", "A_plus_C:all", "union_all:all",
        "genes_300kb", "genes_1mb"]
BCOL = {200: "#8172b2", 500: "#c44e52", 1000: "#4c72b0"}
DCOL = {"P18": "#4c72b0", "P20": "#55a868", "P36": "#8172b2", "P07": "#c44e52", "P08": "#dd8452"}
plt.rcParams.update({"figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False, "font.size": 8,
                     "axes.titlesize": 9, "axes.titlelocation": "left"})


def fig_a_forest(res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    sh = pd.read_csv(f"{res}/main/test_a_shuffle_controls.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11, 11), sharey=True)
    for ax, d in zip(axes, ("loss", "gain")):
        for j, name in enumerate(SETS):
            y0 = len(SETS) - j
            for k, b in enumerate((200, 500, 1000)):
                y = y0 + (k - 1) * 0.25
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == d) & (t.stratum == "main")]
                if r.empty or not np.isfinite(r.fold.iloc[0]):
                    continue
                r = r.iloc[0]
                s = sh[(sh.feature_set == name) & (sh.bin_kb == b) & (sh.direction == ("loss" if d == "loss" else "all"))].fold
                if len(s):
                    ax.plot([s.quantile(.025), s.quantile(.975)], [y, y], color="#d0d0d0", lw=5, solid_capstyle="butt", zorder=0)
                ax.plot([r.fold_lo, r.fold_hi], [y, y], color=BCOL[b], lw=1)
                ax.scatter(r.fold, y, s=14, color=BCOL[b], zorder=3, label=f"{b} kb" if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":")
        ax.set_xlim(0, 6); ax.set_title(f"Interstitial {'losses' if d == 'loss' else 'gains'} (E2–E4)")
        ax.set_xlabel("fold (observed / expected)")
    axes[0].set_yticks(range(len(SETS), 0, -1), SETS)
    axes[0].legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    return fig


def fig_a_strata(res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    sel = ["macheret2020:all", "macheret2020:HS68", "humcfs:all", "genes_1mb"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 5), sharey=True)
    for i, d in enumerate(("loss", "all")):
        for k, b in enumerate((200, 500, 1000)):
            ax = axes[i, k]
            for k, name in enumerate(sel):
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == d) & t.stratum.isin(["E1", "E2", "E3", "E4"])]
                r = r.set_index("stratum").reindex(["E1", "E2", "E3", "E4"])
                x = np.arange(4) + (k - 1.5) * 0.12  # dodge the sets so the CIs do not overlap
                ax.errorbar(x, r.fold, yerr=[r.fold - r.fold_lo, r.fold_hi - r.fold], fmt="o-", ms=3, lw=.8,
                            capsize=2, elinewidth=.8, label=name)
            ax.axhline(1, color="#666666", lw=.7, ls=":")
            ax.set_xticks(range(4), ["E1\n0.2–0.3", "E2\n0.3–0.4", "E3\n0.4–0.6", "E4\n≥0.6"])
            ax.set_title(f"{b} kb · {'losses' if d == 'loss' else 'all directions'}")
            ax.set_ylim(0, 8)
    axes[0, 0].set_ylabel("fold (95% CI)"); axes[1, 0].set_ylabel("fold (95% CI)")
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout()
    return fig


def fig_a_sensitivity(res):
    runs = [("main", f"{res}/main/test_a.tsv"), ("no cuts", f"{res}/nocuts/test_a.tsv"),
            ("same cells", f"{res}/samecells/test_a.tsv")] + \
           [(f"without {d}", f"{res}/main/test_a_without_{d}.tsv") for d in ("P07", "P08", "P18", "P20", "P36")]
    sel = ["macheret2020:all", "macheret2020:HS68", "ji2020:all", "humcfs:all", "genes_1mb"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), sharey=True)
    for ax, b in zip(axes, (200, 500, 1000)):
        for j, name in enumerate(sel):
            for i, (lab, f) in enumerate(runs):
                t = pd.read_csv(f, sep="\t")
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == "loss") & (t.stratum == "main")]
                if r.empty:
                    continue
                y = len(sel) - j + (i - 3.5) * 0.1
                r = r.iloc[0]
                ax.plot([r.fold_lo, r.fold_hi], [y, y], color=plt.cm.tab10(i), lw=.8)
                ax.scatter(r.fold, y, s=10, color=plt.cm.tab10(i), label=lab if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":")
        ax.set_title(f"{b} kb · losses (E2–E4)"); ax.set_xlabel("fold (95% CI)"); ax.set_xlim(0, 8)
    axes[0].set_yticks(range(len(sel), 0, -1), sel)
    axes[2].legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    return fig


def fig_b_irr(res):
    t = pd.read_csv(f"{res}/main/test_b.tsv", sep="\t")
    sets = t.feature_set.unique().tolist()
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), sharey=True)
    for ax, b in zip(axes, (200, 500, 1000)):
        for j, name in enumerate(sets):
            for k, m in enumerate(("unadjusted", "adjusted")):
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == "loss") & (t.stratum == "main") & (t.model == m)]
                if r.empty:
                    continue
                r = r.iloc[0]; y = len(sets) - j + (k - .5) * .3
                col = "#999999" if m == "unadjusted" else "#c44e52"
                ax.plot([r.irr_lo, r.irr_hi], [y, y], color=col, lw=1)
                ax.scatter(r.irr, y, s=12, color=col, label=("unadjusted" if m == "unadjusted" else "adjusted") if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":"); ax.set_xscale("log"); ax.set_xlim(.1, 20)
        ax.set_title(f"{b} kb · losses (E2–E4)"); ax.set_xlabel("IRR (log scale)")
    axes[0].set_yticks(range(len(sets), 0, -1), sets)
    axes[2].legend(fontsize=7)
    fig.tight_layout()
    return fig


def fig_c_donors(res):
    t = pd.read_csv(f"{res}/main/test_c.tsv", sep="\t")
    sel = ["macheret2020:all", "macheret2020:HS68", "ji2020:all", "humcfs:all"]
    donors = ["pooled", "P07", "P08", "P18", "P20", "P36"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, b in zip(axes, (200, 500, 1000)):
        for j, name in enumerate(sel):
            for i, dn in enumerate(donors):
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == "loss") & (t.donor == dn)]
                if r.empty or not np.isfinite(r.ratio.iloc[0]):
                    continue
                r = r.iloc[0]; y = len(sel) - j + (i - 2.5) * 0.12
                col = "#222222" if dn == "pooled" else DCOL[dn]
                ax.plot([r.ratio_lo, r.ratio_hi], [y, y], color=col, lw=.8)
                ax.scatter(r.ratio, y, s=10, color=col, label=dn if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":"); ax.set_xlim(0, 8)
        ax.set_title(f"{b} kb · losses (E2–E4), per breakend"); ax.set_xlabel("observed / expected")
    axes[0].set_yticks(range(len(sel), 0, -1), sel)
    axes[2].legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    return fig


def fig_regions(res):
    be = pd.read_csv(f"{res}/main/breakends_500kb.tsv", sep="\t", dtype={"chromosome": str})
    be = be[(be.type == "interstitial") & (be.direction == "loss") & be.stratum.isin(["E2", "E3", "E4"])]
    cfs = pd.read_csv(f"{ROOT}/data/features/cfs/cfs_catalogs_hg19.tsv", sep="\t")
    cfs["chromosome"] = cfs.chrom.str.replace("chr", "")
    genes = pd.read_csv(f"{ROOT}/data/features/genes/large_genes_300kb.tsv", sep="\t", dtype={"chromosome": str})
    rt = pd.read_csv(f"{res}/main/usable_bins_rt_500kb.tsv", sep="\t", dtype={"chromosome": str})
    regions = [("FHIT", "3", 56e6, 66e6), ("LRP1B", "2", 138e6, 148e6), ("IMMP2L", "7", 106e6, 116e6)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6))
    for ax, (g, c, s, e) in zip(axes, regions):
        y = 0
        rr = rt[(rt.chromosome == c) & (rt.start >= s) & (rt.end <= e)]
        ax.plot((rr.start + rr.end) / 2e6, 7.4 + rr.rt / 40, color="#55a868", lw=1)
        ax.text(s / 1e6, 7.4 + 90 / 40, "replication timing (up = early)", fontsize=6, color="#55a868", va="top")
        for k, (lab, sub) in enumerate((("MiDAS HS68", cfs[(cfs.catalog == "macheret2020") & (cfs.cell_type == "HS68")]),
                                        ("MiDAS U2OS+HeLa", cfs[(cfs.catalog == "macheret2020") & (cfs.cell_type != "HS68")]),
                                        ("HumCFS", cfs[cfs.catalog == "humcfs"]),
                                        ("genes > 300 kb", genes))):
            yy = 6.6 - k * .5
            x = sub[(sub.chromosome == c) & (sub.end >= s) & (sub.start <= e)]
            for r in x.itertuples():
                ax.plot([max(r.start, s) / 1e6, min(r.end, e) / 1e6], [yy, yy], lw=4, color="#444444" if k < 3 else "#b8962e",
                        solid_capstyle="butt")
            ax.text(s / 1e6, yy + .15, lab, fontsize=6)
        x = be[(be.chromosome == c) & (be.hi >= s) & (be.lo <= e)].sort_values("pos")
        for i, r in enumerate(x.itertuples()):
            yy = 4.6 - i * 0.12
            ax.plot([r.lo / 1e6, r.hi / 1e6], [yy, yy], color=DCOL.get(r.donor, "k"), lw=2)
        ax.set_xlim(s / 1e6, e / 1e6); ax.set_yticks([]); ax.set_ylim(min(3.5, 4.4 - len(x) * .12), 9.8)
        ax.set_title(f"{g} · chr{c}:{s / 1e6:.0f}–{e / 1e6:.0f} Mb · {len(x)} loss breakends")
        ax.set_xlabel(f"chr{c} (Mb, hg19)")
    for dn, col in DCOL.items():
        axes[2].plot([], [], color=col, lw=2, label=dn)
    axes[2].legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    return fig


FIGS = {"fig_a_forest": fig_a_forest, "fig_a_strata": fig_a_strata, "fig_a_sensitivity": fig_a_sensitivity,
        "fig_b_irr": fig_b_irr, "fig_c_donors": fig_c_donors, "fig_regions": fig_regions}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for name, f in FIGS.items():
        fig = f(a.res)
        fig.savefig(f"{a.out}/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
        print("wrote", name)


if __name__ == "__main__":
    main()

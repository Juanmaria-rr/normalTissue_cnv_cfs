#!/usr/bin/env python3
"""Figures of the CFS enrichment analysis (pilot), written to --out as PNG.

  fig_a_forest      test A fold (95% bootstrap CI) per feature set, interstitial events,
                    main stratum (E2-E4), losses, gains and all directions, three bins; grey bar =
                    2.5-97.5% range of the fold over RT-matched shuffles of that set and direction;
                    hollow marker = expected < 2 events (fold dominated by counting noise)
  fig_a_strata      test A fold (95% bootstrap CI) per deviation stratum (E1-E4), losses, gains and all directions
  fig_a_sensitivity test A fold (95% bootstrap CI) in every robustness variant (robustness.py):
                    main, no cuts, same cells, cells500, refined, linkage 1 / 5 Mb, collapse on
                    E2-E4, leave-one-donor-out
  fig_top_donors    per-donor observed / expected (exact Poisson 95% CI) of the strongest test A
                    row (top_row(): smallest permutation p among main-stratum CFS rows)
  fig_b_irr         test B IRR (95% CI, conservative SE), crude / adjusted /
                    adjusted + genes > 300 kb
  fig_c_donors      test C ratio (95% cell bootstrap CI; one-sided when O = 0) per donor
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

DLAB = {"loss": "losses", "gain": "gains", "all": "all directions"}
DIRS = ("loss", "gain", "all")
CONTROLS = ["genes_300kb", "genes_1mb"]
MCOL = {"crude": "#999999", "adjusted": "#c44e52", "adjusted_genes300": "#4c72b0"}


def top_row(t):
    """Strongest test A row: main stratum, CFS sets (controls excluded), smallest p_perm then p_enrich."""
    x = t[(t.stratum == "main") & ~t.feature_set.isin(CONTROLS) & t.p_perm.notna()]
    return x.sort_values(["p_perm", "p_enrich"]).iloc[0]


def poisson_ci(o, alpha=.05):
    from scipy.stats import chi2
    lo = chi2.ppf(alpha / 2, 2 * o) / 2 if o > 0 else 0.0
    hi = chi2.ppf(1 - alpha / 2, 2 * (o + 1)) / 2
    return lo, hi


def fig_a_forest(res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    sh = pd.read_csv(f"{res}/main/test_a_shuffle_controls.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(15, 11), sharey=True)
    for ax, d in zip(axes, DIRS):
        for j, name in enumerate(SETS):
            y0 = len(SETS) - j
            for k, b in enumerate((200, 500, 1000)):
                y = y0 + (k - 1) * 0.25
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == d) & (t.stratum == "main")]
                if r.empty or not np.isfinite(r.fold.iloc[0]):
                    continue
                r = r.iloc[0]
                s = sh[(sh.feature_set == name) & (sh.bin_kb == b) & (sh.direction == d)].fold
                if len(s):
                    ax.plot([s.quantile(.025), s.quantile(.975)], [y, y], color="#d0d0d0", lw=5, solid_capstyle="butt", zorder=0)
                ax.plot([r.fold_lo, r.fold_hi], [y, y], color=BCOL[b], lw=1)
                low = r.expected < 2
                ax.scatter(r.fold, y, s=14, zorder=3, label=f"{b} kb" if j == 0 else None,
                           facecolor="white" if low else BCOL[b], edgecolor=BCOL[b], linewidth=1)
        ax.axvline(1, color="#666666", lw=.7, ls=":")
        ax.set_xlim(0, 6); ax.set_title(f"Interstitial {DLAB[d]} (E2–E4)")
        ax.set_xlabel("fold (observed / expected)")
    axes[0].set_yticks(range(len(SETS), 0, -1), SETS)
    axes[0].legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    return fig


def fig_a_strata(res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    sel = ["macheret2020:all", "macheret2020:HS68", "humcfs:all", "genes_1mb"]
    fig, axes = plt.subplots(3, 3, figsize=(11, 7.5), sharey=True)
    for i, d in enumerate(DIRS):
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
            ax.set_title(f"{b} kb · {DLAB[d]}")
            ax.set_ylim(0, 8)
    for i in range(3):
        axes[i, 0].set_ylabel("fold (95% CI)")
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout()
    return fig


def fig_a_sensitivity(res):
    v = pd.read_csv(f"{res}/robustness_variants.tsv", sep="\t")
    runs = list(dict.fromkeys(v.variant))
    sel = ["macheret2020:all", "macheret2020:HS68", "ji2020:all", "humcfs:all", "consensusA:core", "genes_1mb"]
    fig, axs = plt.subplots(3, 3, figsize=(12, 15), sharey=True)
    for (row, d), (col, b) in [((i, d), (k, b)) for i, d in enumerate(DIRS) for k, b in enumerate((200, 500, 1000))]:
        ax = axs[row, col]
        for j, name in enumerate(sel):
            for i, lab in enumerate(runs):
                r = v[(v.variant == lab) & (v.feature_set == name) & (v.bin_kb == b) & (v.direction == d)]
                if r.empty or not np.isfinite(r.fold.iloc[0]):
                    continue
                y = len(sel) - j + (i - (len(runs) - 1) / 2) * 0.06
                r = r.iloc[0]
                ax.plot([r.fold_lo, r.fold_hi], [y, y], color=plt.cm.tab20(i), lw=.8)
                ax.scatter(r.fold, y, s=8, color=plt.cm.tab20(i), label=lab if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":")
        ax.set_title(f"{b} kb · {DLAB[d]} (E2–E4)"); ax.set_xlabel("fold (95% CI)"); ax.set_xlim(0, 8)
    for i in range(3):
        axs[i, 0].set_yticks(range(len(sel), 0, -1), sel)
    handles, labels = [], []
    for ax in axs.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h); labels.append(l)
    fig.legend(handles, labels, fontsize=7, loc="lower center", ncol=7, bbox_to_anchor=(.5, -.02))
    fig.tight_layout(rect=(0, .03, 1, 1))
    return fig


def fig_top_donors(res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    top = top_row(t)
    pr = pd.read_csv(f"{res}/main/test_a_event_probs.tsv", sep="\t", dtype={"chromosome": str})
    x = pr[(pr.bin_kb == top.bin_kb) & (pr.feature_set == top.feature_set) & (pr.type == "interstitial") & pr.main_def & pr.p.notna()]
    if top.direction != "all":
        x = x[x.direction == top.direction]
    g = x.groupby("donor").agg(O=("o", "sum"), E=("p", "sum"), n=("o", "size"))
    g.loc["pooled"] = [g.O.sum(), g.E.sum(), g.n.sum()]
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    for i, (dn, r) in enumerate(g.iloc[::-1].iterrows()):
        lo, hi = poisson_ci(int(r.O))
        col = "#222222" if dn == "pooled" else DCOL.get(dn, "k")
        ax.plot([lo / r.E, hi / r.E], [i, i], color=col, lw=1.2)
        ax.scatter(r.O / r.E, i, s=22, color=col, zorder=3)
        ax.text(10.3, i, f"{int(r.O)} / {r.E:.2f}  (n = {int(r.n)})", va="center", fontsize=7)
    ax.set_yticks(range(len(g)), list(g.index[::-1]))
    ax.axvline(1, color="#666666", lw=.7, ls=":"); ax.set_xlim(0, 10)
    ax.set_xlabel("observed / expected (exact Poisson 95% CI)")
    ax.set_title(f"{top.feature_set} · {int(top.bin_kb)} kb · interstitial {DLAB[top.direction]} (E2–E4) · per donor   [O / E (events)]")
    fig.tight_layout()
    return fig


def fig_b_irr(res):
    t = pd.read_csv(f"{res}/main/test_b.tsv", sep="\t")
    sets = t.feature_set.unique().tolist()
    fig, axs = plt.subplots(3, 3, figsize=(12, 15), sharey=True)
    for (row, d), (col, b) in [((i, d), (k, b)) for i, d in enumerate(DIRS) for k, b in enumerate((200, 500, 1000))]:
        ax = axs[row, col]
        for j, name in enumerate(sets):
            for k, m in enumerate(MCOL):
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == d) & (t.stratum == "main") & (t.model == m)]
                if r.empty or not np.isfinite(r.irr.iloc[0]):
                    continue
                r = r.iloc[0]; y = len(sets) - j + (k - 1) * .25
                ax.plot([r.irr_lo, r.irr_hi], [y, y], color=MCOL[m], lw=1)
                ax.scatter(r.irr, y, s=10, color=MCOL[m], label=m if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":"); ax.set_xscale("log"); ax.set_xlim(.1, 20)
        ax.set_title(f"{b} kb · {DLAB[d]} (E2–E4)"); ax.set_xlabel("IRR (log scale)")
    for i in range(3):
        axs[i, 0].set_yticks(range(len(sets), 0, -1), sets, fontsize=7)
    axs[0, 2].legend(fontsize=7)
    fig.tight_layout()
    return fig


def fig_c_donors(res):
    t = pd.read_csv(f"{res}/main/test_c.tsv", sep="\t")
    sel = ["macheret2020:all", "macheret2020:HS68", "ji2020:all", "humcfs:all"]
    donors = ["pooled", "P07", "P08", "P18", "P20", "P36"]
    fig, axs = plt.subplots(3, 3, figsize=(12, 10), sharey=True)
    for (row, d), (col, b) in [((i, d), (k, b)) for i, d in enumerate(DIRS) for k, b in enumerate((200, 500, 1000))]:
        ax = axs[row, col]
        for j, name in enumerate(sel):
            for i, dn in enumerate(donors):
                r = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == d) & (t.donor == dn)]
                if r.empty or not np.isfinite(r.ratio.iloc[0]):
                    continue
                r = r.iloc[0]; y = len(sel) - j + (i - 2.5) * 0.12
                col = "#222222" if dn == "pooled" else DCOL[dn]
                ax.plot([r.ratio_lo, r.ratio_hi], [y, y], color=col, lw=.8)
                ax.scatter(r.ratio, y, s=10, color=col, label=dn if j == 0 else None)
        ax.axvline(1, color="#666666", lw=.7, ls=":"); ax.set_xlim(0, 8)
        ax.set_title(f"{b} kb · {DLAB[d]} (E2–E4), per breakend"); ax.set_xlabel("observed / expected")
    for i in range(3):
        axs[i, 0].set_yticks(range(len(sel), 0, -1), sel)
    axs[0, 2].legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    return fig


def fig_regions(res):
    be = pd.read_csv(f"{res}/main/breakends_500kb.tsv", sep="\t", dtype={"chromosome": str})
    be = be[(be.type == "interstitial") & be.direction.isin(["loss", "gain"]) & be.stratum.isin(["E2", "E3", "E4"])]
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
            ax.plot([r.lo / 1e6, r.hi / 1e6], [yy, yy], color=DCOL.get(r.donor, "k"), lw=2, ls="-" if r.direction == "loss" else ":")
        ax.set_xlim(s / 1e6, e / 1e6); ax.set_yticks([]); ax.set_ylim(min(3.5, 4.4 - len(x) * .12), 9.8)
        nl, ng = int((x.direction == "loss").sum()), int((x.direction == "gain").sum())
        ax.set_title(f"{g} · chr{c}:{s / 1e6:.0f}–{e / 1e6:.0f} Mb · {nl} loss, {ng} gain breakends")
        ax.set_xlabel(f"chr{c} (Mb, hg19)")
    for dn, col in DCOL.items():
        axes[2].plot([], [], color=col, lw=2, label=dn)
    axes[2].plot([], [], color="#666666", lw=2, ls="-", label="loss")
    axes[2].plot([], [], color="#666666", lw=2, ls=":", label="gain")
    axes[2].legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    return fig


FIGS = {"fig_a_forest": fig_a_forest, "fig_a_strata": fig_a_strata, "fig_a_sensitivity": fig_a_sensitivity,
        "fig_top_donors": fig_top_donors, "fig_b_irr": fig_b_irr, "fig_c_donors": fig_c_donors, "fig_regions": fig_regions}


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

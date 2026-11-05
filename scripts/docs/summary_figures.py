#!/usr/bin/env python3
"""Summary figures of the CFS enrichment pilot for the results deck (16:9 PNG, 2400 x 1250 px).

  sum_test_a      test A at 500 kb, main stratum, losses and gains: fold with event-bootstrap
                  95% CI per set; grey band = 2.5-97.5% of the fold over 5,000 shuffled sets;
                  filled marker = p perm. < 0.05; hollow = expected < 2 events
  sum_top         strongest test A row (top_row()): per-donor observed / expected (exact Poisson
                  95% CI) and its fold in every robustness variant at the three bins
  sum_concordance 500 kb losses, every CFS set: test A fold vs test B adjusted IRR and vs test C
                  pooled O/E (log axes)
  sum_power       500 kb losses: projected smallest fold detectable with 80% power at 49 donors
                  (alpha 0.05 and global alpha) next to the pilot fold, sets ordered by expected events

Inputs: results/enrichment/pilot/main/test_a.tsv, test_a_shuffle_controls.tsv,
test_a_event_probs.tsv, test_b.tsv, test_c.tsv; robustness_variants.tsv; power_scaleup.tsv.

Usage: python3 scripts/docs/summary_figures.py --res results/enrichment/pilot --out docs/slides/figs
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from enrichment_figures import poisson_ci, top_row  # noqa: E402

INK, MUTED, RULE, BG = "#1B2A33", "#5B6B75", "#D5DCDF", "#F6F7F4"
LOSS, GAIN, ACCENT = "#C2571A", "#2F6DB5", "#0F6E73"
DONOR = {"P07": "#C2571A", "P08": "#2F6DB5", "P18": "#0F6E73", "P20": "#8A5A9E", "P36": "#B08A1E"}
CONTROLS = ["genes_300kb", "genes_1mb"]
SETS = ["macheret2020:all", "macheret2020:HS68", "macheret2020:HeLa", "macheret2020:U2OS", "ji2020:all", "ji2020:HBE",
        "ji2020:HeLa", "ji2020:U2OS", "humcfs:all", "letallec2013:all", "letallec2013:epithelial", "letallec2013:MCF10A",
        "brison2019:all", "sarni2020:all", "wilson2015:all", "simpson2021:all", "consensusA:core", "consensusA:reproducible",
        "consensusA:all", "breast_predicted:all", "A_plus_C:all", "union_all:all", "genes_300kb", "genes_1mb"]


def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 200, "font.size": 11, "font.family": "sans-serif",
                         "font.sans-serif": ["IBM Plex Sans", "Helvetica Neue", "Arial", "DejaVu Sans"],
                         "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK, "text.color": INK,
                         "axes.titlelocation": "left", "axes.titlesize": 13, "axes.titleweight": "bold",
                         "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                         "grid.color": RULE, "grid.linewidth": .6})
    return plt


def sum_test_a(plt, res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    sh = pd.read_csv(f"{res}/main/test_a_shuffle_controls.tsv", sep="\t")
    t = t[(t.stratum == "main") & (t.bin_kb == 500)]
    sh = sh[sh.bin_kb == 500]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.25), sharey=True)
    for ax, (d, col, lab) in zip(axes, [("loss", LOSS, "Losses"), ("gain", GAIN, "Gains")]):
        for j, s in enumerate(SETS):
            y = len(SETS) - j
            r = t[(t.feature_set == s) & (t.direction == d)]
            if r.empty or not np.isfinite(r.fold.iloc[0]):
                continue
            r = r.iloc[0]
            f = sh[(sh.feature_set == s) & (sh.direction == d)].fold.dropna()
            if len(f):
                ax.plot([f.quantile(.025), f.quantile(.975)], [y, y], color=RULE, lw=7, solid_capstyle="butt", zorder=0)
            ax.plot([r.fold_lo, min(r.fold_hi, 7.9)], [y, y], color=col, lw=1.4)
            sig = np.isfinite(r.p_perm) and r.p_perm < .05
            ax.scatter(r.fold, y, s=42, zorder=3, facecolor=col if sig else BG, edgecolor=col, linewidth=1.4,
                       alpha=.45 if r.expected < 2 else 1)
            ax.text(8.05, y, f"{int(r.observed)}/{r.expected:.1f}", va="center", fontsize=9, color=MUTED)
        ax.axvline(1, color=MUTED, lw=.8, ls=(0, (3, 3)))
        ax.set_xlim(0, 8); ax.set_xlabel("fold = observed / expected events")
        ax.set_title(f"{lab} · interstitial E2–E4 · 500 kb", color=col)
        ax.grid(axis="x")
    axes[0].set_yticks(range(len(SETS), 0, -1), SETS, fontsize=10)
    axes[1].text(8.05, len(SETS) + .9, "O/E", fontsize=9, color=MUTED)
    fig.tight_layout()
    return fig


def sum_top(plt, res):
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    top = top_row(t)
    pr = pd.read_csv(f"{res}/main/test_a_event_probs.tsv", sep="\t", dtype={"chromosome": str})
    x = pr[(pr.bin_kb == top.bin_kb) & (pr.feature_set == top.feature_set) & (pr.type == "interstitial") & pr.main_def & pr.p.notna()]
    if top.direction != "all":
        x = x[x.direction == top.direction]
    g = x.groupby("donor").agg(O=("o", "sum"), E=("p", "sum"), n=("o", "size"))
    g.loc["pooled"] = [g.O.sum(), g.E.sum(), g.n.sum()]
    v = pd.read_csv(f"{res}/robustness_variants.tsv", sep="\t")
    v = v[(v.feature_set == top.feature_set) & (v.direction == top.direction)]
    runs = list(dict.fromkeys(v.variant))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6.25), gridspec_kw={"width_ratios": [1, 1.25]})
    for i, (dn, r) in enumerate(g.iloc[::-1].iterrows()):
        lo, hi = poisson_ci(int(r.O))
        col = INK if dn == "pooled" else DONOR.get(dn, INK)
        a1.plot([lo / r.E, min(hi / r.E, 11.9)], [i, i], color=col, lw=1.6)
        a1.scatter(r.O / r.E, i, s=60, color=col, zorder=3)
        a1.text(12.2, i, f"{int(r.O)} / {r.E:.2f}", va="center", fontsize=10, color=MUTED)
    a1.set_yticks(range(len(g)), list(g.index[::-1]))
    a1.axvline(1, color=MUTED, lw=.8, ls=(0, (3, 3))); a1.set_xlim(0, 12); a1.grid(axis="x")
    a1.set_xlabel("observed / expected (exact Poisson 95% CI)")
    a1.set_title("Per donor (500 kb)")
    bcol = {200: "#8C9AA3", 500: LOSS, 1000: INK}
    for i, lab in enumerate(runs[::-1]):
        for k, b in enumerate((200, 500, 1000)):
            r = v[(v.variant == lab) & (v.bin_kb == b)]
            if r.empty or not np.isfinite(r.fold.iloc[0]):
                continue
            r = r.iloc[0]; y = i + (k - 1) * .22
            a2.plot([r.fold_lo, min(r.fold_hi, 9.9)], [y, y], color=bcol[b], lw=1.2)
            a2.scatter(r.fold, y, s=26, color=bcol[b], zorder=3, label=f"{'1 Mb' if b == 1000 else f'{b} kb'}" if i == 0 else None)
    a2.set_yticks(range(len(runs)), runs[::-1])
    a2.axvline(1, color=MUTED, lw=.8, ls=(0, (3, 3))); a2.set_xlim(0, 10); a2.grid(axis="x")
    a2.set_xlabel("fold (event-bootstrap 95% CI)")
    a2.set_title("Robustness variants, three bins")
    a2.legend(frameon=False, fontsize=10, loc="lower right")
    fig.tight_layout()
    return fig


def sum_concordance(plt, res):
    ta = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    tb = pd.read_csv(f"{res}/main/test_b.tsv", sep="\t")
    tc = pd.read_csv(f"{res}/main/test_c.tsv", sep="\t")
    a = ta[(ta.stratum == "main") & (ta.bin_kb == 500) & (ta.direction == "loss")].set_index("feature_set")
    b = tb[(tb.stratum == "main") & (tb.bin_kb == 500) & (tb.direction == "loss") & (tb.model == "adjusted")].set_index("feature_set")
    c = tc[(tc.donor == "pooled") & (tc.bin_kb == 500) & (tc.direction == "loss")].set_index("feature_set")
    label = {"macheret2020:HS68": (-150, 14), "macheret2020:all": (30, -40), "consensusA:core": (40, -70), "union_all:all": (-120, 60)}
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.25))
    for ax, other, col, lab in ((axes[0], b, "irr", "Test B: adjusted IRR"), (axes[1], c, "ratio", "Test C: pooled O/E (per breakend)")):
        common = [s for s in a.index if s in other.index and np.isfinite(a.loc[s, "fold"]) and np.isfinite(other.loc[s, col])
                  and a.loc[s, "fold"] > 0 and other.loc[s, col] > 0]
        xs, ys = a.loc[common, "fold"].values, other.loc[common, col].values
        ctrl = np.array([s in CONTROLS for s in common])
        ax.scatter(xs[~ctrl], ys[~ctrl], s=46, color=ACCENT, edgecolor=BG, linewidth=1, zorder=3, label="CFS set")
        ax.scatter(xs[ctrl], ys[ctrl], s=46, facecolor=BG, edgecolor=MUTED, linewidth=1.2, zorder=3, label="gene-size control")
        for s, xv, yv in zip(common, xs, ys):
            if s in label:
                ax.annotate(s, (xv, yv), xytext=label[s], textcoords="offset points", fontsize=9, color=INK,
                            arrowprops=dict(arrowstyle="-", color=MUTED, lw=.7, shrinkA=0, shrinkB=4))
        ax.set_xscale("log"); ax.set_yscale("log")
        lim = (.3, 6)
        ax.plot(lim, lim, color=MUTED, lw=.8, ls=(0, (3, 3)))
        ax.axhline(1, color=RULE, lw=.8); ax.axvline(1, color=RULE, lw=.8)
        ax.set_xlim(*lim); ax.set_ylim(*lim)
        ax.set_xticks([.5, 1, 2, 4], ["0.5", "1", "2", "4"]); ax.set_yticks([.5, 1, 2, 4], ["0.5", "1", "2", "4"])
        ax.set_xlabel("Test A: fold (events)"); ax.set_ylabel(lab)
        ax.set_title(f"A vs {lab.split(':')[0][-1]} · 500 kb losses")
    axes[0].legend(frameon=False, fontsize=10, loc="upper left")
    fig.tight_layout()
    return fig


def sum_power(plt, res):
    p = pd.read_csv(f"{res}/power_scaleup.tsv", sep="\t")
    p = p[(p.bin_kb == 500) & (p.direction == "loss") & p.feature_set.isin(SETS)].sort_values("expected_full")
    fig, ax = plt.subplots(figsize=(12, 6.25))
    for i, r in enumerate(p.itertuples()):
        ax.plot([r.min_fold_80, r.min_fold_80_global], [i, i], color=RULE, lw=3, zorder=0)
        ax.scatter(r.min_fold_80, i, s=34, color=MUTED, zorder=2, label="detectable, α = 0.05" if i == 0 else None)
        ax.scatter(r.min_fold_80_global, i, s=34, color=INK, zorder=2, label="detectable, global α" if i == 0 else None)
        ax.scatter(r.fold_pilot, i, s=60, marker="D", facecolor=BG, edgecolor=LOSS, linewidth=1.5, zorder=3,
                   label="pilot fold" if i == 0 else None)
        ax.text(4.55, i, f"{r.expected_full:.0f}", va="center", fontsize=9, color=MUTED)
    ax.set_yticks(range(len(p)), p.feature_set)
    ax.axvline(1, color=MUTED, lw=.8, ls=(0, (3, 3))); ax.set_xlim(0, 4.5); ax.grid(axis="x")
    ax.text(4.55, len(p) - .2, "E (49 donors)", fontsize=9, color=MUTED)
    ax.set_xlabel("fold (smallest detectable with 80% power at 49 donors; pilot estimate)")
    ax.set_title("Projected power · 500 kb losses · sets ordered by expected events")
    ax.legend(frameon=False, fontsize=10, loc="lower right")
    fig.tight_layout()
    return fig


FIGS = {"sum_test_a": sum_test_a, "sum_top": sum_top, "sum_concordance": sum_concordance, "sum_power": sum_power}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    plt = setup()
    os.makedirs(a.out, exist_ok=True)
    for name, f in FIGS.items():
        fig = f(plt, a.res)
        fig.savefig(f"{a.out}/{name}.png")
        plt.close(fig)
        print("wrote", name)


if __name__ == "__main__":
    main()

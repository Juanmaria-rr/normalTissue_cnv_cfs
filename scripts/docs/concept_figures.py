#!/usr/bin/env python3
"""Figures that explain replication timing / rt_rho and noise_excess with pilot data.

Used by the technical report (PNG files) and by notebooks/01_calibration_exploration.ipynb
(the notebook imports the functions and shows the figures inline).

  fig_rt_regions      RT along chr3:50-70 Mb (FHIT / FRA3B) at 500 kb, and mean RT of known
                      early and late regions
  fig_sphase          a theoretical mid-S cell on the same window vs a real high-rt_rho cell;
                      mean bin ratio per RT decile for three real cells
  fig_noise_toy       simulated bins (Poisson vs overdispersed) and the distribution of the
                      standardised neighbour differences d
  fig_noise_depth     P18 at 200 kb: raw neighbour SD vs noise_excess against depth
  fig_effective_reads depth_ok cells at 500 kb: reads vs effective reads (reads / noise_excess^2)

Inputs: results/pilot_cn/qc, results/pilot_cn/notebook (example-cell profiles), the RT
annotation of sWGS_relativeCN.
Usage: python3 scripts/docs/concept_figures.py --out docs/technical_report
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
QC = f"{ROOT}/results/pilot_cn/qc"
NB = f"{ROOT}/results/pilot_cn/notebook"
RT_FILE = f"{ROOT}/external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv"
DONORS = ["P18", "P20", "P36", "P07", "P08"]
DCOL = {"P18": "#4c72b0", "P20": "#55a868", "P36": "#8172b2", "P07": "#c44e52", "P08": "#dd8452"}
EARLY, LATE, THR = "#d08c2e", "#3f6f9e", "#b8962e"
S_CELL, G1_CELL, BIAS_CELL = "P07_chip1_co_r65c30", "P18_chip1_dna_r13c3", "P20_chip1_dna_r55c28"
plt.rcParams.update({"figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False,
                     "font.size": 9, "axes.titlesize": 9, "axes.titlelocation": "left"})

_rt = None


def rt50():
    global _rt
    if _rt is None:
        _rt = pd.read_csv(RT_FILE, sep="\t", usecols=["seqnames", "start", "end", "gc", "replicationTime"])
        _rt["seqnames"] = _rt.seqnames.astype(str)
    return _rt


def rt_binned(b):
    r = rt50().assign(bstart=lambda x: ((x.start - 1) // (b * 1000)) * b * 1000 + 1)
    return r.groupby(["seqnames", "bstart"]).replicationTime.mean()


def profile(cell, b):
    p = pd.read_csv(f"{NB}/profiles_{b}kb_bins.tsv.gz", sep="\t", dtype={"chromosome": str})
    p = p[(p.cell_id == cell) & p.use & p.chromosome.isin([str(i) for i in range(1, 23)])].copy()
    p["ratio"] = p.copynumber / np.nanmedian(p.segmented)
    p["rt"] = rt_binned(b).reindex(list(zip(p.chromosome, p.start))).values
    return p


def fig_rt_regions():
    r = rt_binned(500)
    w = r.loc["3"]
    w = w[(w.index >= 50e6) & (w.index < 70e6)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 3), gridspec_kw={"width_ratios": [2.2, 1]})
    x = (w.index + 250e3) / 1e6
    a.plot(x, w.values, color="#333333", lw=1.2)
    a.fill_between(x, w.values, 47, where=w.values >= 47, color=EARLY, alpha=.25, lw=0)
    a.fill_between(x, w.values, 47, where=w.values < 47, color=LATE, alpha=.25, lw=0)
    a.axvspan(59.735, 61.251, color="#999999", alpha=.18, lw=0)
    a.text(60.5, 80, "FHIT\n(FRA3B)", ha="center", va="top", fontsize=8)
    a.axhline(47, color="#999999", lw=.6, ls=":")
    a.set_ylim(10, 85); a.set_xlabel("chr3 (Mb, hg19)"); a.set_ylabel("RT (high = early)")
    a.set_title("A · Replication timing in 500 kb bins around FRA3B")
    regs = [("GAPDH", "12", 6.4e6, 6.8e6), ("TP53", "17", 7.4e6, 7.7e6), ("FHIT", "3", 59.7e6, 61.3e6),
            ("WWOX", "16", 78.1e6, 79.3e6), ("PARK2", "6", 161.7e6, 163.2e6)]
    r50 = rt50()
    vals = [r50[(r50.seqnames == c) & (r50.start >= s) & (r50.end <= e)].replicationTime.mean() for _, c, s, e in regs]
    cols = [EARLY if v > 47 else LATE for v in vals]
    b.barh([g for g, *_ in regs][::-1], vals[::-1], color=cols[::-1], height=.6)
    for i, v in enumerate(vals[::-1]):
        b.text(v + 1, i, f"{v:.0f}", va="center", fontsize=8)
    b.axvline(47, color="#999999", lw=.6, ls=":"); b.set_xlim(0, 90); b.set_xlabel("mean RT")
    b.set_title("B · Known regions")
    fig.tight_layout()
    return fig


def fig_sphase():
    rt = rt_binned(500).loc["3"]
    rt = rt[(rt.index >= 50e6) & (rt.index < 70e6)]
    x = (rt.index + 250e3) / 1e6
    med = rt50().replicationTime.median()
    copies = np.where(rt.values > med, 4, 2)
    sim = copies / 3
    real = profile(S_CELL, 500)
    rw = real[(real.chromosome == "3") & (real.start >= 50e6) & (real.start < 70e6)]
    fig = plt.figure(figsize=(11, 5.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.6, 1.6, 1])
    a = fig.add_subplot(gs[0, :2])
    a.step(x, sim, where="mid", color="#333333", lw=1.4, label="theoretical cell in mid S phase")
    a.scatter((rw.start + 250e3) / 1e6, rw.ratio, s=8, color=DCOL["P07"], label=f"{S_CELL} (bins)")
    for y in (0.7, 1.3):
        a.axhline(y, color=THR, ls="--", lw=.8)
    a.axvspan(59.735, 61.251, color="#999999", alpha=.18, lw=0)
    a.text(60.5, 1.62, "FHIT", ha="center", fontsize=8)
    a.set_ylim(0.3, 1.7); a.set_ylabel("ratio relative to the cell"); a.set_xlabel("chr3 (Mb)")
    a.legend(fontsize=7, loc="lower left")
    a.set_title("A · An S-phase cell \"gains\" early and \"loses\" late regions: false breakend at the edge of FRA3B")
    b = fig.add_subplot(gs[1, :2])
    b.scatter(rt.values, sim + np.random.default_rng(1).normal(0, .02, len(sim)), s=10, color="#333333", label="theoretical")
    b.scatter(rw.rt, rw.ratio, s=10, color=DCOL["P07"], label="real")
    b.axvline(med, color="#999999", ls=":", lw=.8)
    b.set_xlabel("bin RT"); b.set_ylabel("ratio"); b.legend(fontsize=7)
    b.set_title("B · The same stretch, ratio vs RT")
    c = fig.add_subplot(gs[:, 2])
    for cell, col, lab in ((S_CELL, DCOL["P07"], "P07 co · S phase"), (BIAS_CELL, DCOL["P20"], "P20 · batch bias"),
                           (G1_CELL, DCOL["P18"], "P18 · G1")):
        p = profile(cell, 500).dropna(subset=["rt", "ratio"])
        m = p.groupby(pd.qcut(p.rt, 10, labels=False)).ratio.mean()
        rho = p[["ratio", "rt"]].corr("spearman").iloc[0, 1]
        c.plot(m.index + 1, m.values, "o-", ms=3, color=col, label=f"{lab} (ρ = {rho:.2f})")
    c.axhline(1, color="#999999", lw=.6); c.set_ylim(.85, 1.2)
    c.set_xlabel("RT decile (1 late … 10 early)"); c.set_ylabel("mean ratio")
    c.legend(fontsize=7, loc="upper left"); c.set_title("C · Decile signature")
    fig.tight_layout()
    return fig


def fig_noise_toy(lam=200, n=2000, seed=1):
    # lam = 200: with small integer counts the MAD is coarse (Poisson at lam = 50 reads 0.89);
    # real corrected counts are continuous, so this does not apply to the pilot cells
    rng = np.random.default_rng(seed)
    pois = rng.poisson(lam, n)
    phi = 2.25  # variance = phi * lambda: negative binomial with size = lam / (phi - 1)
    size = lam / (phi - 1)
    od = rng.negative_binomial(size, size / (size + lam), n)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3), gridspec_kw={"width_ratios": [1.4, 1.4, 1.2]})
    for ax, x, col, lab in ((axes[0], pois, "#4c72b0", "Poisson (variance = λ)"),
                            (axes[1], od, "#c44e52", "overdispersed (variance = 2.25 λ)")):
        ax.scatter(np.arange(200), x[:200], s=5, color=col)
        ax.axhline(lam, color="#333333", lw=.8)
        ax.set_ylim(0, 2.2 * lam); ax.set_xlabel("bin"); ax.set_ylabel("reads")
        ax.set_title(f"{'A' if ax is axes[0] else 'B'} · {lab}")
    bins = np.linspace(-6, 6, 61)
    for x, col, lab in ((pois, "#4c72b0", "Poisson"), (od, "#c44e52", "overdispersed")):
        d = np.diff(x) / np.sqrt(2 * lam)
        ne = 1.4826 * np.median(np.abs(d - np.median(d)))
        axes[2].hist(d, bins=bins, density=True, histtype="step", color=col, lw=1.4, label=f"{lab}: noise_excess = {ne:.2f}")
    g = np.linspace(-6, 6, 200)
    axes[2].plot(g, np.exp(-g ** 2 / 2) / np.sqrt(2 * np.pi), color="#999999", ls=":", label="N(0, 1)")
    axes[2].set_xlabel("d = (x[i+1] − x[i]) / √(s[i+1] + s[i])"); axes[2].legend(fontsize=7)
    axes[2].set_title("C · Standardised differences")
    fig.tight_layout()
    return fig


def fig_noise_depth():
    q = pd.read_csv(f"{QC}/cn_qc_200kb.tsv", sep="\t")
    p = q[(q.donor == "P18") & (q.reads > 20000)].copy()
    lam = p.reads / p.used_bins
    p["raw_sd"] = np.sqrt(2 / lam) * p.noise_excess
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 3.2))
    a.scatter(p.reads / 1e3, p.raw_sd, s=4, alpha=.5, color=DCOL["P18"])
    xs = np.linspace(p.reads.min(), p.reads.max(), 200)
    a.plot(xs / 1e3, np.sqrt(2 / (xs / p.used_bins.iloc[0])), color="#333333", lw=1, label="pure Poisson: √(2/λ)")
    a.set_xscale("log"); a.set_xlabel("usable reads (thousands, log)"); a.set_ylabel("raw neighbour SD (ratio)")
    a.legend(fontsize=7); a.set_title("A · Raw SD falls with depth")
    b.scatter(p.reads / 1e3, p.noise_excess, s=4, alpha=.5, color=DCOL["P18"])
    b.axhline(1, color="#333333", lw=1)
    b.set_xscale("log"); b.set_ylim(.8, 1.6); b.set_xlabel("usable reads (thousands, log)"); b.set_ylabel("noise_excess")
    b.set_title("B · noise_excess stays close to 1 at any depth")
    for cell, lab in (("P18_chip1_dna_r0c50", "r0c50"), ("P18_chip1_dna_r36c10", "r36c10")):
        r = p[p.cell_id == cell]
        if len(r):
            a.annotate(lab, (r.reads.iloc[0] / 1e3, r.raw_sd.iloc[0]), xytext=(6, 6), textcoords="offset points", fontsize=7)
            b.annotate(lab, (r.reads.iloc[0] / 1e3, r.noise_excess.iloc[0]), xytext=(6, 6), textcoords="offset points", fontsize=7)
    fig.tight_layout()
    return fig


def fig_effective_reads():
    q = pd.read_csv(f"{QC}/cn_qc_500kb.tsv", sep="\t")
    ok = q[q.depth_ok == 1].copy()
    ok["eff"] = ok.reads / ok.noise_excess ** 2
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for d in DONORS:
        g = ok[ok.donor == d]
        ax.scatter(g.reads / 1e3, g.eff / 1e3, s=3, alpha=.4, color=DCOL[d], label=d)
    lim = [ok.reads.min() / 1e3 * .9, ok.reads.max() / 1e3 * 1.1]
    ax.plot(lim, lim, color="#333333", lw=.8)
    for v, lab in ((145.286, "500 kb threshold"), (64.022, "1 Mb threshold")):
        ax.axhline(v, color=THR, ls="--", lw=.8); ax.text(lim[1], v, f" {lab}", va="center", fontsize=7)
    for cell in ("P18_chip1_dna_r13c3", "P07_chip1_co_r65c30", "P07_chip1_co_r30c35"):
        r = ok[ok.cell_id == cell]
        if len(r):
            ax.scatter(r.reads / 1e3, r.eff / 1e3, s=40, facecolors="none", edgecolors="k")
            ax.annotate(cell.replace("_chip1_", " "), (r.reads.iloc[0] / 1e3, r.eff.iloc[0] / 1e3), xytext=(8, -2),
                        textcoords="offset points", fontsize=7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("usable reads (thousands)")
    ax.set_ylabel("effective reads = reads / noise_excess² (thousands)")
    ax.legend(fontsize=7, markerscale=3, loc="upper left")
    ax.set_title("Cells with depth_ok at 500 kb: how much useful information they carry")
    fig.tight_layout()
    return fig


FIGS = {"fig_rt_regions": fig_rt_regions, "fig_sphase": fig_sphase, "fig_noise_toy": fig_noise_toy,
        "fig_noise_depth": fig_noise_depth, "fig_effective_reads": fig_effective_reads}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    plt.switch_backend("Agg")  # set here, not at import, so the notebook keeps inline output
    for name, f in FIGS.items():
        fig = f()
        fig.savefig(os.path.join(a.out, f"{name}.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("wrote", name)


if __name__ == "__main__":
    main()

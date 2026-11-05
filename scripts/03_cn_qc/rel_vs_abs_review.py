#!/usr/bin/env python3
"""Review PDFs (English) of relative-ratio vs integer-CN breakends, and a curated example set.

Modes:
  --select    N random cells per donor among cells with depth_ok at 500 kb (seed 1)
              -> <out>/review_cells.txt
  --examples  curated "good example" breakend cells from the whole pilot
              -> <out>/example_cells.txt, <out>/example_cells.tsv (cell, donor, event, why)
  --report    render the cells of --cells-file into PDFs (one page per cell, 3 bins) and
              write the concordance tables. Export the cells first with export_relcn.R
              (--cells <file> --out <out>/profiles_<bin>kb).

Breakends, per cell and bin, by two routes on the same segmentation (autosomes, use bins):
  relative  rel = segment / median of the cell's segments; gain > 1.3, loss < 0.7
            (call_breakends.R rule)
  absolute  best unfiltered fit (lowest clonality over ploidy x purity in
            <abscn>/<donor>_<bin>kb/<cell>_resfitting.txt), converted with the
            fitting_absCN formula CN = (seg / seqdepth - 2(1 - p)) / p,
            seqdepth = mean(copynumber) / (P p + 2 (1 - p)), and rounded; any change of
            integer CN is a breakend
  In both, consecutive bins with the same state are merged, segments < 1 Mb dropped, and a
  breakend is a state change between consecutive retained segments of a chromosome. A
  breakend is "shared" if the other route has one on the same chromosome within 2 bins.

Example selection (--examples), from the call_breakends.R outputs:
  cells with depth_ok at 500 kb, noise_excess <= P18 p99 at 500 kb, rt_rho <= 0.10;
  a 500 kb breakend is "reproduced" if the same cell has a breakend on the same chromosome
  within 3 Mb at 1 Mb (and at 200 kb when the cell has depth_ok there);
  events = reproduced breakends clustered per donor and chromosome within 3 Mb, so a
  clone contributes one event; per event the cell with the lowest noise_excess is kept;
  a breakend whose interval [lo, hi] (which spans any gap of filtered bins) overlaps an
  approximate pericentromeric block +/- 1 Mb (centromere + adjacent heterochromatin, e.g.
  chr1 121.5-142.6 Mb) is labelled arm-level, others interstitial;
  per donor up to --per-donor interstitial events plus one arm-level event, preferring
  events seen in more cells.

Limitations: mean(copynumber) uses exported autosomal use bins (may differ slightly from
the upstream object); rounding flips state at x.5; pericentromeric blocks are approximate.

Usage:
  python3 scripts/03_cn_qc/rel_vs_abs_review.py --select   --qc results/pilot_cn/qc --out results/pilot_cn/review
  python3 scripts/03_cn_qc/rel_vs_abs_review.py --examples --qc results/pilot_cn/qc --out results/pilot_cn/review
  python3 scripts/03_cn_qc/rel_vs_abs_review.py --report --qc results/pilot_cn/qc --abscn results/pilot_cn/absCN \
      --out results/pilot_cn/review --cells-file results/pilot_cn/review/review_cells.txt --prefix review
  python3 scripts/03_cn_qc/rel_vs_abs_review.py --report --qc results/pilot_cn/qc --abscn results/pilot_cn/absCN \
      --out results/pilot_cn/review --cells-file results/pilot_cn/review/example_cells.txt --prefix example_breakends \
      --single-pdf --profiles examples_profiles
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

BINS = [200, 500, 1000]
GAIN, LOSS, MIN_MB, MATCH_BINS = 1.3, 0.7, 1.0, 2
# approximate hg19 pericentromeric blocks (centromere + adjacent heterochromatin; for
# acrocentrics the whole p arm), Mb, from UCSC cytoBand acen/gvar/stalk bands
PERICEN = {1: (121.5, 142.6), 2: (90.5, 96.0), 3: (87.9, 93.9), 4: (48.2, 52.7), 5: (46.1, 50.7),
           6: (58.5, 62.6), 7: (57.9, 61.9), 8: (43.1, 48.1), 9: (38.0, 71.0), 10: (38.0, 42.9),
           11: (51.0, 55.8), 12: (33.3, 38.2), 13: (0, 19.5), 14: (0, 19.1), 15: (0, 20.7),
           16: (34.6, 46.5), 17: (22.2, 25.8), 18: (15.0, 18.9), 19: (24.2, 28.6), 20: (25.6, 29.8),
           21: (0, 14.3), 22: (0, 16.3)}
C_SHARED, C_REL, C_ABS = "#2e8b57", "#c0392b", "#2c6fbb"


def select(qc, out, n, seed):
    q = pd.read_csv(f"{qc}/cn_qc_500kb.tsv", sep="\t")
    ok = q[q.depth_ok == 1]
    cells = pd.concat([g.sample(min(n, len(g)), random_state=seed) for _, g in ok.groupby("donor")])
    os.makedirs(out, exist_ok=True)
    cells.cell_id.to_csv(f"{out}/review_cells.txt", index=False, header=False)
    print(len(cells), "cells;", cells.donor.value_counts().to_dict())


def examples(qc, out, per_donor):
    q = {b: pd.read_csv(f"{qc}/cn_qc_{b}kb.tsv", sep="\t").set_index("cell_id") for b in BINS}
    be = {b: pd.read_csv(f"{qc}/breakends_{b}kb_breakends.tsv", sep="\t", dtype={"chromosome": str}) for b in BINS}
    q5 = q[500]
    n99 = q5[(q5.donor == "P18") & (q5.depth_ok == 1)].noise_excess.quantile(.99)
    good = q5[(q5.depth_ok == 1) & (q5.noise_excess <= n99) & (q5.rt_rho <= 0.10)].index
    b5 = be[500][be[500].cell_id.isin(good)].copy()

    def near(b, cell, chrom, pos):
        x = be[b][(be[b].cell_id == cell) & (be[b].chromosome == chrom)]
        return bool((np.abs(x.pos - pos) <= 3e6).any())

    keep = []
    for r in b5.itertuples():
        ok1 = near(1000, r.cell_id, r.chromosome, r.pos)
        ok2 = near(200, r.cell_id, r.chromosome, r.pos) if q[200].loc[r.cell_id].depth_ok == 1 else True
        keep.append(ok1 and ok2)
    b5 = b5[keep].copy()
    b5["donor"] = b5.cell_id.str.split("_").str[0]
    b5["chrn"] = b5.chromosome.astype(int)
    pc = b5.chrn.map(PERICEN)
    # the breakend interval [lo, hi] spans any gap of filtered bins, so a breakend whose flanking
    # segments sit on opposite sides of the pericentromeric block overlaps it
    b5["arm_level"] = [(l / 1e6 <= bh + 1) and (h / 1e6 >= bl - 1) for l, h, (bl, bh) in zip(b5.lo, b5.hi, pc)]
    b5["noise"] = b5.cell_id.map(q5.noise_excess)
    b5 = b5.sort_values(["donor", "chrn", "pos"])
    ev, k = [], 0
    for (_, _), g in b5.groupby(["donor", "chrn"]):
        last = None
        for r in g.itertuples():
            if last is None or r.pos - last > 3e6:
                k += 1
            ev.append(k); last = r.pos
    b5["event"] = ev
    rows = []
    for e, g in b5.groupby("event"):
        best = g.sort_values("noise").iloc[0]
        rows.append(dict(cell_id=best.cell_id, donor=best.donor, chromosome=best.chromosome, pos_mb=round(best.pos / 1e6, 1),
                         change=f"{best.left_state}->{best.right_state}", arm_level=bool(best.arm_level),
                         cells_with_event=g.cell_id.nunique(), noise_excess=best.noise))
    ev = pd.DataFrame(rows).sort_values(["donor", "cells_with_event"], ascending=[True, False])
    inter = ev[~ev.arm_level].groupby("donor").head(per_donor)
    arm = ev[ev.arm_level].groupby("donor").head(1)
    pick = pd.concat([inter, arm]).drop_duplicates("cell_id").sort_values(["donor", "arm_level"]).copy()
    pick["why"] = np.where(pick.arm_level, "arm-level", "interstitial") + ", reproduced at 1 Mb (and 200 kb if deep); " \
        + pick.cells_with_event.astype(str) + " cell(s) of the donor carry it"
    os.makedirs(out, exist_ok=True)
    pick.to_csv(f"{out}/example_cells.tsv", sep="\t", index=False)
    pick.cell_id.to_csv(f"{out}/example_cells.txt", index=False, header=False)
    print(f"noise cut {n99:.3f}; {len(ev)} candidate events; {len(pick)} examples")
    print(pick[["cell_id", "chromosome", "pos_mb", "change", "arm_level", "cells_with_event"]].to_string(index=False))


def best_fit(abscn, cell, b):
    f = f"{abscn}/{cell.split('_')[0]}_{b}kb/{cell}_resfitting.txt"
    if not os.path.exists(f):
        return None
    r = pd.read_csv(f, sep="\t")
    r = r.loc[r.clonality.idxmin()]
    return float(r.ploidy), float(r.purity)


def segments(p, state):
    brk = np.r_[True, (p.chrn.values[1:] != p.chrn.values[:-1]) | (state[1:] != state[:-1])]
    s = pd.DataFrame({"chrn": p.chrn.values, "start": p.start.values, "end": p.end.values,
                      "state": state, "g": np.cumsum(brk), "i": np.arange(len(p))})
    s = s.groupby("g").agg(chrn=("chrn", "first"), start=("start", "min"), end=("end", "max"),
                           state=("state", "first"), i0=("i", "min"))
    return s[(s.end - s.start + 1) / 1e6 >= MIN_MB].reset_index(drop=True)


def breakends(s):
    idx = np.flatnonzero((s.chrn.values[1:] == s.chrn.values[:-1]) & (s.state.values[1:] != s.state.values[:-1]))
    return pd.DataFrame({"chrn": s.chrn.values[idx], "i": s.i0.values[idx + 1], "pos": s.start.values[idx + 1]})


def calls(p, ploidy, purity):
    seg = p.segmented.values
    rel = seg / np.nanmedian(seg)
    rstate = np.where(rel > GAIN, 1, np.where(rel < LOSS, -1, 0))
    seqdepth = np.nanmean(p.copynumber.values) / (ploidy * purity + 2 * (1 - purity))
    cn = (seg / seqdepth - 2 * (1 - purity)) / purity
    icn = np.round(cn)
    rb, ab = breakends(segments(p, rstate)), breakends(segments(p, icn.astype(int)))

    def match(a, b):
        return np.array([((b.chrn == r.chrn) & (np.abs(b.i - r.i) <= MATCH_BINS)).any() for r in a.itertuples()], bool)
    return rel, cn, icn, rb, ab, match(rb, ab), match(ab, rb)


def load_profiles(out, b, cells, pref="profiles"):
    p = pd.read_csv(f"{out}/{pref}_{b}kb_bins.tsv.gz", sep="\t", dtype={"chromosome": str})
    p = p[p.cell_id.isin(cells) & p.use & p.chromosome.isin([str(i) for i in range(1, 23)])]
    p = p[np.isfinite(p.segmented) & np.isfinite(p.copynumber)].copy()
    p["chrn"] = p.chromosome.astype(int)
    return {c: g.sort_values(["chrn", "start"]).reset_index(drop=True) for c, g in p.groupby("cell_id")}


def chrom_axis(ax, p):
    edges = np.r_[0, np.flatnonzero(np.diff(p.chrn.values)) + 1, len(p)]
    for e in edges[1:-1]:
        ax.axvline(e - .5, color="#e2e2e2", lw=.5)
    mids = (edges[:-1] + edges[1:]) / 2
    ax.set_xticks(mids, [str(c) for c in p.chrn.values[edges[:-1]]], fontsize=6)
    ax.tick_params(axis="x", length=0)


def be_label(p, r):
    return f"chr{r.chrn}:{p.start.values[r.i] / 1e6:.1f}"


def cell_page(c, prof, qcs, abscn, note=None):
    fig, axes = plt.subplots(6, 1, figsize=(11.7, 8.3), gridspec_kw={"height_ratios": [1.3, 1] * 3})
    rows = []
    for k, b in enumerate(BINS):
        ar, ai = axes[2 * k], axes[2 * k + 1]
        p, fit = prof[b].get(c), best_fit(abscn, c, b)
        if p is None or fit is None:
            ar.text(.5, .5, f"{b} kb: no data", transform=ar.transAxes, ha="center"); ar.axis("off"); ai.axis("off"); continue
        P, pur = fit
        rel, cn, icn, rb, ab, rs, as_ = calls(p, P, pur)
        q = qcs[b].loc[c]
        rows.append(dict(cell_id=c, donor=c.split("_")[0], bin_kb=b, depth_ok=int(q.depth_ok), noise_excess=q.noise_excess,
                         rt_rho=q.rt_rho, ploidy_best=P, purity_best=pur, n_rel=len(rb), n_abs=len(ab),
                         rel_shared=int(rs.sum()), abs_shared=int(as_.sum()),
                         rel_breakends=";".join(be_label(p, r) for r in rb.itertuples())))
        x = np.arange(len(p))
        ar.scatter(x, np.clip(p.copynumber / np.nanmedian(p.segmented), 0, 3), s=1, c="#a7b0b8", lw=0)
        ar.plot(x, np.clip(rel, 0, 3), color="#222222", lw=.9)
        for y in (LOSS, GAIN):
            ar.axhline(y, color="#b8962e", lw=.6, ls="--")
        ai.scatter(x, np.clip(cn, 0, 6), s=1, c="#a7b0b8", lw=0)
        ai.plot(x, np.clip(icn, 0, 6), color="#222222", lw=.9)
        ai.axhline(round(P), color="#999999", lw=.6, ls=":")
        for ax in (ar, ai):
            chrom_axis(ax, p); ax.set_xlim(-.5, len(p) - .5)
        ar.set_xticklabels([])
        # band of +/- 1 bin around the boundary, drawn under the data so the segment step stays visible
        w = max(1.0, len(p) / 400)
        for r, sh in zip(rb.itertuples(), rs):
            ar.axvspan(r.i - .5 - w, r.i - .5 + w, color=C_SHARED if sh else C_REL, alpha=.35, lw=0, zorder=0)
        for r, sh in zip(ab.itertuples(), as_):
            ai.axvspan(r.i - .5 - w, r.i - .5 + w, color=C_SHARED if sh else C_ABS, alpha=.35, lw=0, zorder=0)
        ar.set_ylim(0, 3); ai.set_ylim(0, 6)
        ar.set_ylabel("ratio", fontsize=7); ai.set_ylabel("integer CN", fontsize=7)
        lab = ", ".join(be_label(p, r) + " Mb" for r in rb.itertuples()) or "none"
        ar.set_title(f"{b} kb · {int(q.reads / 1e3)}k usable reads · depth_ok {int(q.depth_ok)} · noise_excess {q.noise_excess:.2f} · "
                     f"rt_rho {q.rt_rho:.2f} · relative: {len(rb)} breakend(s) [{lab}]", fontsize=7, loc="left")
        ai.set_title(f"absolute (best unfiltered fit): ploidy {P:.1f}, purity {pur:.2f} · {len(ab)} breakend(s), "
                     f"{as_.sum()} shared with relative", fontsize=7, loc="left")
    fig.text(.01, .99, c, fontsize=10, fontweight="bold", va="top")
    fig.text(.01, .965, (note + "   ·   " if note else "") +
             "breakend bands: green = both routes · red = relative only · blue = absolute only", fontsize=7.5, va="top")
    fig.tight_layout(rect=(0, 0, 1, .945))
    return fig, rows


def legend_page(title, lines):
    fig = plt.figure(figsize=(11.7, 8.3))
    fig.text(.06, .92, title, fontsize=15, fontweight="bold")
    y = .85
    for ln in lines:
        fig.text(.06, y, ln, fontsize=9, va="top")
        y -= .035 * (1 + ln.count("\n"))
    return fig


def table_page(t):
    fig = plt.figure(figsize=(11.7, 8.3))
    fig.text(.06, .93, "Selected examples", fontsize=15, fontweight="bold")
    fig.text(.06, .895, "Cells with depth at 500 kb, noise_excess <= P18 p99, rt_rho <= 0.10; breakend reproduced at 1 Mb "
             "(and 200 kb when deep enough); one cell per event (clones counted once).\nArm-level = breakend interval "
             "overlaps the pericentromeric block; interstitial otherwise. 'Cells' = cells of the donor carrying the event.",
             fontsize=8, va="top")
    cols = ["page", "cell", "donor", "breakend (500 kb)", "change", "type", "cells"]
    rows = [[str(i + 3), r.cell_id, r.donor, f"chr{r.chromosome}:{r.pos_mb} Mb", r.change.replace("->", " → "),
             "arm-level" if r.arm_level else "interstitial", str(r.cells_with_event)] for i, r in enumerate(t.itertuples())]
    ax = fig.add_axes([.06, .05, .88, .78]); ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=cols, loc="upper left", cellLoc="left", colLoc="left",
                  colWidths=[.06, .24, .08, .2, .18, .12, .07])
    tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.25)
    return fig


METHOD = [
    "Data: Lin et al. 2024 (PRJNA1003661), single-cell DNA-seq of normal breast epithelium; pilot donors P18, P20, P36, P07, P08.",
    "Reads: R1 only, cutadapt (Tn5 mosaic end), bwa mem to hg19, duplicates marked; usable = MAPQ >= 37, non-duplicate.",
    "Relative CN: QDNAseq (GC + mappability correction, CBS segmentation, undo.SD = 1) at 200 kb, 500 kb and 1 Mb bins.",
    "Each bin block, top panel: ratio = bin (grey) or segment (black) value / median of the cell's segments.\n"
    "Dashed lines = 0.7 / 1.3, the gain/loss thresholds of the paper. Retained segments >= 1 Mb.",
    "Bottom panel: integer copy number from the best unfiltered absolute fit (fitting_absCN grid, ploidy 1.4-8, purity 0.9-1),\n"
    "rounded; dotted line = rounded ploidy.",
    "Breakend = change of state between consecutive retained segments of a chromosome. Bands: green = found by both routes\n"
    "(within 2 bins), red = relative only, blue = absolute (integer) only. Chromosome numbers on the x axis.",
    "depth_ok = usable reads >= 15 x 2 x usable bins / 0.91 (diploid threshold of fitting_absCN): 396k / 145k / 64k reads.",
    "noise_excess = robust SD of neighbour-bin differences divided by the Poisson expectation (1 = pure counting noise).",
    "rt_rho = Spearman correlation between bin relative CN and replication timing (> 0 suggests an S-phase cell).",
    "Thresholds are not yet calibrated; calls are preliminary.",
]


def report(qc, abscn, out, cells_file, prefix, notes=None, single_pdf=False, profiles="profiles", table=None):
    cells = pd.read_csv(cells_file, header=None)[0].tolist()
    prof = {b: load_profiles(out, b, set(cells), profiles) for b in BINS}
    qcs = {b: pd.read_csv(f"{qc}/cn_qc_{b}kb.tsv", sep="\t").set_index("cell_id") for b in BINS}
    rows, figs = [], {}
    for c in cells:
        fig, r = cell_page(c, prof, qcs, abscn, (notes or {}).get(c))
        rows += r
        figs.setdefault("all" if single_pdf else c.split("_")[0], []).append(fig)
    for key, fl in figs.items():
        name = f"{out}/{prefix}.pdf" if single_pdf else f"{out}/{prefix}_{key}.pdf"
        title = "Example breakends, pilot cells" if single_pdf else f"Random review cells, donor {key} (n = {len(fl)})"
        with PdfPages(name) as pdf:
            lp = legend_page(title, METHOD); pdf.savefig(lp); plt.close(lp)
            if single_pdf and table is not None:
                tp = table_page(table); pdf.savefig(tp); plt.close(tp)
            for f in fl:
                pdf.savefig(f); plt.close(f)
        print("wrote", name)
    t = pd.DataFrame(rows)
    t.to_csv(f"{out}/{prefix}_concordance.tsv", sep="\t", index=False)
    s = t[t.depth_ok == 1].groupby(["bin_kb", "donor"]).agg(
        cells=("cell_id", "size"), rel=("n_rel", "sum"), abs=("n_abs", "sum"),
        rel_shared=("rel_shared", "sum"), abs_shared=("abs_shared", "sum")).reset_index()
    s.to_csv(f"{out}/{prefix}_concordance_summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--select", action="store_true")
    ap.add_argument("--examples", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--qc", required=True)
    ap.add_argument("--abscn")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells-file")
    ap.add_argument("--prefix", default="review")
    ap.add_argument("--single-pdf", action="store_true")
    ap.add_argument("--profiles", default="profiles", help="prefix of the export_relcn.R outputs in --out")
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--per-donor", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.select:
        select(a.qc, a.out, a.n, a.seed)
    if a.examples:
        examples(a.qc, a.out, a.per_donor)
    if a.report:
        notes = None
        ex = f"{a.out}/example_cells.tsv"
        if a.single_pdf and os.path.exists(ex):
            e = pd.read_csv(ex, sep="\t")
            e = e.set_index("cell_id").loc[pd.read_csv(a.cells_file, header=None)[0]].reset_index()
            notes = {r.cell_id: f"example: chr{r.chromosome} ~{r.pos_mb} Mb, {r.change.replace('->', ' → ')}, {r.why}" for r in e.itertuples()}
        report(a.qc, a.abscn, a.out, a.cells_file, a.prefix, notes, a.single_pdf, a.profiles, e if notes else None)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Data funnel of the pilot: how many cells (and breakends / events) survive each criterion.

Cell stages (per donor; stages 4-11 per bin):
   1 wells          wells of the donor in ENA (cells.tsv, one per cell_id)
   2 imaging        imaging qc_flag == Good, cell_flag == yes, not a control well (pilot selection)
   3 aligned        aligned and level-1 QC computed
   4 mapping        mapping rate >= 90% (alarm only; not used to exclude)
   5 fitted         >= 20,000 usable reads (pureSamples_absCN --total.reads)
   6 depth_ok       mean corrected reads per usable bin / 2 >= 15 (recalibrate_depth.py)
   7 noise          + noise_excess <= 1.5
   8 rt             + rt_rho <= 0.10  (= pass_cuts, the main analysis set)
   9 any breakend   + >= 1 breakend, main definition (E2-E4)
  10 interstitial   + >= 1 interstitial breakend (E2-E4)
  11 interst. loss  + >= 1 interstitial loss breakend (E2-E4)
Breakend / event stages (main set, per bin): breakends E1-E4 -> E2-E4 -> interstitial ->
collapsed events -> losses -> events overlapping MiDAS (Macheret union) -> overlapping
MiDAS HS68.
Each stage is a subset of the previous one (cumulative), except 4, which is reported
alongside 3 and does not remove cells.

Outputs: <res>/funnel_cells.tsv, <res>/funnel_breakends.tsv, and a figure in --fig-dir.

Usage: python3 scripts/06_enrichment/funnel.py --res results/enrichment/pilot --cells data/metadata/cells.tsv \
    --align data/pilot/align_qc_pilot.tsv --fig-dir docs/enrichment
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DONORS = ["P07", "P08", "P18", "P20", "P36"]
BINS = [200, 500, 1000]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--align", required=True)
    ap.add_argument("--fig-dir", required=True)
    a = ap.parse_args()

    meta = pd.read_csv(a.cells, sep="\t", low_memory=False)
    meta = meta[meta.donor.isin(DONORS)].drop_duplicates("cell_id")
    img = meta[(meta.qc_flag == "Good") & (meta.cell_flag == "yes") & (meta.control_well == 0)]
    al = pd.read_csv(a.align, sep="\t")
    mc = pd.read_csv(f"{a.res}/master_cells.tsv", sep="\t")

    rows = []
    for d in DONORS + ["total"]:
        sel = (lambda x: x) if d == "total" else (lambda x, d=d: x[x.donor == d])
        common = [("1 wells", len(sel(meta))), ("2 imaging", len(sel(img))), ("3 aligned", len(sel(al))),
                  ("4 mapping >= 90%", int((sel(al).mapped_pct >= 90).sum()))]
        for b in BINS:
            c = sel(mc[mc.bin_kb == b])
            st = c[c.depth_ok == 1]
            n1 = st[st.noise_excess <= 1.5]
            r1 = n1[n1.rt_rho <= 0.10]
            stages = common + [("5 fitted (>= 20k reads)", len(c)), ("6 depth_ok", len(st)),
                               ("7 noise_excess <= 1.5", len(n1)), ("8 rt_rho <= 0.10", len(r1)),
                               ("9 >= 1 breakend (E2-E4)", int((r1.breakends_main > 0).sum())),
                               ("10 >= 1 interstitial", int(((r1.interstitial_loss + r1.interstitial_gain) > 0).sum())),
                               ("11 >= 1 interstitial loss", int((r1.interstitial_loss > 0).sum()))]
            for s, n in stages:
                rows.append(dict(donor=d, bin_kb=b, stage=s, cells=n))
    fc = pd.DataFrame(rows)
    fc.to_csv(f"{a.res}/funnel_cells.tsv", sep="\t", index=False)

    brows = []
    be = pd.read_csv(f"{a.res}/master_breakends.tsv", sep="\t", dtype={"chromosome": str})
    ev = pd.read_csv(f"{a.res}/master_events.tsv", sep="\t", dtype={"chromosome": str})
    for b in BINS:
        x = be[(be.set == "main") & (be.bin_kb == b)]
        e = ev[(ev.set == "main") & (ev.bin_kb == b)]
        m = x[x.stratum.isin(["E2", "E3", "E4"])]
        i = m[m.type == "interstitial"]
        ei = e[(e.type == "interstitial") & e.main_def]
        el = ei[ei.direction == "loss"]
        for s, n in [("breakends E1-E4", len(x)), ("breakends E2-E4", len(m)), ("interstitial breakends", len(i)),
                     ("interstitial events (collapsed)", len(ei)), ("interstitial loss events", len(el)),
                     ("loss events in MiDAS (Macheret union)", int(el["ov:macheret2020:all"].sum())),
                     ("loss events in MiDAS HS68", int(el["ov:macheret2020:HS68"].sum()))]:
            brows.append(dict(bin_kb=b, stage=s, n=n))
    fb = pd.DataFrame(brows)
    fb.to_csv(f"{a.res}/funnel_breakends.tsv", sep="\t", index=False)

    tot = fc[fc.donor == "total"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
    for ax, b in zip(axes, BINS):
        t = tot[(tot.bin_kb == b) & (tot.stage != "4 mapping >= 90%")]
        y = np.arange(len(t))[::-1]
        ax.barh(y, t.cells, color="#0f6e73", height=.65)
        for k, (yy, v) in enumerate(zip(y, t.cells)):
            pct = "" if k == 0 else f" ({100 * v / t.cells.iloc[1]:.0f} %)"   # % of imaging-selected cells
            ax.text(v, yy, f" {v:,}".replace(",", ".") + pct, va="center", fontsize=7)
        ax.set_yticks(y, t.stage); ax.set_xscale("log"); ax.set_xlim(50, 4e4)
        ax.set_title("1 Mb" if b == 1000 else f"{b} kb", loc="left"); ax.set_xlabel("cells (log scale)")
    fig.tight_layout()
    os.makedirs(a.fig_dir, exist_ok=True)
    fig.savefig(f"{a.fig_dir}/fig_funnel.png", dpi=150, bbox_inches="tight")
    print(tot.pivot(index="stage", columns="bin_kb", values="cells").to_string())
    print(fb.pivot(index="stage", columns="bin_kb", values="n").to_string())


if __name__ == "__main__":
    main()

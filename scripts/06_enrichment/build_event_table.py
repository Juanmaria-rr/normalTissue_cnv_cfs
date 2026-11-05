#!/usr/bin/env python3
"""Breakend and event tables for the CFS enrichment tests, per bin.

Cell sets (one output folder each):
  main       depth_ok, noise_excess <= 1.5, rt_rho <= 0.10 (decisions 2026-11-04)
  nocuts     depth_ok only (sensitivity)
  samecells  main cuts passed in all three bins (sensitivity: the same cells in every bin)
  cells500   the main cells of 500 kb, called at every bin whatever their depth there, so
             bins differ in resolution only, not in which cells enter
  refined    main-set 500 kb breakends refined to 200 kb in the same cell: a 200 kb breakend
             of the same direction whose position lies within [lo - 500 kb, hi + 500 kb]
             (nearest one); interval from the 200 kb call, deviation and stratum from the
             500 kb call; breakends without a 200 kb match are left out (events_200kb.tsv only;
             refinement_summary.tsv gives the matched fraction)
depth_ok is the recalibrated rule (mean corrected reads per used bin / 2 >= 15, see
cn_qc_metrics.R / recalibrate_depth.py).
Breakends (per cell), from the run-length segments of call_breakends.R:
  - state per segment with the permissive rule: gain if ratio > 1.2, loss if < 0.8;
    consecutive segments with the same state are merged (ratio = length-weighted mean);
    merged segments with < 1 Mb of usable bins (n_bins x bin size) are dropped, so a
    segment cannot pass the length filter by spanning a gap of filtered bins
  - breakend = boundary between consecutive retained segments of a chromosome with
    different state; interval [lo, hi] = start of the last bin of the left segment to end
    of the first bin of the right one (one bin either side plus any gap)
  - deviation = max |ratio - 1| of the two flanks; stratum E1 [0.2, 0.3), E2 [0.3, 0.4),
    E3 [0.4, 0.6), E4 >= 0.6; main_def = E2-E4 (the paper's 1.3 / 0.7)
  - direction = loss / gain from the flank with the larger deviation; "mixed" when one
    flank is a gain and the other a loss
  - type = arm_level if [lo, hi] overlaps a pericentromeric block +/- 1 Mb, else interstitial
  - terminal = the altered flank (the one with the larger deviation) is the first or last
    retained segment of the chromosome, i.e. the alteration runs to the chromosome end
  - rt = mean replication timing of usable bins overlapping [lo, hi]; rt_quintile from
    the genome-wide distribution of usable-bin RT at that bin size
Events (per donor): breakends of the same chromosome, direction and type whose positions
  lie within LINK_BP (default 2 Mb) of the previous one (single linkage on sorted positions)
  form one event. Event pos / lo / hi = medians of its breakends; deviation = median;
  n_cells = distinct cells; terminal = majority of its breakends. A clone therefore counts
  once; two independent breaks closer than LINK_BP also merge. --collapse-on main links only
  E2-E4 breakends (E1 breakends then cannot chain events together).
Robustness options: --link-bp, --collapse-on all|main, --sets (subset of cell sets), --tag
  (suffix of the output folders, e.g. main_link1mb).

Outputs (<out>/<set><tag>/): breakends_<bin>kb.tsv (with event_id), events_<bin>kb.tsv,
usable_bins_rt_<bin>kb.tsv; <out>/summary<tag>.tsv (events per bin x donor x type x direction
x stratum).

Usage: python3 scripts/06_enrichment/build_event_table.py --qc results/pilot_cn/qc \
    --profiles results/pilot_cn/review --rt external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv \
    --out results/enrichment/pilot
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from genome_annot import interval_rt, overlaps_pericen, rt_per_bin, usable_bins  # noqa: E402

BINS = [200, 500, 1000]
T, MIN_BP, LINK_BP = 0.2, 1_000_000, 2_000_000
STRATA = [(0.2, 0.3, "E1"), (0.3, 0.4, "E2"), (0.4, 0.6, "E3"), (0.6, np.inf, "E4")]
CUTS = {"main": dict(noise=1.5, rt=0.10), "nocuts": dict(noise=np.inf, rt=np.inf),
        "samecells": dict(noise=1.5, rt=0.10), "cells500": dict(noise=1.5, rt=0.10)}
REFINE_BIN = 200


def stratum(d):
    for lo, hi, s in STRATA:
        if lo <= d < hi:
            return s
    return None


def cell_breakends(s, bin_bp):
    s = s.reset_index(drop=True)
    st = np.where(s.rel > 1 + T, 1, np.where(s.rel < 1 - T, -1, 0))
    ln = (s.end - s.start + 1).values
    nb = s.n_bins.values
    brk = np.r_[True, (s.chromosome.values[1:] != s.chromosome.values[:-1]) | (st[1:] != st[:-1])]
    g = np.cumsum(brk)
    m = pd.DataFrame({"g": g, "chromosome": s.chromosome.values, "start": s.start.values, "end": s.end.values,
                      "st": st, "w": ln, "nb": nb, "rw": s.rel.values * ln})
    m = m.groupby("g").agg(chromosome=("chromosome", "first"), start=("start", "min"), end=("end", "max"),
                           st=("st", "first"), w=("w", "sum"), nb=("nb", "sum"), rw=("rw", "sum")).reset_index(drop=True)
    m["rel"] = m.rw / m.w
    m = m[m.nb * bin_bp >= MIN_BP].reset_index(drop=True)
    first = m.chromosome.ne(m.chromosome.shift()).values
    last = m.chromosome.ne(m.chromosome.shift(-1)).values
    out = []
    for i in range(len(m) - 1):
        a, b = m.iloc[i], m.iloc[i + 1]
        if a.chromosome != b.chromosome or a.st == b.st:
            continue
        lo, hi = a.end - bin_bp + 1, b.start + bin_bp - 1
        da, db = abs(a.rel - 1), abs(b.rel - 1)
        if a.st != 0 and b.st != 0 and a.st != b.st:
            direction = "mixed"
        else:
            alt = a if da >= db else b
            direction = "gain" if alt.rel > 1 else "loss"
        dev = max(da, db)
        terminal = bool(first[i]) if da >= db else bool(last[i + 1])
        out.append(dict(chromosome=a.chromosome, pos=int(round((lo + hi) / 2)), lo=int(lo), hi=int(hi),
                        left_rel=round(a.rel, 4), right_rel=round(b.rel, 4), deviation=round(dev, 4),
                        stratum=stratum(dev), direction=direction, terminal=terminal))
    return out


def collapse(be, link_bp=LINK_BP):
    be = be.sort_values(["donor", "chromosome", "direction", "type", "pos"]).reset_index(drop=True)
    eid, k, prev = [], 0, None
    for r in be.itertuples():
        key = (r.donor, r.chromosome, r.direction, r.type)
        if prev is None or key != prev[0] or r.pos - prev[1] > link_bp:
            k += 1
        eid.append(k)
        prev = (key, r.pos)
    be["event_id"] = eid
    ev = be.groupby("event_id").agg(donor=("donor", "first"), chromosome=("chromosome", "first"),
                                    direction=("direction", "first"), type=("type", "first"),
                                    pos=("pos", "median"), lo=("lo", "median"), hi=("hi", "median"),
                                    deviation=("deviation", "median"), n_breakends=("cell_id", "size"),
                                    n_cells=("cell_id", "nunique"), rt=("rt", "median"),
                                    terminal=("terminal", "mean")).reset_index()
    ev["terminal"] = ev.terminal >= 0.5
    ev["stratum"] = ev.deviation.map(stratum)
    ev["main_def"] = ev.stratum.isin(["E2", "E3", "E4"])
    for c in ("pos", "lo", "hi"):
        ev[c] = ev[c].round().astype(int)
    return be, ev


def annotate(be, bins):
    be["type"] = ["arm_level" if overlaps_pericen(r.chromosome, r.lo, r.hi) else "interstitial" for r in be.itertuples()]
    be["rt"] = [interval_rt(bins, r.chromosome, r.lo, r.hi) for r in be.itertuples()]
    return be


def call_cells(seg, cells, bin_bp):
    rows = []
    for c, s in seg[seg.cell_id.isin(cells)].groupby("cell_id", sort=False):
        for r in cell_breakends(s, bin_bp):
            r.update(cell_id=c, donor=c.split("_")[0])
            rows.append(r)
    return pd.DataFrame(rows)


def refine(be500, be200):
    """Nearest 200 kb breakend of the same cell, chromosome and direction within [lo - 500 kb, hi + 500 kb]."""
    g200 = {k: g for k, g in be200.groupby(["cell_id", "chromosome", "direction"])}
    rows, matched = [], []
    for r in be500.itertuples():
        g = g200.get((r.cell_id, r.chromosome, r.direction))
        hit = None
        if g is not None:
            k = g[(g.pos >= r.lo - 500_000) & (g.pos <= r.hi + 500_000)]
            if len(k):
                hit = k.iloc[(k.pos - r.pos).abs().argmin()]
        matched.append(hit is not None)
        if hit is not None:
            rows.append(dict(chromosome=r.chromosome, pos=int(hit.pos), lo=int(hit.lo), hi=int(hit.hi),
                             left_rel=r.left_rel, right_rel=r.right_rel, deviation=r.deviation, stratum=r.stratum,
                             direction=r.direction, terminal=r.terminal, cell_id=r.cell_id, donor=r.donor))
    return pd.DataFrame(rows), np.array(matched)


def write_set(d, b, be, ev, bins):
    os.makedirs(d, exist_ok=True)
    be.to_csv(f"{d}/breakends_{b}kb.tsv", sep="\t", index=False)
    ev.to_csv(f"{d}/events_{b}kb.tsv", sep="\t", index=False)
    bins.to_csv(f"{d}/usable_bins_rt_{b}kb.tsv", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--profiles", required=True, help="dir with profiles_<bin>kb_bins.tsv.gz (export_relcn.R)")
    ap.add_argument("--rt", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--link-bp", type=int, default=LINK_BP)
    ap.add_argument("--collapse-on", default="all", choices=["all", "main"])
    ap.add_argument("--sets", default="main,nocuts,samecells,cells500,refined")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    sets = a.sets.split(",")

    qc = {b: pd.read_csv(f"{a.qc}/cn_qc_{b}kb.tsv", sep="\t") for b in BINS}
    passing = {b: set(q.cell_id[(q.depth_ok == 1) & (q.noise_excess <= 1.5) & (q.rt_rho <= 0.10)]) for b, q in qc.items()}
    same = set.intersection(*passing.values())
    print("cells passing in all three bins:", len(same))
    summ, bins_rt, seg_all, be_main500 = [], {}, {}, None
    for b in BINS:
        bin_bp = b * 1000
        bins = usable_bins(f"{a.profiles}/profiles_{b}kb_bins.tsv.gz")
        bins["rt"] = rt_per_bin(a.rt, bins, bin_bp)
        bins_rt[b] = bins
        qedges = np.nanquantile(bins.rt, [.2, .4, .6, .8])
        q = qc[b]
        seg = pd.read_csv(f"{a.qc}/breakends_{b}kb_segments.tsv", sep="\t", dtype={"chromosome": str})
        seg_all[b] = seg
        for name, cut in CUTS.items():
            if name not in sets and not (name == "main" and "refined" in sets):
                continue
            if name == "cells500":
                cells = qc[500][qc[500].cell_id.isin(passing[500])]
            else:
                cells = q[(q.depth_ok == 1) & (q.noise_excess <= cut["noise"]) & (q.rt_rho <= cut["rt"])]
            if name == "samecells":
                cells = cells[cells.cell_id.isin(same)]
            be = annotate(call_cells(seg, set(cells.cell_id), bin_bp), bins)
            if name == "main" and b == 500:
                be_main500 = be.copy()
            if name not in sets:
                continue
            if a.collapse_on == "main":
                be = be[be.stratum.isin(["E2", "E3", "E4"])]
            be, ev = collapse(be, a.link_bp)
            ev["rt_quintile"] = np.searchsorted(qedges, ev.rt) + 1
            ev.loc[ev.rt.isna(), "rt_quintile"] = np.nan
            write_set(f"{a.out}/{name}{a.tag}", b, be, ev, bins)
            sm = ev.groupby(["donor", "type", "direction", "stratum"]).size().rename("events").reset_index()
            sm.insert(0, "bin_kb", b); sm.insert(0, "set", name + a.tag)
            sm["cells_in_set"] = sm.donor.map(cells.groupby("donor").size())
            summ.append(sm)
            print(f"{name}{a.tag} {b} kb: {len(cells)} cells, {len(be)} breakends, {len(ev)} events "
                  f"({(ev.type == 'interstitial').sum()} interstitial, {ev.main_def.sum()} main_def)")

    if "refined" in sets:
        b = REFINE_BIN
        bins = bins_rt[b]
        be200 = call_cells(seg_all[b], set(be_main500.cell_id), b * 1000)
        rb, matched = refine(be_main500, be200)
        rb = annotate(rb, bins)
        be_main500["refined"] = matched
        rs = be_main500[be_main500.stratum.isin(["E2", "E3", "E4"])].groupby(["direction", "type"]).refined.agg(["size", "sum", "mean"])
        rs.columns = ["breakends_500kb", "refined", "fraction"]
        if a.collapse_on == "main":
            rb = rb[rb.stratum.isin(["E2", "E3", "E4"])]
        rb, ev = collapse(rb, a.link_bp)
        ev["rt_quintile"] = np.searchsorted(np.nanquantile(bins.rt, [.2, .4, .6, .8]), ev.rt) + 1
        ev.loc[ev.rt.isna(), "rt_quintile"] = np.nan
        d = f"{a.out}/refined{a.tag}"
        write_set(d, b, rb, ev, bins)
        rs.reset_index().to_csv(f"{d}/refinement_summary.tsv", sep="\t", index=False)
        print("refined (E2-E4 breakends of main 500 kb matched at 200 kb):")
        print(rs.round(3).to_string())
    pd.concat(summ).to_csv(f"{a.out}/summary{a.tag}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

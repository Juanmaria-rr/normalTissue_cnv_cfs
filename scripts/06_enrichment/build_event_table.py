#!/usr/bin/env python3
"""Breakend and event tables for the CFS enrichment tests, per bin.

Cells (two sets):
  main    depth_ok, noise_excess <= 1.5, rt_rho <= 0.10 (decisions 2026-11-04)
  nocuts     depth_ok only (sensitivity)
  samecells  main cuts passed in all three bins (sensitivity: the same cells in every bin)
Breakends (per cell), from the run-length segments of call_breakends.R:
  - state per segment with the permissive rule: gain if ratio > 1.2, loss if < 0.8;
    consecutive segments with the same state are merged (ratio = length-weighted mean);
    merged segments < 1 Mb are dropped
  - breakend = boundary between consecutive retained segments of a chromosome with
    different state; interval [lo, hi] = start of the last bin of the left segment to end
    of the first bin of the right one (one bin either side plus any gap)
  - deviation = max |ratio - 1| of the two flanks; stratum E1 [0.2, 0.3), E2 [0.3, 0.4),
    E3 [0.4, 0.6), E4 >= 0.6; main_def = E2-E4 (the paper's 1.3 / 0.7)
  - direction = loss / gain from the flank with the larger deviation; "mixed" when one
    flank is a gain and the other a loss
  - type = arm_level if [lo, hi] overlaps a pericentromeric block +/- 1 Mb, else interstitial
  - rt = mean replication timing of usable bins overlapping [lo, hi]; rt_quintile from
    the genome-wide distribution of usable-bin RT at that bin size
Events (per donor): breakends of the same chromosome, direction and type whose positions
  lie within 2 Mb of the previous one (single linkage on sorted positions) form one event.
  Event pos / lo / hi = medians of its breakends; deviation = median; n_cells = distinct
  cells. A clone therefore counts once; two independent breaks closer than 2 Mb also merge.

Outputs (<out>/<set>/): breakends_<bin>kb.tsv (with event_id), events_<bin>kb.tsv,
summary.tsv (events per bin x donor x type x direction x stratum).

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
        "samecells": dict(noise=1.5, rt=0.10)}


def stratum(d):
    for lo, hi, s in STRATA:
        if lo <= d < hi:
            return s
    return None


def cell_breakends(s, bin_bp):
    s = s.reset_index(drop=True)
    st = np.where(s.rel > 1 + T, 1, np.where(s.rel < 1 - T, -1, 0))
    ln = (s.end - s.start + 1).values
    brk = np.r_[True, (s.chromosome.values[1:] != s.chromosome.values[:-1]) | (st[1:] != st[:-1])]
    g = np.cumsum(brk)
    m = pd.DataFrame({"g": g, "chromosome": s.chromosome.values, "start": s.start.values, "end": s.end.values,
                      "st": st, "w": ln, "rw": s.rel.values * ln})
    m = m.groupby("g").agg(chromosome=("chromosome", "first"), start=("start", "min"), end=("end", "max"),
                           st=("st", "first"), w=("w", "sum"), rw=("rw", "sum")).reset_index(drop=True)
    m["rel"] = m.rw / m.w
    m = m[(m.end - m.start + 1) >= MIN_BP].reset_index(drop=True)
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
        out.append(dict(chromosome=a.chromosome, pos=int(round((lo + hi) / 2)), lo=int(lo), hi=int(hi),
                        left_rel=round(a.rel, 4), right_rel=round(b.rel, 4), deviation=round(dev, 4),
                        stratum=stratum(dev), direction=direction))
    return out


def collapse(be):
    be = be.sort_values(["donor", "chromosome", "direction", "type", "pos"]).reset_index(drop=True)
    eid, k, prev = [], 0, None
    for r in be.itertuples():
        key = (r.donor, r.chromosome, r.direction, r.type)
        if prev is None or key != prev[0] or r.pos - prev[1] > LINK_BP:
            k += 1
        eid.append(k)
        prev = (key, r.pos)
    be["event_id"] = eid
    ev = be.groupby("event_id").agg(donor=("donor", "first"), chromosome=("chromosome", "first"),
                                    direction=("direction", "first"), type=("type", "first"),
                                    pos=("pos", "median"), lo=("lo", "median"), hi=("hi", "median"),
                                    deviation=("deviation", "median"), n_breakends=("cell_id", "size"),
                                    n_cells=("cell_id", "nunique"), rt=("rt", "median")).reset_index()
    ev["stratum"] = ev.deviation.map(stratum)
    ev["main_def"] = ev.stratum.isin(["E2", "E3", "E4"])
    for c in ("pos", "lo", "hi"):
        ev[c] = ev[c].round().astype(int)
    return be, ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--profiles", required=True, help="dir with profiles_<bin>kb_bins.tsv.gz (export_relcn.R)")
    ap.add_argument("--rt", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # samecells: cells that pass the main cuts in all three bins, so bins compare the same cells
    passing = []
    for b in BINS:
        q = pd.read_csv(f"{a.qc}/cn_qc_{b}kb.tsv", sep="\t")
        passing.append(set(q.cell_id[(q.depth_ok == 1) & (q.noise_excess <= 1.5) & (q.rt_rho <= 0.10)]))
    same = set.intersection(*passing)
    print("cells passing in all three bins:", len(same))
    summ = []
    for b in BINS:
        bin_bp = b * 1000
        bins = usable_bins(f"{a.profiles}/profiles_{b}kb_bins.tsv.gz")
        bins["rt"] = rt_per_bin(a.rt, bins, bin_bp)
        qedges = np.nanquantile(bins.rt, [.2, .4, .6, .8])
        q = pd.read_csv(f"{a.qc}/cn_qc_{b}kb.tsv", sep="\t")
        seg = pd.read_csv(f"{a.qc}/breakends_{b}kb_segments.tsv", sep="\t", dtype={"chromosome": str})
        for name, cut in CUTS.items():
            cells = q[(q.depth_ok == 1) & (q.noise_excess <= cut["noise"]) & (q.rt_rho <= cut["rt"])]
            if name == "samecells":
                cells = cells[cells.cell_id.isin(same)]
            rows = []
            for c, s in seg[seg.cell_id.isin(cells.cell_id)].groupby("cell_id", sort=False):
                for r in cell_breakends(s, bin_bp):
                    r.update(cell_id=c, donor=c.split("_")[0])
                    rows.append(r)
            be = pd.DataFrame(rows)
            be["type"] = ["arm_level" if overlaps_pericen(r.chromosome, r.lo, r.hi) else "interstitial" for r in be.itertuples()]
            be["rt"] = [interval_rt(bins, r.chromosome, r.lo, r.hi) for r in be.itertuples()]
            be, ev = collapse(be)
            ev["rt_quintile"] = np.searchsorted(qedges, ev.rt) + 1
            ev.loc[ev.rt.isna(), "rt_quintile"] = np.nan
            d = f"{a.out}/{name}"
            os.makedirs(d, exist_ok=True)
            be.to_csv(f"{d}/breakends_{b}kb.tsv", sep="\t", index=False)
            ev.to_csv(f"{d}/events_{b}kb.tsv", sep="\t", index=False)
            bins.to_csv(f"{d}/usable_bins_rt_{b}kb.tsv", sep="\t", index=False)
            sm = ev.groupby(["donor", "type", "direction", "stratum"]).size().rename("events").reset_index()
            sm.insert(0, "bin_kb", b); sm.insert(0, "set", name)
            sm["cells_in_set"] = sm.donor.map(cells.groupby("donor").size())
            summ.append(sm)
            print(f"{name} {b} kb: {len(cells)} cells, {len(be)} breakends, {len(ev)} events "
                  f"({(ev.type == 'interstitial').sum()} interstitial, {ev.main_def.sum()} main_def)")
    pd.concat(summ).to_csv(f"{a.out}/summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

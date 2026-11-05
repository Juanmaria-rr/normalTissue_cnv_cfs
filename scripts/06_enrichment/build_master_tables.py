#!/usr/bin/env python3
"""Master tables of the pilot enrichment analysis: cells, breakends and events with every
annotation used by the tests, in one place.

  master_cells.tsv      one row per cell x bin: metadata (donor, assay, instrument, read
                        length, runs), alignment QC (dup %), CN QC (reads, depth_ok,
                        noise_excess, rt_rho, ploidy_best, fit_disagree, ploidy_flag),
                        pass_cuts (depth_ok & noise <= 1.5 & rt_rho <= 0.10), and counts of
                        breakends / interstitial losses / gains (main definition E2-E4)
  master_breakends.tsv  one row per breakend x cell x bin (main and nocuts sets, column
                        `set`): interval, flank ratios, deviation, stratum, direction, type,
                        RT, rt_quintile, event_id, and one 0/1 column per feature set
                        (`ov:<set>`) = interval overlaps the set
  master_events.tsv     one row per collapsed event x bin x set, with `ov:<set>` and the
                        exact null overlap probability `p:<set>` used by tests A
Feature sets: as in test_a_permutation.py. Also writes master_tables.xlsx (cells and events;
breakends are too many for a sheet in the full dataset, kept as TSV) when openpyxl exists.

Usage: python3 scripts/06_enrichment/build_master_tables.py --res results/enrichment/pilot \
    --qc results/pilot_cn/qc --align data/pilot/align_qc_pilot.tsv --cells data/metadata/cells.tsv \
    --cfs data/features/cfs/cfs_catalogs_hg19.tsv --genes data/features/genes
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from test_a_permutation import Null, feature_sets, overlaps  # noqa: E402

BINS = [200, 500, 1000]
SETS = ["main", "nocuts"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--align", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--genes", required=True)
    a = ap.parse_args()
    fs = feature_sets(a.cfs, a.genes)

    meta = pd.read_csv(a.cells, sep="\t", low_memory=False)
    meta = meta.groupby("cell_id").agg(instrument=("instrument", lambda x: "|".join(sorted(set(x)))),
                                       layout=("layout", lambda x: "|".join(sorted(set(x)))),
                                       read_length=("read_length", lambda x: "|".join(map(str, sorted(set(x))))),
                                       runs=("run_accession", "size")).reset_index()
    al = pd.read_csv(a.align, sep="\t")[["cell_id", "reads_in", "mapped_pct", "dup_pct"]]
    cells = []
    for b in BINS:
        q = pd.read_csv(f"{a.qc}/cn_qc_{b}kb.tsv", sep="\t")
        q = q.drop(columns=[c for c in ("instrument", "layout", "read_length") if c in q.columns])
        q["pass_cuts"] = ((q.depth_ok == 1) & (q.noise_excess <= 1.5) & (q.rt_rho <= 0.10)).astype(int)
        be = pd.read_csv(f"{a.res}/nocuts/breakends_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        be["main_def"] = be.stratum.isin(["E2", "E3", "E4"])
        cnt = be[be.main_def].groupby("cell_id").agg(
            breakends_main=("pos", "size"),
            interstitial_loss=("direction", lambda x: int(((x == "loss") & (be.loc[x.index, "type"] == "interstitial")).sum())),
            interstitial_gain=("direction", lambda x: int(((x == "gain") & (be.loc[x.index, "type"] == "interstitial")).sum())))
        q = q.merge(cnt, left_on="cell_id", right_index=True, how="left")
        for c in ("breakends_main", "interstitial_loss", "interstitial_gain"):
            q[c] = q[c].fillna(0).astype(int)
        cells.append(q)
    cells = pd.concat(cells).merge(meta, on="cell_id", how="left").merge(al, on="cell_id", how="left")
    front = ["cell_id", "donor", "bin_kb", "coassay", "instrument", "layout", "read_length", "runs", "reads_in",
             "dup_pct", "reads", "depth_ok", "noise_excess", "rt_rho", "pass_cuts"]
    cells = cells[front + [c for c in cells.columns if c not in front]]
    cells.to_csv(f"{a.res}/master_cells.tsv", sep="\t", index=False)
    print("master_cells:", cells.shape)

    allbe, allev = [], []
    for s in SETS:
        for b in BINS:
            be = pd.read_csv(f"{a.res}/{s}/breakends_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
            ev = pd.read_csv(f"{a.res}/{s}/events_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
            bins = pd.read_csv(f"{a.res}/{s}/usable_bins_rt_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
            null = Null(bins)
            edges = np.nanquantile(bins.rt, [.2, .4, .6, .8])
            be["rt_quintile"] = np.searchsorted(edges, be.rt) + 1
            be["main_def"] = be.stratum.isin(["E2", "E3", "E4"])
            for name, feat in fs.items():
                be[f"ov:{name}"] = [int(overlaps(feat, c, lo, hi)[0]) for c, lo, hi in zip(be.chromosome, be.lo, be.hi)]
                ov, pp = [], []
                for r in ev.itertuples():
                    ov.append(int(overlaps(feat, r.chromosome, r.lo, r.hi)[0]))
                    if r.type != "interstitial" or not np.isfinite(r.rt_quintile):
                        pp.append(np.nan); continue
                    x, y = null.candidates(r.chromosome, int(r.rt_quintile), int(r.hi - r.lo))
                    pp.append(overlaps(feat, r.chromosome, x, y).mean() if len(x) else np.nan)
                ev[f"ov:{name}"] = ov
                ev[f"p:{name}"] = np.round(pp, 5)
            be.insert(0, "bin_kb", b); be.insert(0, "set", s)
            ev.insert(0, "bin_kb", b); ev.insert(0, "set", s)
            allbe.append(be); allev.append(ev)
            print(s, b, "done")
    be = pd.concat(allbe); ev = pd.concat(allev)
    be.to_csv(f"{a.res}/master_breakends.tsv", sep="\t", index=False)
    ev.to_csv(f"{a.res}/master_events.tsv", sep="\t", index=False)
    print("master_breakends:", be.shape, "master_events:", ev.shape)
    try:
        with pd.ExcelWriter(f"{a.res}/master_tables.xlsx") as x:
            cells.to_excel(x, sheet_name="cells", index=False)
            ev.to_excel(x, sheet_name="events", index=False)
            be[be.set == "main"].to_excel(x, sheet_name="breakends_main", index=False)
        print("wrote master_tables.xlsx")
    except ImportError:
        print("openpyxl not available: no xlsx")


if __name__ == "__main__":
    main()

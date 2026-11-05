#!/usr/bin/env python3
"""Select a panel of example cells per donor for visual inspection of the calibration questions.

Per donor, from cells with depth_ok at 500 kb (seed 1), each cell tagged with why it was
picked (a cell can carry several tags):
  random      8 random cells
  batch       4 random cells per batch (assay x instrument x read length)
  noisy       6 highest noise_excess
  rt          6 highest rt_rho
  altered     8 random cells with >= 1 non-neutral retained segment (default rule)
  shallow     4 random cells that pass depth at 1 Mb but not at 500 kb
Writes <out>.txt (one cell_id per line, input of export_relcn.R --cells) and <out>.tsv
(cell_id, donor, tags).

Usage: python3 scripts/03_cn_qc/select_example_cells.py --qc results/pilot_cn/qc \
    --out results/pilot_cn/notebook/example_cells
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    q5 = pd.read_csv(os.path.join(a.qc, "cn_qc_500kb.tsv"), sep="\t")
    q1 = pd.read_csv(os.path.join(a.qc, "cn_qc_1000kb.tsv"), sep="\t")
    c5 = pd.read_csv(os.path.join(a.qc, "breakends_500kb_cells.tsv"), sep="\t")
    q5 = q5.merge(c5[["cell_id", "n_gain_seg", "n_loss_seg"]], on="cell_id", how="left")
    q5["altered"] = (q5.n_gain_seg + q5.n_loss_seg) > 0
    q5["batch"] = q5.coassay.astype(str) + "|" + q5.instrument + "|" + q5.read_length.astype(str)
    ok1 = set(q1.cell_id[q1.depth_ok == 1])

    tags = {}
    def add(ids, tag):
        for c in ids:
            tags.setdefault(c, []).append(tag)

    for d, g in q5.groupby("donor"):
        ok = g[g.depth_ok == 1]
        add(ok.sample(min(8, len(ok)), random_state=a.seed).cell_id, "random")
        for _, b in ok.groupby("batch"):
            add(b.sample(min(4, len(b)), random_state=a.seed).cell_id, "batch")
        add(ok.nlargest(6, "noise_excess").cell_id, "noisy")
        add(ok.nlargest(6, "rt_rho").cell_id, "rt")
        alt = ok[ok.altered]
        add(alt.sample(min(8, len(alt)), random_state=a.seed).cell_id, "altered")
        sh = g[(g.depth_ok == 0) & g.cell_id.isin(ok1)]
        add(sh.sample(min(4, len(sh)), random_state=a.seed).cell_id, "shallow")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    out = pd.DataFrame([(c, c.split("_")[0], ",".join(t)) for c, t in tags.items()],
                       columns=["cell_id", "donor", "tags"]).sort_values("cell_id")
    out.to_csv(a.out + ".tsv", sep="\t", index=False)
    out.cell_id.to_csv(a.out + ".txt", index=False, header=False)
    print(len(out), "cells;", out.groupby("donor").size().to_dict())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Calibration tables for level-2 QC and breakend parameters (pilot, per bin).

Inputs (per bin): <qc>/cn_qc_<bin>kb.tsv (cn_qc_metrics.R) and
<qc>/breakends_<bin>kb_segments.tsv (call_breakends.R; all run-length segments with their
ratio `rel`, before the minimum-segment filter).

Outputs to --out:
  noise_by_batch.tsv   depth_ok cells per bin x donor x assay x instrument x read length:
                       n, median / p95 / p99 noise_excess. Shows whether excess noise follows
                       the donor or the batch.
  rt_rho_by_donor.tsv  depth_ok cells per bin x donor: median / p95 / p99 rt_rho, and % above
                       candidate cuts (P18 p99 of that bin, 0.05, 0.10).
  threshold_sweep.tsv  depth_ok cells passing the noise cut (P18 p99 of the bin) and
                       rt_rho <= 0.10, per bin x ratio threshold x minimum-segment rule:
                         p18_altered_pct  % of P18 cells with >= 1 retained non-neutral
                                          segment (specificity; paper: 1.0% aneuploid)
                         p08_6q_recall    % of P08 reference clone cells in which a gain
                                          segment overlapping chr6:70-93 Mb is retained
                         altered_pct_<donor>
  cut_options_exclusion.tsv  % of depth_ok cells excluded per bin x option x donor x assay,
                       for candidate noise and rt_rho cuts: global (P18 p99), absolute cap,
                       batch-relative (value / batch median > P18 p99 / P18 median; for rt_rho,
                       value - batch median > P18 p99 - P18 median). Batch = donor x assay x
                       instrument x read length.
  Reference clone: P08 cells that, at 1 Mb with the default rule (1.3/0.7, >= 1 Mb), carry
  a retained gain segment overlapping chr6:70-93 Mb and pass depth at 1 Mb. 1 Mb is the
  bin with most P08 cells passing depth. The reference is defined with the default rule,
  so recall at that rule and bin is 100% by construction.

Ratio thresholds are symmetric around 1 (gain > g, loss < 2 - g). Minimum-segment rules:
  mb1        length >= 1 Mb (current)
  mb1_bins3  length >= 1 Mb and >= 3 bins
  bins5      >= 5 bins

Usage: python3 scripts/03_cn_qc/calibration_tables.py --qc results/pilot_cn/qc \
    --out results/pilot_cn/calibration
"""
import argparse
import os

import numpy as np
import pandas as pd

BINS = [200, 500, 1000]
GAINS = [1.2, 1.25, 1.3, 1.35, 1.4]
RULES = {
    "mb1": lambda s: s.len_mb >= 1,
    "mb1_bins3": lambda s: (s.len_mb >= 1) & (s.n_bins >= 3),
    "bins5": lambda s: s.n_bins >= 5,
}
RT_MAX = 0.10


def load(qc, b):
    q = pd.read_csv(os.path.join(qc, f"cn_qc_{b}kb.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(qc, f"breakends_{b}kb_segments.tsv"), sep="\t",
                    dtype={"chromosome": str})
    s["len_mb"] = (s.end - s.start + 1) / 1e6
    return q, s


def altered_cells(s, gain, rule):
    k = s[RULES[rule](s)]
    alt = k[(k.rel > gain) | (k.rel < 2 - gain)]
    return set(alt.cell_id)


def gain_6q_cells(s, gain, rule):
    k = s[RULES[rule](s) & (s.chromosome == "6") & (s.rel > gain)
          & (s.start < 93e6) & (s.end > 70e6)]
    return set(k.cell_id)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    data = {b: load(a.qc, b) for b in BINS}
    q1, s1 = data[1000]
    ok1 = set(q1.cell_id[(q1.depth_ok == 1) & (q1.donor == "P08")])
    ref = gain_6q_cells(s1, 1.3, "mb1") & ok1
    print(f"P08 6q reference clone: {len(ref)} cells")

    nb, rt, sw = [], [], []
    for b, (q, s) in data.items():
        ok = q[q.depth_ok == 1].copy()
        p18 = ok[ok.donor == "P18"]
        k_noise = p18.noise_excess.quantile(0.99)
        rt99 = p18.rt_rho.quantile(0.99)

        for key, g in ok.groupby(["donor", "coassay", "instrument", "read_length"]):
            nb.append(dict(bin_kb=b, donor=key[0], coassay=key[1], instrument=key[2],
                           read_length=key[3], n=len(g),
                           median=g.noise_excess.median(), p95=g.noise_excess.quantile(.95),
                           p99=g.noise_excess.quantile(.99)))
        for d, g in ok.groupby("donor"):
            rt.append(dict(bin_kb=b, donor=d, n=len(g), median=g.rt_rho.median(),
                           p95=g.rt_rho.quantile(.95), p99=g.rt_rho.quantile(.99),
                           pct_gt_p18p99=100 * (g.rt_rho > rt99).mean(),
                           pct_gt_0_05=100 * (g.rt_rho > 0.05).mean(),
                           pct_gt_0_10=100 * (g.rt_rho > 0.10).mean()))

        keep = ok[(ok.noise_excess <= k_noise) & (ok.rt_rho <= RT_MAX)]
        ss = s[s.cell_id.isin(set(keep.cell_id))]
        ref_b = ref & set(keep.cell_id)
        for gthr in GAINS:
            for rule in RULES:
                alt = altered_cells(ss, gthr, rule)
                r = dict(bin_kb=b, gain=gthr, loss=round(2 - gthr, 2), min_seg=rule,
                         noise_cut=round(k_noise, 3), n_cells=len(keep),
                         p18_altered_pct=100 * np.mean([c in alt for c in keep.cell_id[keep.donor == "P18"]]),
                         p08_6q_ref_cells=len(ref_b),
                         p08_6q_recall=100 * len(gain_6q_cells(ss, gthr, rule) & ref_b) / max(1, len(ref_b)))
                for d, g in keep.groupby("donor"):
                    r[f"altered_pct_{d}"] = 100 * np.mean([c in alt for c in g.cell_id])
                sw.append(r)

    # exclusion under candidate cut options; batch = donor x assay x instrument x read length
    ex = []
    for b, (q, s) in data.items():
        ok = q[q.depth_ok == 1].copy()
        ok["batch"] = ok.donor + "|" + ok.coassay.astype(str) + "|" + ok.instrument + "|" + ok.read_length.astype(str)
        p18 = ok[ok.donor == "P18"]
        n99, n50 = p18.noise_excess.quantile(.99), p18.noise_excess.median()
        r99, r50 = p18.rt_rho.quantile(.99), p18.rt_rho.median()
        bmed_n = ok.groupby("batch").noise_excess.transform("median")
        bmed_r = ok.groupby("batch").rt_rho.transform("median")
        opts = {
            "noise_global_p18p99": ok.noise_excess > n99,
            "noise_cap_1.5": ok.noise_excess > 1.5,
            "noise_batch_relative": ok.noise_excess / bmed_n > n99 / n50,
            "rt_p18p99": ok.rt_rho > r99,
            "rt_fixed_0.10": ok.rt_rho > 0.10,
            "rt_batch_centred": ok.rt_rho - bmed_r > r99 - r50,
        }
        for name, flag in opts.items():
            for (d, c), g in flag.groupby([ok.donor, ok.coassay]):
                ex.append(dict(bin_kb=b, option=name, donor=d, coassay=c, n=len(g),
                               excluded_pct=100 * g.mean()))
    pd.DataFrame(ex).round(2).to_csv(os.path.join(a.out, "cut_options_exclusion.tsv"), sep="\t", index=False)

    for name, rows in (("noise_by_batch", nb), ("rt_rho_by_donor", rt), ("threshold_sweep", sw)):
        pd.DataFrame(rows).round(3).to_csv(os.path.join(a.out, f"{name}.tsv"), sep="\t", index=False)
    print("written:", a.out)


if __name__ == "__main__":
    main()

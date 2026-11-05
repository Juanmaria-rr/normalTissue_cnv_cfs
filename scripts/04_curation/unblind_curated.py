#!/usr/bin/env python3
"""Unblind a curated batch and write the inputs for downstream use.

Reads curation/data/<batch>/curated_fits.csv (written by the app; sample_id = blinded
code) and curation/keys/<batch>_key.tsv, and writes:
  <out>/<batch>_curated.tsv   cell_id, bin_kb, purity, ploidy, star_rating, notes, saved_at
  <out>/bestfitting.txt       sampleID, ploidy, purity (tab), the file bestsolution_absCN.R
                              reads from its own directory (/opt/fitting_absCN in the image)
Fails if a curated code is not in the key.

Usage:
  python3 scripts/04_curation/unblind_curated.py --batch test20_200kb \
      --data-dir curation/data --key-dir curation/keys --out results/curation
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--key-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fits = pd.read_csv(os.path.join(a.data_dir, a.batch, "curated_fits.csv"), dtype={"sample_id": str})
    key = pd.read_csv(os.path.join(a.key_dir, f"{a.batch}_key.tsv"), sep="\t", dtype={"sample": str})
    m = fits.merge(key, left_on="sample_id", right_on="sample", how="left")
    if m.cell_id.isna().any():
        raise SystemExit(f"codes not in key: {', '.join(m.loc[m.cell_id.isna(), 'sample_id'])}")

    os.makedirs(a.out, exist_ok=True)
    cols = ["cell_id", "bin_kb", "purity", "ploidy", "star_rating", "notes"] + \
           (["saved_at"] if "saved_at" in m.columns else [])
    m[cols].sort_values("cell_id").to_csv(os.path.join(a.out, f"{a.batch}_curated.tsv"), sep="\t", index=False)
    m.rename(columns={"cell_id": "sampleID"})[["sampleID", "ploidy", "purity"]].to_csv(
        os.path.join(a.out, "bestfitting.txt"), sep="\t", index=False)
    print(f"{len(m)} curated cells -> {a.out}")


if __name__ == "__main__":
    main()

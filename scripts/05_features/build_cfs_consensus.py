#!/usr/bin/env python3
"""Derived CFS catalogues: consensus A (evidence tiers), A + C and the union of everything.

Rules (fixed in WORKLOG 2026-11-05 (2) before testing):
  method classes  MiDAS (macheret2020, ji2020); RT-delay (brison2019, sarni2020); FISH
                  (letallec2013); cytogenetic (humcfs); CNV (wilson2015); cancer-SV (simpson2021)
  regions         merged overlapping intervals (no padding) of MiDAS, RT-delay, CNV and
                  cancer-SV sites; FISH and cytogenetic intervals add evidence by overlap only
  core            >= 1 replication class (MiDAS, RT-delay) and >= 1 breakage class (FISH,
                  cytogenetic, CNV, cancer-SV)
  reproducible    not core, >= 2 distinct cell types among the region's sources
  single          the rest (reported, not tested)
  A               core + reproducible
  A+C             A union breast_predicted (catalogue C)
  union           every site of every catalogue (autosomes and X), including C
Inputs: data/features/cfs/cfs_catalogs_hg19.tsv, data/features/cfs/breast_predicted_cfs_hg19.tsv.
Outputs (--out): cfs_consensus_regions.tsv (all regions with evidence columns and tier) and
cfs_derived_hg19.tsv (catalogue rows: consensusA with cell_type = tier, breast_predicted,
A_plus_C, union_all), read by the enrichment tests together with cfs_catalogs_hg19.tsv.

Usage: python3 scripts/05_features/build_cfs_consensus.py --cfs data/features/cfs/cfs_catalogs_hg19.tsv \
    --breast data/features/cfs/breast_predicted_cfs_hg19.tsv --out data/features/cfs
"""
import argparse

import pandas as pd

CLASS = {"macheret2020": "MiDAS", "ji2020": "MiDAS", "brison2019": "RT-delay", "sarni2020": "RT-delay",
         "letallec2013": "FISH", "humcfs": "cytogenetic", "wilson2015": "CNV", "simpson2021": "cancer-SV"}
BOUNDARY = {"MiDAS", "RT-delay", "CNV", "cancer-SV"}
REPLICATION = {"MiDAS", "RT-delay"}
BREAKAGE = {"FISH", "cytogenetic", "CNV", "cancer-SV"}
COLS = ["catalog", "cell_type", "site_id", "name", "chrom", "start", "end", "source_build", "known_cfs", "note"]


def merge_regions(d):
    out = []
    for c, g in d.sort_values(["chrom", "start"]).groupby("chrom"):
        cur = None
        for r in g.itertuples():
            if cur is None or r.start > cur["end"]:
                if cur:
                    out.append(cur)
                cur = dict(chrom=c, start=r.start, end=r.end, members=[r.Index])
            else:
                cur["end"] = max(cur["end"], r.end); cur["members"].append(r.Index)
        out.append(cur)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfs", required=True)
    ap.add_argument("--breast", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    c = pd.read_csv(a.cfs, sep="\t")
    c["cell_type"] = c.cell_type.fillna("n/a")  # one Le Tallec row has no cell line
    c["cls"] = c.catalog.map(CLASS)
    b = pd.read_csv(a.breast, sep="\t")

    bnd = c[c.cls.isin(BOUNDARY)]
    ann = c[~c.cls.isin(BOUNDARY)]
    rows = []
    for k, r in enumerate(merge_regions(bnd)):
        m = bnd.loc[r["members"]]
        o = ann[(ann.chrom == r["chrom"]) & (ann.start < r["end"]) & (ann.end > r["start"])]
        cls = set(m.cls) | set(o.cls)
        cells = set(m.catalog + ":" + m.cell_type)
        tier = "core" if (cls & REPLICATION and cls & BREAKAGE) else "reproducible" if len(cells) >= 2 else "single"
        genes = sorted({n for n in m["name"].fillna("").astype(str) if n and n != "nan" and not n.startswith(("INTERGENIC", "SDR_"))})
        rows.append(dict(region_id=f"R{k + 1}", chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                         size_kb=round((r["end"] - r["start"]) / 1e3, 1), tier=tier,
                         classes=",".join(sorted(cls)), n_classes=len(cls),
                         sources=",".join(sorted(cells | set(o.catalog + ":" + o.cell_type))), n_cell_types=len(cells),
                         known_cfs=",".join(sorted({x for x in pd.concat([m.known_cfs, o.known_cfs]).dropna().astype(str) if x})),
                         names=";".join(genes[:6])))
    reg = pd.DataFrame(rows)
    reg.to_csv(f"{a.out}/cfs_consensus_regions.tsv", sep="\t", index=False)

    A = reg[reg.tier.isin(["core", "reproducible"])]
    derived = [pd.DataFrame(dict(catalog="consensusA", cell_type=A.tier, site_id=A.region_id, name=A.names,
                                 chrom=A.chrom, start=A.start, end=A.end, source_build="hg19",
                                 known_cfs=A.known_cfs, note=A.classes)),
               b.assign(catalog="breast_predicted")[COLS],
               pd.DataFrame(dict(catalog="A_plus_C", cell_type="consensus A", site_id=A.region_id, name=A.names,
                                 chrom=A.chrom, start=A.start, end=A.end, source_build="hg19", known_cfs=A.known_cfs,
                                 note=A.classes)),
               b.assign(catalog="A_plus_C", cell_type="breast predicted")[COLS],
               pd.concat([c[COLS], b[COLS]]).assign(catalog="union_all", cell_type="all sources")]
    d = pd.concat(derived, ignore_index=True)[COLS]
    d.to_csv(f"{a.out}/cfs_derived_hg19.tsv", sep="\t", index=False)
    print(reg.tier.value_counts().to_string())
    print(reg.groupby("tier").size_kb.sum().div(1e3).round(1).rename("Mb").to_string())
    d["mb"] = (d.end - d.start) / 1e6
    print(d.groupby(["catalog", "cell_type"]).mb.agg(["size", "sum"]).round(1).to_string())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Breast-predicted CFS (catalogue C): large genes accessible in breast epithelium and late-replicating.

Rules (fixed in WORKLOG 2026-11-05 (2) before testing):
  - genes > 300 kb: refGene hg19 spans (data/features/genes/large_genes_300kb.tsv)
  - accessible: promoter accessibility >= median over all refGene genes; promoter
    accessibility = max over the gene's TSSs of Tn5 insertions in TSS +/- 1 kb, pooled over
    the scATAC fragments of GSE272504 (19 donors, hg19), only the co-assay epithelial
    barcodes of the authors' map (HBCA_aneuploidy_coassay_barcode_cellname_dict.tsv, matched
    within each donor's file), in
    insertions per million
  - late: mean replication timing of the 50 kb windows inside the gene < genome-wide median
Insertions are counted at both fragment ends (Tn5 cut sites) into 100 bp bins per
chromosome, cached in <atac>/insertions_100bp.npz.

Outputs: <out>/breast_predicted_cfs_hg19.tsv (catalogue rows, same columns as
cfs_catalogs_hg19.tsv) and <out>/breast_gene_accessibility.tsv (all genes: TSS accessibility,
large / accessible / late flags).
Limitations: pooled over luminal and basal lineages (no lineage labels in the map);
promoter accessibility is a proxy for transcription; the RT annotation is not from breast.

Usage: python3 scripts/05_features/build_breast_predicted_cfs.py --atac data/features/atac/raw \
    --barcodes <dir>/HBCA_aneuploidy_coassay_barcode_cellname_dict.tsv --refgene data/features/genes/raw/refGene_hg19.txt.gz \
    --genes data/features/genes --rt external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv \
    --out data/features/cfs
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

BIN = 100
CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
COLS = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts",
        "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]


def insertions(atac, barcodes):
    cache = f"{atac}/insertions_100bp.npz"
    if os.path.exists(cache):
        z = np.load(cache)
        return {c: z[c] for c in z.files if c != "_total"}, int(z["_total"])
    bc = pd.read_csv(barcodes, sep="\t")
    keep_by_donor = {d: set(g.barcode) for d, g in bc.groupby("Patient_ID")}  # barcodes repeat across chips
    counts = {c: np.zeros(250_000_000 // BIN + 1, np.int32) for c in CHR}
    total, kept_frag = 0, 0
    for f in sorted(glob.glob(f"{atac}/*fragments.tsv.gz")):
        keep = keep_by_donor.get(os.path.basename(f).split("_")[1].split(".")[0], set())
        for ch in pd.read_csv(f, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3],
                              names=["chrom", "start", "end", "barcode"], chunksize=2_000_000):
            ch = ch[ch.barcode.isin(keep) & ch.chrom.isin(CHR)]
            kept_frag += len(ch)
            for c, g in ch.groupby("chrom"):
                for pos in (g.start.values, g.end.values - 1):
                    pos = pos[pos >= 0]  # a few fragments have negative coordinates (P21 file)
                    b = pos // BIN
                    counts[c] += np.bincount(b, minlength=len(counts[c]))[:len(counts[c])].astype(np.int32)
            total += 2 * len(ch)
        print(os.path.basename(f), "fragments kept so far:", kept_frag)
    np.savez_compressed(cache, _total=total, **counts)
    return counts, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atac", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--refgene", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--rt", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    counts, total = insertions(a.atac, a.barcodes)
    cum = {c: np.concatenate([[0], np.cumsum(v, dtype=np.int64)]) for c, v in counts.items()}

    def window(c, s, e):  # insertions in [s, e)
        if c not in cum:
            return 0
        i, j = max(s // BIN, 0), min(e // BIN + 1, len(cum[c]) - 1)
        return int(cum[c][j] - cum[c][i])

    r = pd.read_csv(a.refgene, sep="\t", header=None, names=COLS, usecols=["chrom", "strand", "txStart", "txEnd", "name2"])
    r = r[r.chrom.isin(CHR[:-1])]
    r["tss"] = np.where(r.strand == "+", r.txStart, r.txEnd)
    r["acc"] = [window(c, t - 1000, t + 1000) for c, t in zip(r.chrom, r.tss)]
    acc = r.groupby(["name2", "chrom"]).acc.max().reset_index()
    acc["acc_cpm"] = acc.acc / total * 1e6
    med = acc.acc_cpm.median()

    rt = pd.read_csv(a.rt, sep="\t", usecols=["seqnames", "start", "end", "replicationTime"])
    rt["seqnames"] = rt.seqnames.astype(str)
    rt_med = rt.replicationTime.median()
    big = pd.read_csv(f"{a.genes}/large_genes_300kb.tsv", sep="\t", dtype={"chromosome": str})
    big["chrom"] = "chr" + big.chromosome
    big = big.merge(acc, left_on=["gene", "chrom"], right_on=["name2", "chrom"], how="left")
    big["rt_mean"] = [rt[(rt.seqnames == c) & (rt.start >= s) & (rt.end <= e)].replicationTime.mean()
                      for c, s, e in zip(big.chromosome, big.start, big.end)]
    big["accessible"] = big.acc_cpm >= med
    big["late"] = big.rt_mean < rt_med
    acc.merge(big[["gene", "chrom", "span_kb", "rt_mean", "accessible", "late"]], left_on=["name2", "chrom"],
              right_on=["gene", "chrom"], how="left").drop(columns="gene").to_csv(
        f"{a.out}/breast_gene_accessibility.tsv", sep="\t", index=False)
    c = big[big.accessible & big.late]
    out = pd.DataFrame(dict(catalog="breast_predicted", cell_type="breast epithelium (GSE272504 scATAC, 19 donors)",
                            site_id=c.gene, name=c.gene, chrom=c.chrom, start=c.start, end=c.end, source_build="hg19",
                            known_cfs="", note=[f"span {s:.0f} kb; promoter {p:.1f} CPM; RT {t:.1f}"
                                                for s, p, t in zip(c.span_kb, c.acc_cpm, c.rt_mean)]))
    out.to_csv(f"{a.out}/breast_predicted_cfs_hg19.tsv", sep="\t", index=False)
    print(f"insertions {total:,}; promoter median {med:.2f} CPM; RT median {rt_med:.1f}")
    print(f"genes > 300 kb: {len(big)}; accessible {big.accessible.sum()}; late {big.late.sum()}; both {len(c)} "
          f"({(c.end - c.start).sum() / 1e6:.0f} Mb)")
    for g in ("FHIT", "WWOX", "MACROD2", "LRP1B", "IMMP2L", "PARK2", "PRKN", "CSMD1", "PTPRD", "NAALADL2"):
        x = big[big.gene == g]
        if len(x):
            x = x.iloc[0]
            print(f"  {g}: promoter {x.acc_cpm:.1f} CPM ({'≥' if x.accessible else '<'} median), RT {x.rt_mean:.1f} "
                  f"({'late' if x.late else 'early'}) -> {'in C' if x.accessible and x.late else 'not in C'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Large-gene features (hg19) for the enrichment controls and covariates.

Source: UCSC refGene (hg19). A gene's span is the union of its overlapping transcripts on one
autosome; transcripts with the same name at non-overlapping loci give separate spans. Writes to --out:
  genes_hg19.tsv          gene, chromosome, start (0-based), end, span_kb
  large_genes_300kb.tsv   span > 300 kb (CFS-associated size class, Le Tallec 2013)
  large_genes_1mb.tsv     span > 1 Mb (positive control: large transcription units, Wilson 2015)
Limitation: refGene spans are annotation-based; they do not say whether the gene is
transcribed in breast epithelium, which is what links large genes to fragility.

Usage: python3 scripts/05_features/build_gene_features.py --raw data/features/genes/raw --out data/features/genes
"""
import argparse
import os
import subprocess

import pandas as pd

URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
COLS = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts",
        "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.raw, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    f = os.path.join(a.raw, "refGene_hg19.txt.gz")
    if not os.path.exists(f):
        subprocess.run(["curl", "-sSLf", "-o", f, URL], check=True)
    r = pd.read_csv(f, sep="\t", header=None, names=COLS, usecols=["chrom", "txStart", "txEnd", "name2"])
    r = r[r.chrom.isin([f"chr{i}" for i in range(1, 23)])]
    # transcripts sharing a name can sit at distant loci (paralogues, mis-annotations): merge only
    # overlapping transcripts, so each locus is its own gene span
    r = r.sort_values(["name2", "chrom", "txStart"])
    rows = []
    for (gene, chrom), t in r.groupby(["name2", "chrom"], sort=False):
        s0, e0 = None, None
        for s, e in zip(t.txStart, t.txEnd):
            if s0 is None or s > e0:
                if s0 is not None:
                    rows.append((gene, chrom, s0, e0))
                s0, e0 = s, e
            else:
                e0 = max(e0, e)
        rows.append((gene, chrom, s0, e0))
    g = pd.DataFrame(rows, columns=["gene", "chromosome", "start", "end"])
    g["chromosome"] = g.chromosome.str.replace("chr", "")
    g["span_kb"] = ((g.end - g.start) / 1e3).round(1)
    g = g.sort_values(["chromosome", "start"])
    g.to_csv(f"{a.out}/genes_hg19.tsv", sep="\t", index=False)
    for lab, thr in (("300kb", 300), ("1mb", 1000)):
        x = g[g.span_kb > thr]
        x.to_csv(f"{a.out}/large_genes_{lab}.tsv", sep="\t", index=False)
        print(f"genes > {lab}: {len(x)} ({x.span_kb.sum() / 1e3:.0f} Mb)")
    print(len(g), "genes")


if __name__ == "__main__":
    main()

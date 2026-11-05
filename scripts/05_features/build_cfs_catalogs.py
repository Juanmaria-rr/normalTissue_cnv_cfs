#!/usr/bin/env python3
"""Download and harmonise common fragile site (CFS) catalogues into one hg19 table.

Catalogues (decision 2026-11-03: all analysed in parallel):
  macheret2020  MiDAS-seq domains, Macheret et al. 2020 Cell Res 30:997 (PMC7784693),
                Table S1. U2OS, HeLa, HS68 (normal fibroblasts). hg19, 10 kb bins
                (bin b = bp b*10000+1 .. (b+1)*10000).
  ji2020        MiDAS-seq sites (MDS), Ji et al. 2020 Cell Res 30:1009, Tables S1-S3.
                U2OS, HeLa, HBE (bronchial epithelium). hg38, 100 kb resolution; lifted to
                hg19. Tables are PDFs, parsed with `pdftotext -layout`.
  humcfs        HumCFS (Kumar et al. 2019 BMC Genomics), GFF3. Cytogenetic band
                coordinates in hg38 (e.g. FRA3B 3p14.2 = 58.6-63.7 Mb), lifted to hg19.
                Only common, aphidicolin-induced sites are kept. Band-level resolution:
                FRA16D (16q23.2) spans 79.2-81.7 Mb and misses most of WWOX
                (78.1-79.2 Mb).
  brison2019    aphidicolin significantly delayed regions (Repli-seq), JEFF lymphoblastoid cells,
                Brison et al. 2019 Nat Commun 10:5693, Supplementary Data 1, hg19.
  sarni2020     fragility-signature candidate genes (APH RT delay + transcription + TAD boundary),
                BJ-hTERT fibroblasts, Sarni et al. 2020 Nat Commun 11:3613, Supplementary Table 3,
                hg38 lifted to hg19. Candidates from a signature, not all validated as fragile.
  wilson2015    aphidicolin-induced CNV hotspot regions (>= 1 APH CNV), 090 and HF1 human
                fibroblasts, Wilson et al. 2015 Genome Res 25:189, Table S3A, hg19.
  simpson2021   the 21 fragile sites most altered by SVs in PCAWG cancers, Simpson et al. 2021
                Commun Biol 4:567 Table 1 (from Li et al. 2020 ED Fig. 9B), GRCh38 lifted to hg19.
  letallec2013  MCF10A (breast epithelium) sites of Le Tallec et al. 2013 Cell Rep 4:420.
                Table S4 (mmc4.xlsx, CFS mapped at the molecular level, 45 human CFS in
                several cell types, only 3 in MCF10A), hg19. Not downloadable
                automatically; read from --letallec when given.

Liftover (hg38 -> hg19, UCSC chain): start and end are lifted separately. An interval is
dropped if an end fails to lift, lands on another chromosome, or its length changes by
> 20%. Dropped intervals are listed in the log.

Outputs (--out):
  cfs_catalogs_hg19.tsv  catalog, cell_type, site_id, name, chrom, start (0-based), end,
                         source_build, known_cfs, note
  cfs_catalogs_summary.tsv  per catalog x cell type: n, total Mb, median/max size Mb
Raw downloads go to --raw (git-ignored).

Requires: pandas, xlrd, pyliftover, pdftotext (poppler).
Usage: python3 scripts/05_features/build_cfs_catalogs.py --raw data/features/cfs/raw \
    --out data/features/cfs [--letallec letallec_mcf10a.tsv]
"""
import argparse
import os
import re
import subprocess
import zipfile

import pandas as pd
from pyliftover import LiftOver

URLS = {
    "macheret_zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7784693/supplementaryFiles",
    "ji_S1": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0357-y/MediaObjects/41422_2020_357_MOESM8_ESM.pdf",
    "ji_S2": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0357-y/MediaObjects/41422_2020_357_MOESM9_ESM.pdf",
    "ji_S3": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0357-y/MediaObjects/41422_2020_357_MOESM10_ESM.pdf",
    "humcfs_gff": "https://webs.iiitd.edu.in/raghava/humcfs/humcfs_fragile.gff3.zip",
    "chain": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    "brison_zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6911102/supplementaryFiles",
    "sarni_zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7367836/supplementaryFiles",
    "wilson_zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4315293/supplementaryFiles",
}
AUTOSOMES_X = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
COLS = ["catalog", "cell_type", "site_id", "name", "chrom", "start", "end",
        "source_build", "known_cfs", "note"]


def fetch(url, path):
    # curl rather than urllib: python.org builds on macOS ship without CA certificates
    if not os.path.exists(path):
        subprocess.run(["curl", "-sSLf", "-A", "Mozilla/5.0", "-o", path, url], check=True)
    return path


def macheret(raw):
    z = fetch(URLS["macheret_zip"], os.path.join(raw, "macheret2020_supp.zip"))
    with zipfile.ZipFile(z) as zf:
        name = [n for n in zf.namelist() if n.endswith("MOESM9_ESM.xls")][0]
        xls = zf.extract(name, raw)
    v = pd.read_excel(xls, header=None)
    rows, cell = [], None
    for _, r in v.iterrows():
        lab = str(r[0])
        if lab.endswith(" cells"):
            cell = lab.replace(" cells", "")
            continue
        if cell is None or not str(r[1]).endswith("-pk"):
            continue
        chrom = f"chr{int(r[3])}" if str(r[3]).isdigit() else f"chr{r[3]}"
        rows.append(dict(catalog="macheret2020", cell_type=cell, site_id=f"{cell}_{int(r[0])}",
                         name=str(r[17]), chrom=chrom, start=int(r[4]) * 10000,
                         end=(int(r[6]) + 1) * 10000, source_build="hg19", known_cfs="",
                         note=str(r[1])))
    return pd.DataFrame(rows)


def ji(raw, lo):
    rows = []
    for key, cell in (("ji_S1", "U2OS"), ("ji_S2", "HeLa"), ("ji_S3", "HBE")):
        pdf = fetch(URLS[key], os.path.join(raw, f"ji2020_{key[-2:]}.pdf"))
        txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], check=True,
                             capture_output=True, text=True).stdout
        for line in txt.splitlines():
            # columns differ between tables (U2OS adds width, RT, origin density, gene length,
            # nascent RNA): take coordinates by position, the gene as the first token with a
            # letter after them, and known CFS as the FRA tokens
            m = re.match(r"\s*(\d+)\s+(chr\w+)\s+(\d+)\s+(\d+)\s+(.+?)\s*$", line)
            if not m:
                continue
            rest = m[5].split()
            gene = next((t for t in rest if re.search(r"[A-Za-z]", t) and not t.startswith(("FRA", "NA", "and"))), "")
            rows.append(dict(catalog="ji2020", cell_type=cell, site_id=f"{cell}_{m[1]}",
                             name=gene, chrom=m[2], start=int(m[3]), end=int(m[4]),
                             source_build="hg38", known_cfs=",".join(re.findall(r"FRA\w+", m[5])),
                             note="MiDAS-seq site (MDS)"))
    return lift(pd.DataFrame(rows), lo)


def humcfs(raw, lo):
    z = fetch(URLS["humcfs_gff"], os.path.join(raw, "humcfs_fragile.gff3.zip"))
    rows = []
    with zipfile.ZipFile(z) as zf:
        for n in zf.namelist():
            if not n.endswith(".gff3") or "__MACOSX" in n:
                continue
            for line in zf.read(n).decode().splitlines():
                f = line.split("\t")
                if len(f) > 9:  # FRAXB row has stray empty columns
                    f = [x for x in f if x != ""]
                if len(f) < 9 or f[2] != "FragileSite":
                    continue
                attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                attr = {k.lower(): v for k, v in attr.items()}
                rows.append(dict(catalog="humcfs", cell_type="cytogenetic (mostly lymphocytes)",
                                 site_id=attr.get("name", ""), name=attr.get("name", ""),
                                 chrom=f[0].replace("chrx", "chrX"), start=int(f[3]) - 1,
                                 end=int(f[4]), source_build="hg38", known_cfs=attr.get("name", ""),
                                 note=f"{attr.get('inducer', '')};{attr.get('frequency', '')};"
                                      f"band={attr.get('band', '')}"))
    d = pd.DataFrame(rows)
    keep = d.note.str.lower().str.contains("aphidicolin") & d.note.str.lower().str.contains("common")
    big = keep & (d.end - d.start > 20e6)  # annotation errors (e.g. FRA10D ends at 1,070 Mb)
    for _, r in d[big].iterrows():
        print(f"  humcfs dropped (> 20 Mb): {r.site_id} {r.chrom}:{r.start}-{r.end}")
    keep &= ~big
    print(f"humcfs: {len(d)} sites, {keep.sum()} common aphidicolin-induced <= 20 Mb kept")
    return lift(d[keep].copy(), lo)


def lift_point(lo, chrom, pos, step):
    # band ends at chromosome ends or in gaps do not lift; move inwards up to 1 Mb
    for k in range(101):
        x = lo.convert_coordinate(chrom, pos + k * step)
        if x:
            return x
    return []


def lift(d, lo):
    out, bad = [], []
    for _, r in d.iterrows():
        a = lift_point(lo, r.chrom, int(r.start), 10000)
        b = lift_point(lo, r.chrom, int(r.end) - 1, -10000)
        ok = bool(a) and bool(b) and a[0][0] == r.chrom and b[0][0] == r.chrom
        if ok:
            s, e = sorted((a[0][1], b[0][1] + 1))
            ok = abs((e - s) - (r.end - r.start)) <= 0.2 * (r.end - r.start)
        if not ok:
            bad.append(f"{r.catalog}:{r.site_id} {r.chrom}:{r.start}-{r.end}")
            continue
        r = r.copy(); r.start, r.end = s, e
        out.append(r)
    for x in bad:
        print("  liftover dropped:", x)
    return pd.DataFrame(out)


def letallec(path):
    # Table S4 (mmc4): CFS mapped at the molecular level, hg19, one row per CFS x cell type;
    # rows without a CFS name continue the previous one (extra large genes). Mouse (MEFs)
    # and chicken (DT40) rows and the 5 sites placed with a single BAC ("BAC", no
    # boundaries) are dropped.
    v = pd.read_excel(path, sheet_name="Table S4", header=0)
    v.columns = ["cfs", "cells", "chrom", "inf", "sup", "genes", "gene_kb", "pct_in_gene", "ref"]
    v["genes"] = v.genes.astype(str).where(v.genes.notna(), "")
    v["grp"] = v.cfs.notna().cumsum()
    genes = v.groupby("grp").genes.agg(lambda g: ",".join(x for x in g if x))
    v = v[v.cfs.notna()].copy()
    v["genes"] = v.grp.map(genes)
    v = v[~v.cells.isin(["MEFs", "DT40"]) & pd.to_numeric(v.inf, errors="coerce").notna()]
    v["cells"] = v.cells.fillna("n/a")
    name = v.cfs.str.replace(" ", "")
    return pd.DataFrame(dict(catalog="letallec2013", cell_type=v.cells.values,
                             site_id=(name + "_" + v.cells.str.replace(r"[./ ]+", "", regex=True)).values,
                             name=name.values, chrom=v.chrom.values,
                             start=v.inf.astype(int).values - 1, end=v.sup.astype(int).values,
                             source_build="hg19", known_cfs=name.values,
                             note=("large_genes=" + v.genes + ";ref=" + v.ref.astype(str)).values))



def supp_file(raw, key, name):
    z = fetch(URLS[key], os.path.join(raw, f"{key}.zip"))
    with zipfile.ZipFile(z) as zf:
        member = [n for n in zf.namelist() if n.endswith(name)][0]
        return zf.extract(member, raw)


def brison(raw):
    # Supplementary Data 1: aphidicolin significantly delayed regions (SDR), JEFF lymphoblastoid
    # cells, Repli-seq, hg19 ("chrN:start-end")
    v = pd.read_excel(supp_file(raw, "brison_zip", "MOESM4_ESM.xlsx"), header=None)
    rows = []
    for i, r in v.iterrows():
        m = re.match(r"(chr\w+):(\d+)-(\d+)", str(r[0]))
        if m:
            rows.append(dict(catalog="brison2019", cell_type="JEFF lymphoblastoid", site_id=f"SDR_{len(rows) + 1}",
                             name=str(r[6]) if pd.notna(r[6]) else "", chrom=m[1], start=int(m[2]), end=int(m[3]),
                             source_build="hg19", known_cfs=str(r[5]) if pd.notna(r[5]) else "",
                             note="aphidicolin significantly delayed region (Repli-seq)"))
    return pd.DataFrame(rows)


def sarni(raw, lo):
    # Supplementary Table 3: fragility-signature candidate genes, BJ-hTERT fibroblasts, hg38
    v = pd.read_excel(supp_file(raw, "sarni_zip", "MOESM6_ESM.xlsx"), header=1)
    v = v[v.Chr.astype(str).str.startswith("chr")]
    d = pd.DataFrame(dict(catalog="sarni2020", cell_type="BJ-hTERT fibroblast", site_id=v["Gene.Name"].astype(str),
                          name=v["Gene.Name"].astype(str), chrom=v.Chr, start=v.Start.astype(int) - 1, end=v.End.astype(int),
                          source_build="hg38", known_cfs="",
                          note="fragility-signature candidate gene (APH RT delay + transcription + TAD)"))
    return lift(d, lo)


def wilson(raw):
    # Table S3A: CNV regions, human fibroblasts 090 and HF1, hg19; kept if >= 1 aphidicolin-induced CNV
    f = supp_file(raw, "wilson_zip", "Table_S3A.xlsx")
    rows = []
    for sheet, cell in (("090", "090 fibroblast"), ("HF1", "HF1 fibroblast")):
        v = pd.read_excel(f, sheet_name=sheet, header=1)
        v = v[v.chrom.astype(str).str.startswith("chr") & (pd.to_numeric(v.APH, errors="coerce") >= 1)]
        for i, r in enumerate(v.itertuples()):
            rows.append(dict(catalog="wilson2015", cell_type=cell, site_id=f"{sheet}_{i + 1}", name=str(r.band),
                             chrom=r.chrom, start=int(r.start), end=int(r.end), source_build="hg19", known_cfs="",
                             note=f"aphidicolin CNV hotspot (APH CNVs = {int(r.APH)})"))
    return pd.DataFrame(rows)



# Simpson, Pye & Whitaker 2021 Commun Biol 4:567, Table 1: fragile sites most altered by SVs in
# PCAWG cancers (transcribed by the authors from Li et al. 2020 Nature, Extended Data Fig. 9B).
# GRCh38.p13 coordinates as printed; ";" = site between two adjacent genes.
SIMPSON2021 = [
    ("FHIT", "chr3", 59747277, 61251459), ("MACROD2", "chr20", 13995369, 16053197),
    ("WWOX", "chr16", 78099430, 79212667), ("IMMP2L", "chr7", 110662644, 111562517),
    ("NAALADL2", "chr3", 174438573, 175810548), ("LRP1B", "chr2", 140231423, 142131016),
    ("PDE4D", "chr5", 58969038, 60522120), ("CCSER1", "chr4", 90127535, 91601913),
    ("DMD", "chrX", 31097677, 33339441), ("PACRG;PARK2", "chr6", 161740845, 163315492),
    ("KIF26B;SMYD3", "chr1", 245353678, 246507312), ("PTPRD", "chr9", 8314246, 10612723),
    ("LSAMP", "chr3", 115802363, 117139389), ("AUTS2", "chr7", 69598296, 70793506),
    ("RBFOX1", "chr16", 5239802, 7713340), ("CSMD1", "chr8", 2935353, 4994972),
    ("PRKG1", "chr10", 50990888, 52298423), ("DIAPH2", "chrX", 96684663, 97604997),
    ("NEGR1", "chr1", 71395943, 72282539), ("GPC6", "chr13", 93226807, 94408020),
    ("CTNNA3", "chr10", 65912523, 67696195),
]


def simpson(lo):
    d = pd.DataFrame([dict(catalog="simpson2021", cell_type="PCAWG pan-cancer SV hotspots", site_id=g, name=g,
                           chrom=c, start=s - 1, end=e, source_build="hg38", known_cfs="",
                           note="most SV-altered fragile site (Li et al. 2020 ED Fig. 9B via Simpson 2021 Table 1)")
                      for g, c, s, e in SIMPSON2021])
    return lift(d, lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--letallec")
    a = ap.parse_args()
    os.makedirs(a.raw, exist_ok=True); os.makedirs(a.out, exist_ok=True)

    chain = fetch(URLS["chain"], os.path.join(a.raw, "hg38ToHg19.over.chain.gz"))
    lo = LiftOver(chain)
    parts = [macheret(a.raw), ji(a.raw, lo), humcfs(a.raw, lo), brison(a.raw), sarni(a.raw, lo), wilson(a.raw), simpson(lo)]
    if a.letallec:
        parts.append(letallec(a.letallec))
    d = pd.concat(parts, ignore_index=True)[COLS]
    d = d[d.chrom.isin(AUTOSOMES_X)].sort_values(["catalog", "cell_type", "chrom", "start"])
    d.to_csv(os.path.join(a.out, "cfs_catalogs_hg19.tsv"), sep="\t", index=False)

    d["mb"] = (d.end - d.start) / 1e6
    s = d.groupby(["catalog", "cell_type"]).agg(n=("mb", "size"), total_mb=("mb", "sum"),
                                                 median_mb=("mb", "median"), max_mb=("mb", "max"))
    s = s.round(2).reset_index()
    s.to_csv(os.path.join(a.out, "cfs_catalogs_summary.tsv"), sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()

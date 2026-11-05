#!/usr/bin/env python3
"""Build docs/enrichment/enrichment_pilot.html: methods and results of the CFS enrichment
analysis on the pilot, generated from the result tables (no hand-typed numbers).

Inputs: results/enrichment/pilot/{main,nocuts,samecells}/ (build_event_table.py,
test_a_permutation.py, test_b_regression.py, test_c_per_donor.py), funnel tables (funnel.py),
results/pilot_cn/qc/cn_qc_<bin>kb.tsv, data/features/cfs/cfs_catalogs_summary.tsv,
data/features/cfs/cfs_derived_hg19.tsv, data/features/genes/, and the figures of
scripts/docs/enrichment_figures.py and funnel.py (same output folder).

Usage: python3 scripts/docs/build_enrichment_html.py --res results/enrichment/pilot --out docs/enrichment
"""
import argparse
import html
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, f"{ROOT}/scripts/06_enrichment")
BINS = [200, 500, 1000]
MAIN_SETS = ["macheret2020:all", "macheret2020:U2OS", "macheret2020:HeLa", "macheret2020:HS68", "ji2020:all",
             "ji2020:HBE", "humcfs:all", "letallec2013:all", "letallec2013:epithelial", "letallec2013:MCF10A",
             "brison2019:all", "sarni2020:all", "wilson2015:all", "simpson2021:all", "consensusA:core",
             "consensusA:reproducible", "consensusA:all", "breast_predicted:all", "A_plus_C:all", "union_all:all",
             "genes_300kb", "genes_1mb"]


def fmt(x, d=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    if isinstance(x, (int, np.integer)):
        return f"{x:,}"
    return f"{x:.{d}f}"


def fp(p):
    if p is None or not np.isfinite(p):
        return "–"
    return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}"


def table(df, cols, num=()):
    h = "".join(f"<th{' class=num' if c in num else ''}>{html.escape(c)}</th>" for c in cols)
    rows = []
    for r in df[cols].itertuples(index=False):
        rows.append("<tr>" + "".join(f"<td{' class=num' if c in num else ''}>{html.escape(str(v))}</td>"
                                       for c, v in zip(cols, r)) + "</tr>")
    return f'<div class="tablebox"><table><thead><tr>{h}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def method(text):
    return f'<div class="method">{text}</div>'


def fig(name, caption):
    return f'<figure><img src="{name}.png" alt="{html.escape(caption[:120])}"><figcaption class="method">{caption}</figcaption></figure>'


def cells_section(res):
    rows = []
    for b in BINS:
        q = pd.read_csv(f"{ROOT}/results/pilot_cn/qc/cn_qc_{b}kb.tsv", sep="\t")
        for d, g in q.groupby("donor"):
            ok = g[g.depth_ok == 1]
            m = ok[(ok.noise_excess <= 1.5) & (ok.rt_rho <= 0.10)]
            rows.append({"bin": f"{b} kb", "donor": d, "fitted": fmt(len(g)), "depth_ok": fmt(len(ok)),
                         "noise > 1.5": fmt(int((ok.noise_excess > 1.5).sum())),
                         "rt_rho > 0.10": fmt(int((ok.rt_rho > 0.10).sum())), "main set": fmt(len(m))})
    return pd.DataFrame(rows)


def events_section(res):
    rows = []
    for b in BINS:
        e = pd.read_csv(f"{res}/main/events_{b}kb.tsv", sep="\t")
        be = pd.read_csv(f"{res}/main/breakends_{b}kb.tsv", sep="\t")
        for (t, d), g in e.groupby(["type", "direction"]):
            r = {"bin": f"{b} kb", "type": t, "direction": d,
                 "breakends": fmt(int(be[(be.type == t) & (be.direction == d)].shape[0])), "events": fmt(len(g))}
            for s in ("E1", "E2", "E3", "E4"):
                r[s] = fmt(int((g.stratum == s).sum()))
            r["in 1 cell"] = fmt(int((g.n_cells == 1).sum()))
            rows.append(r)
    return pd.DataFrame(rows)


def test_a_rows(t, direction, stratum="main"):
    x = t[t.feature_set.isin(MAIN_SETS) & (t.direction == direction) & (t.stratum == stratum)].copy()
    x["order"] = x.feature_set.map({s: i for i, s in enumerate(MAIN_SETS)})
    x = x.sort_values(["order", "bin_kb"])
    return pd.DataFrame({"set": x.feature_set, "bin": x.bin_kb.astype(str) + " kb", "events": x.n_events.map(lambda v: fmt(int(v))),
                         "obs.": x.observed.map(lambda v: fmt(int(v))), "exp.": x.expected.map(fmt),
                         "fold [95% CI]": [f"{fmt(f)} [{fmt(l)}–{fmt(h)}]" for f, l, h in zip(x.fold, x.fold_lo, x.fold_hi)],
                         "p (enrich.)": x.p_enrich.map(fp), "q": x.q_enrich.map(fp)})


def shuffle_rows(res):
    c = pd.read_csv(f"{res}/main/test_a_shuffle_controls.tsv", sep="\t")
    t = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    rows = []
    for name in MAIN_SETS:
        for b in BINS:
            for d in ("loss", "all"):
                s = c[(c.feature_set == name) & (c.bin_kb == b) & (c.direction == d)].fold.dropna()
                o = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == d) & (t.stratum == "main")]
                if s.empty or o.empty:
                    continue
                rows.append({"set": name, "bin": f"{b} kb", "direction": d, "shuffles": fmt(len(s)),
                             "median": fmt(s.median()), "2.5%": fmt(s.quantile(.025)), "97.5%": fmt(s.quantile(.975)),
                             "observed fold": fmt(o.fold.iloc[0])})
    return pd.DataFrame(rows)


def sensitivity_rows(res):
    runs = [("main", f"{res}/main/test_a.tsv"), ("no cuts", f"{res}/nocuts/test_a.tsv"),
            ("same cells", f"{res}/samecells/test_a.tsv")] + \
           [(f"without {d}", f"{res}/main/test_a_without_{d}.tsv") for d in ("P07", "P08", "P18", "P20", "P36")]
    rows = []
    for name in ("macheret2020:all", "macheret2020:HS68", "ji2020:all", "humcfs:all", "genes_1mb"):
        for b in BINS:
            r = {"set": name, "bin": f"{b} kb"}
            for lab, f in runs:
                t = pd.read_csv(f, sep="\t")
                x = t[(t.feature_set == name) & (t.bin_kb == b) & (t.direction == "loss") & (t.stratum == "main")]
                r[lab] = f"{fmt(x.fold.iloc[0])} ({int(x.observed.iloc[0])}/{fmt(x.expected.iloc[0], 1)})" if len(x) else "–"
            rows.append(r)
    return pd.DataFrame(rows)


def overlapping_events(res, name, b, direction="loss"):
    from test_a_permutation import feature_sets, overlaps
    fs = feature_sets(f"{ROOT}/data/features/cfs/cfs_catalogs_hg19.tsv", f"{ROOT}/data/features/genes")
    cfs = pd.read_csv(f"{ROOT}/data/features/cfs/cfs_catalogs_hg19.tsv", sep="\t")
    cfs["chromosome"] = cfs.chrom.str.replace("chr", "")
    cat, ct = name.split(":")
    sites = cfs[(cfs.catalog == cat) & ((cfs.cell_type == ct) | (ct == "all"))]
    e = pd.read_csv(f"{res}/main/events_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
    e = e[(e.type == "interstitial") & (e.direction == direction) & e.main_def]
    rows = []
    for r in e.itertuples():
        if overlaps(fs[name], r.chromosome, r.lo, r.hi)[0]:
            m = sites[(sites.chromosome == r.chromosome) & (sites.start <= r.hi) & (sites.end >= r.lo)]
            rows.append({"donor": r.donor, "position": f"chr{r.chromosome}:{fmt(r.pos / 1e6, 1)} Mb",
                         "interval": f"{fmt(r.lo / 1e6, 2)}–{fmt(r.hi / 1e6, 2)} Mb", "stratum": r.stratum,
                         "ratio (event median)": fmt(1 - r.deviation if direction == "loss" else 1 + r.deviation),
                         "cells": fmt(int(r.n_cells)), "site annotation (Macheret, Table S1)": "; ".join(sorted(set(m.name.astype(str))))})
    return pd.DataFrame(rows)


def test_b_rows(res):
    t = pd.read_csv(f"{res}/main/test_b.tsv", sep="\t")
    x = t[(t.direction == "loss") & (t.stratum == "main")].copy()
    x["order"] = x.feature_set.map({s: i for i, s in enumerate(MAIN_SETS)})
    x = x.sort_values(["order", "bin_kb", "model"], ascending=[True, True, False])
    return pd.DataFrame({"set": x.feature_set, "bin": x.bin_kb.astype(str) + " kb", "model": x.model,
                         "family": x.family, "events": x.n_events.map(lambda v: fmt(int(v))),
                         "IRR [95% CI]": [f"{fmt(i)} [{fmt(l)}–{fmt(h)}]" for i, l, h in zip(x.irr, x.irr_lo, x.irr_hi)],
                         "p": x.p.map(fp), "q": x.q.map(fp)}), t.family.value_counts().to_dict(), len(t)


def test_c_rows(res):
    t = pd.read_csv(f"{res}/main/test_c.tsv", sep="\t")
    x = t[(t.direction == "loss") & t.feature_set.isin(["macheret2020:all", "macheret2020:HS68", "ji2020:all", "humcfs:all", "genes_1mb"])]
    return pd.DataFrame({"set": x.feature_set, "bin": x.bin_kb.astype(str) + " kb", "donor": x.donor,
                         "cells": x.cells.map(lambda v: fmt(int(v))), "breakends": x.breakends.map(lambda v: fmt(int(v))),
                         "obs.": x.observed.map(lambda v: fmt(int(v))), "exp.": x.expected.map(fmt),
                         "O/E [95% CI]": [f"{fmt(r)} [{fmt(l)}–{fmt(h)}]" for r, l, h in zip(x.ratio, x.ratio_lo, x.ratio_hi)],
                         "heterogeneity p": x.het_p.map(fp)})


# How each site set is built: (set names in the tests, source, assay / what was measured, cell type,
# what one site is, processing, script)
CATALOG_METHODS = [
    ("macheret2020:U2OS / HeLa / HS68 / all", "Macheret et al. 2020, Cell Res 30:997, Table S1",
     "MiDAS-seq: EdU labelling of mitotic DNA synthesis after mild aphidicolin; marks regions still under-replicated at mitosis",
     "U2OS (osteosarcoma), HeLa (cervical cancer), HS68 (normal foreskin fibroblasts)",
     "A MiDAS domain (single, double or merged peak), given as 10 kb bins",
     "Native hg19; bin b = bp b·10,000+1…(b+1)·10,000; no filtering", "build_cfs_catalogs.py"),
    ("ji2020:U2OS / HeLa / HBE / all", "Ji et al. 2020, Cell Res 30:1009, Tables S1–S3",
     "MiDAS-seq (MiDASeq) after aphidicolin",
     "U2OS, HeLa, HBE (normal bronchial epithelium)",
     "A MiDAS site (MDS), 100 kb resolution",
     "PDF tables parsed (coordinates by position, nearest gene, FRA names); hg38 → hg19 liftover", "build_cfs_catalogs.py"),
    ("humcfs:all", "Kumar et al. 2019, BMC Genomics (HumCFS database)",
     "Cytogenetics: chromosome gaps/breaks in metaphase spreads, compiled from the literature",
     "Mostly peripheral blood lymphocytes",
     "A fragile-site band (cytogenetic band coordinates)",
     "Only common, aphidicolin-induced sites; 4 sites > 20 Mb dropped as annotation errors; hg38 → hg19 liftover", "build_cfs_catalogs.py"),
    ("letallec2013:<cell line> / all / epithelial", "Le Tallec et al. 2013, Cell Rep 4:420, Table S4",
     "FISH with BAC probes on metaphases after aphidicolin: CFS borders set where ≥ 85% of breaks fall",
     "12 lines incl. MCF10A and CAL51 (breast), HCT116/LoVo/LS174T (colon), lymphocytes, fibroblasts, K562",
     "A molecularly mapped CFS in one cell type (49 distinct CFS)",
     "Native hg19; mouse/chicken rows and 5 single-BAC sites dropped; ‘epithelial’ = HCT116, LoVo, LS174T, MCF10A, CAL51", "build_cfs_catalogs.py"),
    ("brison2019:all", "Brison et al. 2019, Nat Commun 10:5693, Supplementary Data 1",
     "Repli-seq: regions whose replication is significantly delayed by aphidicolin (SDR)",
     "JEFF lymphoblastoid cells",
     "A significantly delayed region", "Native hg19", "build_cfs_catalogs.py"),
    ("sarni2020:all", "Sarni et al. 2020, Nat Commun 11:3613, Supplementary Table 3",
     "Fragility signature: aphidicolin replication-timing delay (Repli-seq) + transcription (Bru-seq) + TAD boundary (Hi-C)",
     "BJ-hTERT fibroblasts (Hi-C from HFF-hTERT)",
     "A candidate gene carrying the signature (not all validated as fragile)",
     "hg38 → hg19 liftover (1 gene lost)", "build_cfs_catalogs.py"),
    ("wilson2015:all", "Wilson et al. 2015, Genome Res 25:189, Table S3A",
     "De novo CNVs in clones after aphidicolin; regions where CNVs cluster",
     "090 and HF1 normal human fibroblasts",
     "A CNV cluster region with ≥ 1 aphidicolin-induced CNV", "Native hg19", "build_cfs_catalogs.py"),
    ("simpson2021:all", "Simpson et al. 2021, Commun Biol 4:567, Table 1 (from Li et al. 2020 Nature, ED Fig. 9B)",
     "Somatic structural variants in cancer genomes (PCAWG): the fragile sites most altered by deletions/duplications",
     "Pan-cancer tumours",
     "One of the 21 most SV-altered fragile sites (gene span, or the region between two genes)",
     "GRCh38 coordinates as printed → hg19 liftover; chrX sites kept but outside the autosomal tests", "build_cfs_catalogs.py"),
    ("consensusA:core / reproducible / all", "Derived from the eight sources above",
     "Evidence integration (rules fixed before testing)",
     "All of the above",
     "A region = merged overlapping intervals of the sequencing / gene-resolution sources (MiDAS, RT-delay, CNV, cancer-SV); FISH and cytogenetic sites add evidence by overlap only",
     "Method classes: MiDAS, RT-delay, FISH, cytogenetic, CNV, cancer-SV. core = ≥ 1 replication class (MiDAS, RT-delay) and ≥ 1 breakage class (FISH, cytogenetic, CNV, cancer-SV); "
     "reproducible = not core and ≥ 2 cell types; single-source regions are listed but not tested; all = core + reproducible", "build_cfs_consensus.py"),
    ("breast_predicted:all", "Derived: genes from refGene; accessibility from GSE272504 (this study's co-assay)",
     "Prediction from CFS determinants: large genes that are open (promoter accessibility, proxy for transcription) and late-replicating",
     "Normal breast epithelium (28,547 co-assay cells of 19 donors, all lineages pooled)",
     "A gene > 300 kb",
     "Promoter accessibility = max over the gene's TSSs of Tn5 insertions in TSS ± 1 kb (per million insertions), ≥ median of all genes; "
     "mean replication timing of the gene < genome median", "build_breast_predicted_cfs.py"),
    ("A_plus_C:all", "Derived", "Union of consensusA:all and breast_predicted:all", "—", "A site of either set", "—", "build_cfs_consensus.py"),
    ("union_all:all", "Derived", "Union of every site of every catalogue above", "—", "Any site", "—", "build_cfs_consensus.py"),
    ("genes_300kb", "refGene hg19 (UCSC)",
     "Control set: gene length only (no fragility measurement). Le Tallec 2013 found CFS in genes > 300 kb",
     "— (annotation)", "A gene locus = union of overlapping transcripts with the same name, span > 300 kb",
     "Transcripts with the same name at non-overlapping loci give separate loci; not filtered by expression", "build_gene_features.py"),
    ("genes_1mb", "refGene hg19 (UCSC)",
     "Control set: gene length only. Wilson 2015 found transcription units > 1 Mb predictive of CFS and CNV hotspots",
     "— (annotation)", "A gene locus with span > 1 Mb", "As genes_300kb", "build_gene_features.py"),
]


THRESHOLDS = [
    ("Cell selection", "Well imaging: qc_flag = Good, cell_flag = yes, not a control well", "Excludes cell", "Authors' metadata"),
    ("Reads", "R1 only; cutadapt (Tn5 mosaic end), reads < 20 bp discarded", "Excludes reads", "Decided"),
    ("Alignment", "bwa mem to hg19; MAPQ ≥ 37 to be counted", "Excludes reads", "QDNAseq default"),
    ("Duplicates", "samtools markdup (runs of a cell concatenated first)", "Excludes reads", "—"),
    ("Mapping rate", "≥ 90%", "Alarm only", "—"),
    ("Minimum reads to fit", "≥ 20,000 usable reads", "Excludes cell", "Decided (pipeline robustness)"),
    ("Depth (depth_ok)", "≥ 15 × 2 × usable bins / 0.91: 396k / 145k / 64k reads (200 kb / 500 kb / 1 Mb)", "Excludes cell (per bin)", "fitting_absCN diploid threshold"),
    ("Noise", "noise_excess ≤ 1.5", "Excludes cell (per bin)", "Decided 2026-11-04"),
    ("S phase", "rt_rho ≤ 0.10", "Excludes cell (per bin)", "Decided 2026-11-04"),
    ("Absolute fit", "ploidy_flag (best ploidy ≥ 2.6), fit_disagree", "Alarm only", "Decided 2026-10-29"),
    ("Segments", "≥ 1 Mb", "Excludes segment", "Decided"),
    ("Segment state", "Ratio > 1.2 / < 0.8 (calling); main definition E2–E4 = 1.3 / 0.7", "Defines breakend", "Paper; strata decided 2026-11-04"),
    ("Type", "Arm-level if the interval overlaps the pericentromeric block ± 1 Mb", "Excluded from tests A–C", "Decided"),
    ("Clonal collapse", "Same donor, chromosome, direction and type, linked within ≤ 2 Mb", "Groups breakends", "Decided"),
]


def funnel_tables(res):
    fc = pd.read_csv(f"{res}/funnel_cells.tsv", sep="\t")
    fb = pd.read_csv(f"{res}/funnel_breakends.tsv", sep="\t")
    order = list(dict.fromkeys(fc.stage))
    tot = fc[fc.donor == "total"].pivot(index="stage", columns="bin_kb", values="cells").reindex(order)
    tot = tot.reset_index().rename(columns={200: "200 kb", 500: "500 kb", 1000: "1 Mb"})
    for c in ("200 kb", "500 kb", "1 Mb"):
        tot[c] = tot[c].map(lambda v: fmt(int(v)))
    don = fc[(fc.bin_kb == 500) & (fc.donor != "total")].pivot(index="stage", columns="donor", values="cells").reindex(order)
    don = don.reset_index()
    for c in don.columns[1:]:
        don[c] = don[c].map(lambda v: fmt(int(v)))
    border = list(dict.fromkeys(fb.stage))
    br = fb.pivot(index="stage", columns="bin_kb", values="n").reindex(border).reset_index().rename(
        columns={200: "200 kb", 500: "500 kb", 1000: "1 Mb"})
    for c in ("200 kb", "500 kb", "1 Mb"):
        br[c] = br[c].map(lambda v: fmt(int(v)))
    return tot, don, br


# plain-language meaning of each funnel stage (section 2.1)
STAGE_DEFS = [
    ("1 wells", "Every well of the 5 pilot donors listed in ENA (one well = one sequenced library, i.e. one putative cell), including control and empty wells."),
    ("2 imaging", "Wells that the authors' imaging scored as holding exactly one good-quality cell (qc_flag = Good, cell_flag = yes) and that are not control wells. These are the pilot cells."),
    ("3 aligned", "Cells whose reads were trimmed, aligned to hg19 and deduplicated (all imaging-selected cells)."),
    ("4 mapping ≥ 90%", "Cells with at least 90% of reads mapped to the genome. An alarm for a failed library; it does not remove cells."),
    ("5 fitted (≥ 20k reads)", "Cells with at least 20,000 usable reads (mapped, MAPQ ≥ 37, not duplicates); below that a profile is mostly empty and is not fitted."),
    ("6 depth_ok", "Cells with enough usable reads to tell one copy apart from noise: at least 15 reads per bin per chromosome copy in a diploid cell, i.e. 15 × 2 × usable bins / 0.91 = 396k, 145k and 64k reads at 200 kb, 500 kb and 1 Mb. Depends on the bin size."),
    ("7 noise_excess ≤ 1.5", "Cells whose bin-to-bin variation is at most 1.5 times what read-count (Poisson) sampling alone would give; noisier profiles produce false segments."),
    ("8 rt_rho ≤ 0.10", "Cells whose copy-number profile does not follow replication timing (Spearman ρ ≤ 0.10); higher values indicate a cell caught in S phase, whose profile imitates gains and losses. This is the main analysis set."),
    ("9 ≥ 1 breakend (E2–E4)", "Main-set cells with at least one copy-number change of the main definition (a segment ≥ 1 Mb with ratio > 1.3 or < 0.7 next to a different state)."),
    ("10 ≥ 1 interstitial", "Of those, cells with at least one breakend inside a chromosome arm (not at the centromere)."),
    ("11 ≥ 1 interstitial loss", "Of those, cells with at least one interstitial breakend whose altered side is a loss."),
]

CSS = """
:root{--bg:#f7f8f9;--surface:#fff;--fg:#1b232b;--muted:#56636f;--rule:#d6dde3;--accent:#0f6e73;--code:#eef1f4;
--font:"Public Sans",-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;--display:"Newsreader",Georgia,serif;--mono:"JetBrains Mono",ui-monospace,Menlo,monospace}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#11161b;--surface:#171e24;--fg:#e3e8ec;--muted:#9aa7b2;--rule:#2c3640;--accent:#5cc3c7;--code:#1d262e;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#11161b;--surface:#171e24;--fg:#e3e8ec;--muted:#9aa7b2;--rule:#2c3640;--accent:#5cc3c7;--code:#1d262e;color-scheme:dark}
body{background:var(--bg);color:var(--fg);font:15px/1.6 var(--font);padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px));padding-block:36px 90px}
main{max-width:1000px;margin:0 auto}h1,h2,h3{font-family:var(--display);font-weight:600;line-height:1.2;text-wrap:balance}
h1{font-size:2.1rem;margin:0}h2{font-size:1.55rem;margin:3rem 0 .6rem;padding-top:.4rem;border-top:1px solid var(--rule)}h3{font-size:1.15rem;margin:1.8rem 0 .4rem}
p,ul,ol{max-width:78ch}code{font:.85em var(--mono);background:var(--code);padding:.1em .35em;border-radius:4px}
.meta{color:var(--muted);font-size:.9rem}.eyebrow{font:600 .74rem/1 var(--mono);letter-spacing:.08em;text-transform:uppercase;color:var(--accent)}
.tablebox{overflow-x:auto;margin:.8rem 0 .3rem;border:1px solid var(--rule);border-radius:6px;background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:.82rem;font-variant-numeric:tabular-nums}th,td{padding:5px 9px;border-bottom:1px solid var(--rule);text-align:left;vertical-align:top;white-space:nowrap}
thead th{background:var(--code);font-weight:600}.num{text-align:right}
.method{font-size:.82rem;color:var(--muted);border-left:2px solid var(--rule);padding:.2rem 0 .2rem .8rem;margin:.4rem 0 1.4rem;max-width:90ch}
.method b{color:var(--fg)}figure{margin:1.2rem 0}figure img{width:100%;border-radius:6px;background:#fff}
pre{font:.8rem/1.5 var(--mono);background:var(--code);padding:.7rem .9rem;border-radius:6px;overflow-x:auto}
nav ol{columns:2;font-size:.88rem}nav a{color:var(--accent)}
.wrap td{white-space:normal;min-width:9rem;max-width:22rem}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = a.res
    ta = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    cat = pd.read_csv(f"{ROOT}/data/features/cfs/cfs_catalogs_summary.tsv", sep="\t")
    der = pd.read_csv(f"{ROOT}/data/features/cfs/cfs_derived_hg19.tsv", sep="\t")
    der["mb"] = (der.end - der.start) / 1e6
    der = der.groupby(["catalog", "cell_type"]).agg(n=("mb", "size"), total_mb=("mb", "sum"), median_mb=("mb", "median")).reset_index()
    g3 = pd.read_csv(f"{ROOT}/data/features/genes/large_genes_300kb.tsv", sep="\t")
    g1 = pd.read_csv(f"{ROOT}/data/features/genes/large_genes_1mb.tsv", sep="\t")
    tb, fams, nb = test_b_rows(res)
    hs = overlapping_events(res, "macheret2020:HS68", 500)

    parts = []
    P = parts.append
    P(f"""<title>CFS enrichment · pilot</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,600&family=Public+Sans:wght@400;600&family=JetBrains+Mono&display=swap">
<style>{CSS}</style><main>
<span class="eyebrow">CNV_cfs · pilot · hg19</span>
<h1>Breakend enrichment at common fragile sites: methods and pilot results</h1>
<p class="meta">Generated on {date.today().isoformat()} by <code>scripts/docs/build_enrichment_html.py</code> from the tables in <code>results/enrichment/pilot/</code>.
Every number in this document is read from those tables.</p>
<nav><ol><li><a href="#scope">Scope</a></li><li><a href="#funnel">Thresholds and funnel</a></li><li><a href="#cells">Cells</a></li><li><a href="#events">Breakends and events</a></li>
<li><a href="#features">Site sets</a></li><li><a href="#testa">Test A</a></li><li><a href="#testb">Test B</a></li>
<li><a href="#testc">Test C</a></li><li><a href="#regions">Regional views</a></li><li><a href="#limitations">Method limitations</a></li>
<li><a href="#repro">Reproducibility</a></li></ol></nav>""")

    P("""<h2 id="scope">1 · Scope</h2>
<p>Question: whether interstitial copy-number breakends of normal breast epithelial cells (Lin et al. 2024, PRJNA1003661) overlap common fragile sites (CFS) more than expected under a positional null stratified by replication timing.</p>
<ul><li>Data: 5-donor pilot (P07, P08, P18, P20, P36); 200 kb, 500 kb and 1 Mb bins analysed in parallel.</li>
<li>All tests are exploratory: there is no pre-registered primary hypothesis (decision of 4 November 2026). BH correction is applied within each site set (family), over bins × directions × strata.</li>
<li>The cell-selection and breakend-calling criteria were fixed before any overlap with CFS was computed. The rules of the derived catalogues (consensus, breast-predicted, A+C, union) were written down before they were tested.</li></ul>""")

    tot, don, br = funnel_tables(res)
    thr = pd.DataFrame(THRESHOLDS, columns=["step", "criterion and threshold", "effect", "origin"])
    P('<h2 id="funnel">2 · Thresholds and data funnel</h2><p>Criteria applied at each step, from the well to the event that enters the tests. '
      'Several act on reads or segments and do not exclude cells; two are alarms only.</p>')
    P(table(thr, list(thr.columns)))
    P(method("<b>Source:</b> decisions recorded in the project log and README; thresholds implemented in <code>scripts/03_cn_qc/cn_qc_metrics.R</code>, "
             "<code>scripts/06_enrichment/build_event_table.py</code> and the HealthOmics workflows."))
    P("<h3>2.1 Cells per stage (5 donors)</h3><p>Each stage keeps the cells of the previous one that meet one more criterion (stage 4 is reported but removes no cells).</p>")
    sd = pd.DataFrame(STAGE_DEFS, columns=["stage", "what it means"])
    P('<div class="tablebox wrap">' + table(sd, list(sd.columns))[len('<div class="tablebox">'):])
    P(table(tot, list(tot.columns), num=["200 kb", "500 kb", "1 Mb"]))
    P(fig("fig_funnel", "<b>What is shown.</b> Cells left after each criterion, cumulative, per bin (log scale); in brackets, % of the imaging-selected cells (stage 2). "
          "Stages 1–4 do not depend on the bin. The “mapping ≥ 90%” stage does not exclude cells and is not drawn. <b>Source:</b> <code>results/enrichment/pilot/funnel_cells.tsv</code> (script <code>scripts/06_enrichment/funnel.py</code>)."))
    P(method("<b>What is counted.</b> Pilot cells (wells). Stages 1–2: ENA and authors' metadata; 3: aligned cells; 5–8: <code>master_cells</code>; "
             "9–11: main-set cells with at least one breakend of the main definition (E2–E4), interstitial, and interstitial loss. "
             "<b>Limitation.</b> Stage 1 counts every well of the donor in ENA, including control and empty wells. <b>Source:</b> <code>data/metadata/cells.tsv</code>, <code>data/pilot/align_qc_pilot.tsv</code>, <code>results/enrichment/pilot/master_cells.tsv</code>, <code>results/enrichment/pilot/funnel_cells.tsv</code> (script <code>scripts/06_enrichment/funnel.py</code>)."))
    P("<h3>2.2 Cells per stage and donor (500 kb)</h3>")
    P(table(don, list(don.columns), num=list(don.columns[1:])))
    P(method("<b>Source:</b> <code>results/enrichment/pilot/funnel_cells.tsv</code> (script <code>scripts/06_enrichment/funnel.py</code>)."))
    P("<h3>2.3 From breakends to events (main set)</h3>")
    P(table(br, list(br.columns), num=["200 kb", "500 kb", "1 Mb"]))
    P(method("<b>What is counted.</b> Breakends of the main-set cells (one row per cell and breakend); then interstitial E2–E4 collapsed events, the loss events, "
             "and how many of these overlap MiDAS domains (Macheret union and HS68 only). <b>Source:</b> <code>results/enrichment/pilot/funnel_breakends.tsv</code>, <code>results/enrichment/pilot/master_breakends.tsv</code>, <code>results/enrichment/pilot/master_events.tsv</code> (script <code>scripts/06_enrichment/funnel.py</code>)."))

    P('<h2 id="cells">3 · Cells</h2><p>Criteria: <code>depth_ok</code> (usable reads ≥ 15 × 2 × usable bins / 0.91), '
      '<code>noise_excess</code> ≤ 1.5 and <code>rt_rho</code> ≤ 0.10, the same in the three bins. The “no cuts” set requires only <code>depth_ok</code>; '
      'the “same cells” set keeps the cells that pass the main cuts in all three bins.</p>')
    P(table(cells_section(res), ["bin", "donor", "fitted", "depth_ok", "noise > 1.5", "rt_rho > 0.10", "main set"],
            num=["fitted", "depth_ok", "noise > 1.5", "rt_rho > 0.10", "main set"]))
    P(method("<b>What is counted.</b> Cells fitted by <code>pureSamples_absCN</code> (≥ 20,000 reads), those passing <code>depth_ok</code>, "
             "how many of these exceed each cut (a cell can exceed both), and those left in the main set. "
             "<b>Limitation.</b> <code>depth_ok</code> depends on the batch (instrument, read length, assay). <b>Source:</b> <code>results/pilot_cn/qc/cn_qc_&lt;bin&gt;kb.tsv</code> (script <code>scripts/03_cn_qc/cn_qc_metrics.R</code>)."))

    P("""<h2 id="events">4 · Breakends and events</h2>
<ol><li><b>States.</b> On the QDNAseq segments (CBS, <code>undo.SD = 1</code>), segment ratio = value / median of the cell's segments; gain if &gt; 1.2, loss if &lt; 0.8 (permissive rule); consecutive segments with the same state are merged (length-weighted ratio) and merged segments &lt; 1 Mb are dropped.</li>
<li><b>Breakend.</b> Boundary between consecutive retained segments of a chromosome with different states. Interval [lo, hi] = from the start of the last bin of the left segment to the end of the first bin of the right one, plus any gap of filtered bins.</li>
<li><b>Deviation and stratum.</b> Deviation = max |ratio − 1| of the two flanks. E1 0.2–0.3; E2 0.3–0.4; E3 0.4–0.6; E4 ≥ 0.6. The main definition (E2–E4) is the paper's 1.3 / 0.7 threshold.</li>
<li><b>Direction.</b> Loss or gain from the flank with the larger deviation; “mixed” when one flank is a gain and the other a loss.</li>
<li><b>Type.</b> Arm-level if [lo, hi] overlaps an approximate pericentromeric block ± 1 Mb (centromere and adjacent heterochromatin, UCSC acen/gvar/stalk bands; whole p arm of acrocentrics); interstitial otherwise. The tests use interstitial events only.</li>
<li><b>Event.</b> Within each donor, breakends of the same chromosome, direction and type linked within ≤ 2 Mb (single linkage on sorted positions). Event position, lo, hi and deviation = medians of its breakends.</li></ol>""")
    P(table(events_section(res), ["bin", "type", "direction", "breakends", "events", "E1", "E2", "E3", "E4", "in 1 cell"],
            num=["breakends", "events", "E1", "E2", "E3", "E4", "in 1 cell"]))
    P(method("<b>What is counted.</b> Main set. “breakends” = breakends of all cells (a clone contributes one per cell); “events” = after collapsing within the donor; "
             "strata and “in 1 cell” are counted on events. <b>Limitation.</b> Linking at 2 Mb can merge two nearby independent breaks of the same donor. <b>Source:</b> <code>results/enrichment/pilot/main/events_&lt;bin&gt;kb.tsv</code>, <code>results/enrichment/pilot/main/breakends_&lt;bin&gt;kb.tsv</code> (script <code>scripts/06_enrichment/build_event_table.py</code>)."))

    cat_all = pd.concat([cat, der], ignore_index=True)
    cat_t = cat_all.assign(n=cat_all.n.map(lambda v: fmt(int(v))), total=cat_all.total_mb.map(lambda v: fmt(v, 1)), med=cat_all.median_mb.map(fmt))
    P('<h2 id="features">5 · Site sets</h2><p>CFS catalogues harmonised to hg19 (<code>scripts/05_features/build_cfs_catalogs.py</code>), '
      'derived catalogues (<code>build_breast_predicted_cfs.py</code>, <code>build_cfs_consensus.py</code>), '
      'and two large-gene sets from refGene (<code>build_gene_features.py</code>; span = union of overlapping transcripts).</p>')
    P("<h3>5.1 How each set is built</h3><p>One row per set as it appears in the test tables (<code>&lt;catalog&gt;:&lt;cell type&gt;</code>, "
      "or <code>:all</code> for the union of a catalogue's cell types). All coordinates are hg19; hg38/GRCh38 sources are converted with the UCSC "
      "hg38ToHg19 chain, lifting each end separately and dropping an interval that changes chromosome or whose length changes by more than 20%. "
      "Every script is in <code>scripts/05_features/</code> and documents its inputs and rules in its header.</p>")
    cm = pd.DataFrame(CATALOG_METHODS, columns=["set(s)", "source", "assay / what is measured", "cell type", "one site is",
                                                 "processing", "script"])
    P('<div class="tablebox wrap">' + table(cm, list(cm.columns))[len('<div class="tablebox">'):])
    P(method("<b>Why three kinds of sets.</b> Measured catalogues (MiDAS, replication delay, FISH, cytogenetics, CNV hotspots, cancer SVs) each capture "
             "one facet of fragility in specific cell types. Derived sets combine them (consensus) or predict fragility in breast epithelium from "
             "known determinants (breast_predicted). Large-gene sets are controls: most known CFS lie in very large genes, so an enrichment in CFS "
             "should be compared with an enrichment in large genes as such; test B also adjusts for the fraction of each window in genes &gt; 300 kb. "
             "<b>Limitation.</b> The gene sets ignore expression; breast_predicted is their expression-aware counterpart. <b>Source:</b> <code>data/features/cfs/cfs_catalogs_hg19.tsv</code>, <code>data/features/cfs/cfs_derived_hg19.tsv</code>, <code>data/features/cfs/cfs_consensus_regions.tsv</code>, <code>data/features/cfs/breast_gene_accessibility.tsv</code>, <code>data/features/genes/large_genes_*.tsv</code> (script <code>scripts/05_features/*.py</code>)."))
    P("<h3>5.2 Size of each set</h3>")
    P(table(cat_t.rename(columns={"cell_type": "cell type / tier", "total": "total Mb", "med": "median Mb"}),
            ["catalog", "cell type / tier", "n", "total Mb", "median Mb"], num=["n", "total Mb", "median Mb"]))
    P(method(f"<b>What is counted.</b> Sites per set after harmonisation (overlapping sites within a set are counted separately, so total Mb can exceed the covered genome). "
             f"Genes: &gt; 300 kb {fmt(len(g3))} loci ({fmt(g3.span_kb.sum() / 1e3, 0)} Mb); &gt; 1 Mb {fmt(len(g1))} loci ({fmt(g1.span_kb.sum() / 1e3, 0)} Mb). <b>Source:</b> <code>data/features/cfs/cfs_catalogs_summary.tsv</code>, <code>data/features/cfs/cfs_derived_hg19.tsv</code>, <code>data/features/genes/large_genes_*.tsv</code> (script <code>build_cfs_catalogs.py, build_cfs_consensus.py, build_gene_features.py</code>)."))

    P("""<h2 id="testa">6 · Test A: positional permutation</h2>
<p><b>Null.</b> Each interstitial event is re-placed on its chromosome, centred on a usable bin of the same replication-timing quintile (RT of the interval; quintiles of the usable bins of the same size), with the same width and without overlapping pericentromeric blocks ± 1 Mb. Because candidate placements are equiprobable, each event's null overlap probability p<sub>i</sub> is computed exactly by enumerating all its candidate placements.</p>
<ul><li>Observed O = events whose [lo, hi] overlaps the set; expected E = Σ p<sub>i</sub>; fold = O / E.</li>
<li>Enrichment p = fraction of 20,000 simulated sums of independent Bernoulli(p<sub>i</sub>) ≥ O (equivalent to 20,000 joint permutations).</li>
<li>Fold 95% CI: 2,000 bootstrap resamples of events (o<sub>i</sub>, p<sub>i</sub>).</li>
<li>q: BH within each set over 3 bins × 3 directions × 5 strata.</li></ul>""")
    P("<h3>6.1 Null calibration: shuffled sets</h3><p>Each set was shuffled 20 times by moving every site to a random usable bin of the same chromosome and RT quintile, keeping its width; "
      "test A (main stratum) was repeated on each version.</p>")
    P(table(shuffle_rows(res), ["set", "bin", "direction", "shuffles", "median", "2.5%", "97.5%", "observed fold"],
            num=["shuffles", "median", "2.5%", "97.5%", "observed fold"]))
    P(method("<b>What is counted.</b> Distribution of the fold over 20 shuffled versions of each set (interstitial events, E2–E4), and the fold of the real set. "
             "<b>Limitation.</b> 20 shuffles give a coarse estimate of the extreme percentiles. <b>Source:</b> <code>results/enrichment/pilot/main/test_a_shuffle_controls.tsv</code>, <code>results/enrichment/pilot/main/test_a.tsv</code> (script <code>scripts/06_enrichment/test_a_permutation.py</code>)."))
    P("<h3>6.2 Results: interstitial losses (E2–E4)</h3>")
    P(table(test_a_rows(ta, "loss"), ["set", "bin", "events", "obs.", "exp.", "fold [95% CI]", "p (enrich.)", "q"],
            num=["events", "obs.", "exp.", "p (enrich.)", "q"]))
    P(method("<b>Source:</b> <code>results/enrichment/pilot/main/test_a.tsv</code> (script <code>scripts/06_enrichment/test_a_permutation.py</code>). Rows: direction = loss, stratum = main."))
    P("<h3>6.3 Results: interstitial gains (E2–E4)</h3>")
    P(table(test_a_rows(ta, "gain"), ["set", "bin", "events", "obs.", "exp.", "fold [95% CI]", "p (enrich.)", "q"],
            num=["events", "obs.", "exp.", "p (enrich.)", "q"]))
    P(method("<b>Source:</b> <code>results/enrichment/pilot/main/test_a.tsv</code> (script <code>scripts/06_enrichment/test_a_permutation.py</code>). Rows: direction = gain, stratum = main."))
    P(fig("fig_a_forest", "<b>What is shown.</b> Test A fold (point) with bootstrap 95% CI (line), per set and bin, for interstitial E2–E4 losses and gains. "
          "Grey bar: 2.5–97.5% percentiles of the fold over the 20 shuffled versions of the set (losses in the left panel; all directions in the right one). "
          "<b>Limitation.</b> Axis truncated at 6. <b>Source:</b> <code>results/enrichment/pilot/main/test_a.tsv</code>, <code>results/enrichment/pilot/main/test_a_shuffle_controls.tsv</code> (script <code>scripts/docs/enrichment_figures.py</code>)."))
    P("<h3>6.4 By deviation stratum</h3>")
    P(fig("fig_a_strata", "<b>What is shown.</b> Test A fold with bootstrap 95% CI per stratum (E1 0.2–0.3; E2 0.3–0.4; E3 0.4–0.6; E4 ≥ 0.6) for four sets, for losses and for all directions; sets are dodged horizontally. "
          "<b>Limitation.</b> E4 has very few events (see the events table); a fold of 0 or a very high fold in E4 reflects counts of 0–2; CIs are truncated at 8. <b>Source:</b> <code>results/enrichment/pilot/main/test_a.tsv</code> (script <code>scripts/docs/enrichment_figures.py</code>)."))
    P("<h3>6.5 Loss events overlapping HS68 MiDAS sites (500 kb)</h3>")
    P(table(hs, list(hs.columns), num=["cells"]))
    P(method("<b>What is counted.</b> Main-set interstitial E2–E4 loss events at 500 kb whose interval overlaps an HS68 MiDAS domain; the annotation is the gene/intergenic column of Macheret et al. 2020 Table S1. <b>Source:</b> <code>results/enrichment/pilot/main/events_500kb.tsv</code>, <code>data/features/cfs/cfs_catalogs_hg19.tsv</code> (script <code>scripts/docs/build_enrichment_html.py</code>)."))
    P("<h3>6.6 Sensitivity</h3><p>Fold (observed / expected) for interstitial E2–E4 losses: main set, without the noise and S-phase cuts, with the same cells in the three bins (those passing the cuts at 200 kb, 500 kb and 1 Mb), and leaving out each donor.</p>")
    sens = sensitivity_rows(res)
    P(table(sens, list(sens.columns)))
    P(method("<b>Source:</b> <code>results/enrichment/pilot/main/test_a.tsv</code>, <code>results/enrichment/pilot/nocuts/test_a.tsv</code>, <code>results/enrichment/pilot/samecells/test_a.tsv</code>, <code>results/enrichment/pilot/main/test_a_without_&lt;donor&gt;.tsv</code> (script <code>scripts/06_enrichment/test_a_permutation.py</code>)."))
    P(fig("fig_a_sensitivity", "<b>What is shown.</b> Test A fold (point) with bootstrap 95% CI (line) for E2–E4 losses for five sets in each sensitivity analysis. "
          "<b>Limitation.</b> The runs share most events, so their CIs are not independent; axis truncated at 8. <b>Source:</b> <code>results/enrichment/pilot/{main,nocuts,samecells}/test_a*.tsv</code> (script <code>scripts/docs/enrichment_figures.py</code>)."))

    P(f"""<h2 id="testb">7 · Test B: per-window regression</h2>
<p>Windows: usable 1 Mb bins outside pericentromeric blocks ± 1 Mb (the same for the three event bin sizes). Rows: donor × window. Response: number of interstitial events whose midpoint falls in the window.</p>
<ul><li>Covariates: overlap with the set (binary); mean RT, GC and mappability (z-scored, from the 50 kb annotation of <code>sWGS_relativeCN</code>); fraction of the window in genes &gt; 300 kb; chromosome and donor as fixed effects; offset = log(fraction of usable 50 kb windows).</li>
<li>Models: unadjusted (set + donor) and adjusted (all covariates).</li>
<li>Family: negative binomial (NB2) with alpha from the Cameron–Trivedi auxiliary regression on the Poisson fit; if alpha ≤ 0, Poisson with HC1 robust errors. In this run: {", ".join(f"{k} = {v}" for k, v in fams.items())} of {nb} models.</li>
<li>q: BH within each set and model.</li></ul>""")
    P(table(tb, list(tb.columns), num=["events", "p", "q"]))
    P(method("<b>Source:</b> <code>results/enrichment/pilot/main/test_b.tsv</code> (script <code>scripts/06_enrichment/test_b_regression.py</code>). Rows: direction = loss, stratum = main."))
    P(fig("fig_b_irr", "<b>What is shown.</b> IRR (log scale) of overlapping each set for interstitial E2–E4 losses, unadjusted (grey) and adjusted for RT, GC, mappability, genes &gt; 300 kb, chromosome and donor (red). "
          "<b>Limitation.</b> The event midpoint is used (not its width) and the overlap per window is binary. <b>Source:</b> <code>results/enrichment/pilot/main/test_b.tsv</code> (script <code>scripts/docs/enrichment_figures.py</code>)."))
    P(method("Most windows have 0 or 1 event per donor; with counts like these the estimated dispersion is below Poisson (alpha ≤ 0) and Poisson with robust errors is used."))

    tc = test_c_rows(res)
    P("""<h2 id="testc">8 · Test C: per breakend and donor</h2>
<p>Unit: each interstitial breakend (E2–E4) of each cell, without collapsing clones, so an event weighs by the number of cells carrying it. For each breakend, o = overlaps the set; p = exact null probability (same chromosome, RT quintile and width, as in test A).</p>
<ul><li>Per donor: O/E with a 95% CI from a cell bootstrap (2,000 resamples).</li>
<li>Pooled: ΣO / ΣE with a cell bootstrap stratified by donor.</li>
<li>Heterogeneity: X² = Σ (O<sub>d</sub> − R·E<sub>d</sub>)² / (R·E<sub>d</sub>), R = pooled O/E, df = donors − 1 (donors with E<sub>d</sub> &lt; 1 left out).</li>
<li><b>Change from the plan.</b> The plan was a mixed logistic model with a donor random effect; with 5 donors that variance is not estimable, so per-donor ratios and a heterogeneity test are reported instead.</li></ul>""")
    P(table(tc, list(tc.columns), num=["cells", "breakends", "obs.", "exp.", "heterogeneity p"]))
    P(method("<b>Source:</b> <code>results/enrichment/pilot/main/test_c.tsv</code> (script <code>scripts/06_enrichment/test_c_per_donor.py</code>)."))
    P(fig("fig_c_donors", "<b>What is shown.</b> Breakend-level O/E (E2–E4 losses) per donor (colour) and pooled (black), with cell-bootstrap 95% CI. "
          "<b>Limitation.</b> Donors with few breakends have very wide CIs or CIs collapsed at 0. <b>Source:</b> <code>results/enrichment/pilot/main/test_c.tsv</code> (script <code>scripts/docs/enrichment_figures.py</code>)."))

    P('<h2 id="regions">9 · Regional views</h2>')
    P(fig("fig_regions", "<b>What is shown.</b> [lo, hi] intervals of E2–E4 loss breakends (each line is a breakend of one main-set cell, coloured by donor) at 500 kb in three regions, "
          "with HS68 and U2OS+HeLa MiDAS domains (Macheret 2020), HumCFS bands, genes &gt; 300 kb and replication timing per bin (green, up = early). <b>Source:</b> <code>results/enrichment/pilot/main/breakends_500kb.tsv</code>, <code>results/enrichment/pilot/main/usable_bins_rt_500kb.tsv</code>, <code>data/features/cfs/cfs_catalogs_hg19.tsv</code>, <code>data/features/genes/large_genes_300kb.tsv</code> (script <code>scripts/docs/enrichment_figures.py</code>)."))

    P("""<h2 id="limitations">10 · Method limitations</h2><ul>
<li>5-donor pilot; event counts range from tens to hundreds depending on the subset.</li>
<li>Resolution: breakend intervals span 0.4 to several Mb, the order of the site sizes.</li>
<li>Replication timing comes from the <code>sWGS_relativeCN</code> pipeline; its cell type is not documented.</li>
<li>Measured CFS catalogues come from other cell types (mostly tumour lines, lymphocytes and fibroblasts); the breast-predicted set is a prediction from accessibility and replication timing, not a fragility measurement.</li>
<li>Clonal collapse links within 2 Mb and requires the same direction and type; it can merge nearby independent events or split one event with different intervals.</li>
<li>Only breaks that leave a copy-number change are detected.</li>
<li>No primary hypothesis: q-values are computed per family and do not control the error across sets.</li></ul>""")

    P("""<h2 id="repro">11 · Reproducibility</h2><pre>
bash scripts/06_enrichment/run_all.sh
# = build_event_table.py (main, nocuts, samecells) -> test_a_permutation.py (main with shuffles, nocuts,
#   samecells, leave-one-donor-out) -> test_b_regression.py -> test_c_per_donor.py -> build_master_tables.py
#   -> funnel.py -> scripts/docs/enrichment_figures.py -> scripts/docs/build_enrichment_html.py
# catalogues: scripts/05_features/build_cfs_catalogs.py, build_breast_predicted_cfs.py, build_cfs_consensus.py</pre>
<p>Seed 1 in every random step. Full tables: <code>results/enrichment/pilot/{main,nocuts,samecells}/test_a*.tsv</code>, <code>test_b.tsv</code>, <code>test_c.tsv</code>, <code>master_*.tsv</code>.</p></main>""")

    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/enrichment_pilot.html", "w") as f:
        f.write("\n".join(parts))
    print("wrote", f"{a.out}/enrichment_pilot.html")


if __name__ == "__main__":
    main()

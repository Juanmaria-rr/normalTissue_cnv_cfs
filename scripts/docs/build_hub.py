#!/usr/bin/env python3
"""Build docs/hub/: one entry page that gathers every HTML report of the project.

Each report is copied into its own folder with its figures. Pages written as fragments (no
<!doctype>, the format the artifact viewer wraps) get a full HTML skeleton. Every page gets a
link back to the index. Notebooks are exported to HTML with nbconvert (outputs as last
executed). Local-only material (PDFs, Excel and TSV tables under results/, git-ignored) is
listed with relative paths that work when the hub is opened from disk.

Usage: python3 scripts/docs/build_hub.py --out docs/hub
"""
import argparse
import html
import os
import re
import shutil
import subprocess
from datetime import datetime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

REPORTS = [
    dict(slug="technical-report", src="docs/technical_report/cnv_cfs_technical_report.html",
         assets="docs/technical_report/*.png", title="Technical report",
         desc="Full methodology from data to breakends, worked numeric examples, decisions and pending steps."),
    dict(slug="enrichment", src="docs/enrichment/enrichment_pilot.html", assets="docs/enrichment/*.png",
         title="Breakend enrichment at CFS (pilot)",
         desc="Thresholds and data funnel, catalogues, tests A–C, sensitivity and regional views; generated from the result tables."),
    dict(slug="tests-explained", src="docs/enrichment/tests_explainer.html", assets=None,
         title="Tests A, B and C explained",
         desc="Teaching walk-through of the three enrichment tests on real pilot data: what each computes, what it answers, its limitations."),
    dict(slug="qc-levels-1-2", src="docs/qc_levels_1_2/qc_levels_1_2.html", assets=None,
         title="scDNA-seq QC: levels 1 and 2",
         desc="Document of 29 October 2026 (reads, alignment, relative and absolute CN). Some later decisions are not reflected; see the technical report."),
    dict(slug="aws-costs", src="docs/aws_costs/aws_costs.html", assets=None, title="AWS costs",
         desc="Daily account cost by service and by project (estimated), up to 28 October 2026."),
]
NOTEBOOKS = [
    dict(slug="notebook-calibration", src="notebooks/01_calibration_exploration.ipynb",
         title="Notebook 01 · QC calibration", desc="Noise, S phase, ratio threshold and minimum segment, per donor, with cell profiles."),
    dict(slug="notebook-enrichment", src="notebooks/02_enrichment_exploration.ipynb",
         title="Notebook 02 · master tables and test reproduction",
         desc="Cell, breakend and event tables; cell and region viewers; step-by-step reproduction of tests A–C."),
]
LOCAL = [
    ("Breakend examples (PDF)", "results/pilot_cn/review/example_breakends.pdf"),
    ("Review of 30 random cells per donor (PDF)", "results/pilot_cn/review/"),
    ("Master tables (Excel)", "results/enrichment/pilot/master_tables.xlsx"),
    ("Master tables and results (TSV)", "results/enrichment/pilot/"),
    ("Harmonised CFS catalogues (hg19)", "data/features/cfs/"),
    ("Chronological log (WORKLOG.md)", "WORKLOG.md"),
    ("Project status (README.md)", "README.md"),
]
# Project structure guide: (path, in git, what it holds). Kept by hand; update when folders change.
STRUCTURE = [
    ("README.md", "yes", "Short orientation: purpose, status, documents, layout, how to reproduce, environment, known limitations."),
    ("WORKLOG.md", "no (local only)", "Chronological log (newest first): logic, changes and results of every working session, and every decision with its date."),
    ("scripts/00_metadata/", "yes", "ENA and paper metadata parsing (cells.tsv), data inventory, usable-bin counts, pilot cell selection."),
    ("scripts/01_staging/", "yes", "EC2 helper scripts that copied the SRA files from the public us-east-1 bucket to eu-west-1 (see docs/AWS_RUNBOOK.md)."),
    ("scripts/02_align/", "yes", "HealthOmics samplesheets, level-1 QC summary (reads, mapping, duplicates), FASTQ duplicate check."),
    ("scripts/03_cn_qc/", "yes", "Level-2 QC metrics (depth_ok, noise_excess, rt_rho, fit flags), ratio-based breakend calling, calibration tables, example-cell selection, relative-vs-absolute review PDFs, absolute-fit stability."),
    ("scripts/04_curation/", "yes", "Export of relative CN for the curation app, batch building with blinding, unblinding of curated fits (manual-fit audit)."),
    ("scripts/05_features/", "yes", "Feature catalogues for the enrichment tests: CFS catalogues (8 sources, hg19), breast-predicted CFS from scATAC, consensus / A+C / union, large genes, power estimate."),
    ("scripts/06_enrichment/", "yes", "Event tables (breakends, clonal collapse, strata), tests A (positional permutation), B (per-window regression) and C (per donor), master tables, data funnel; run_all.sh rebuilds everything."),
    ("scripts/aws/", "yes", "AWS cost report by service and project, and its HTML page."),
    ("scripts/docs/", "yes", "Builders of the documents and figures: enrichment document, concept figures, QC document, this hub."),
    ("workflows/", "yes", "Nextflow workflows run on AWS HealthOmics: build_reference (bwa index), align (FASTQ to BAM), copynumber (QDNAseq relative CN + absolute fits)."),
    ("containers/", "yes", "Dockerfiles of the align and copynumber images (tools pinned to the group pipelines), build and push to ECR."),
    ("external/", "yes (submodules)", "Group pipelines pinned to a commit: scDNAseq, sWGS_relativeCN (QDNAseq, bin and RT annotations), fitting_absCN, CopyNumberFittingShiny."),
    ("curation/", "code yes; data no", "Vendored curation app (Shiny for Python) for the manual-fit audit, with blinding keys; batches in curation/data are not versioned."),
    ("enrichment_script/", "yes", "The original enrichment script and feature files as received (kept for reference; replaced by scripts/06_enrichment)."),
    ("notebooks/", "yes", "Executed exploration notebooks (01 calibration, 02 enrichment master tables and test reproduction) and the scripts that build them."),
    ("docs/technical_report/", "yes", "Technical report (methodology, worked examples, decisions) and its figures."),
    ("docs/enrichment/", "yes", "Enrichment methods and results document (generated) and its figures."),
    ("docs/qc_levels_1_2/", "yes", "QC levels 1–2 document (29 October 2026) with its data tables and figures."),
    ("docs/aws_costs/", "yes", "AWS cost page and its daily tables."),
    ("docs/AWS_RUNBOOK.md", "yes", "Manual AWS procedure: commands, expected results, failure modes and cost."),
    ("docs/hub/", "no (generated)", "This documents hub (scripts/docs/build_hub.py)."),
    ("data/metadata/", "partly", "cells.tsv (one row per ENA run), donors, usable bins per bin size, data inventory, paper tables."),
    ("data/pilot/", "partly", "Pilot cell lists, samplesheets, alignment QC tables, staging manifests, duplicate check."),
    ("data/features/cfs/", "tables yes; raw no", "Harmonised CFS catalogues in hg19 (cfs_catalogs_hg19.tsv), derived sets (cfs_derived_hg19.tsv), consensus regions, breast gene accessibility; raw/ holds the downloaded supplements."),
    ("data/features/genes/", "lists yes; raw no", "Large-gene lists from refGene hg19 (> 300 kb, > 1 Mb)."),
    ("data/features/atac/", "no", "GSE272504 scATAC fragment files (1.8 GB), the authors' barcode map and the insertion-count cache."),
    ("results/pilot_align/", "no", "Alignment QC outputs of the pilot."),
    ("results/pilot_cn/relCN/, absCN/", "no", "Copy-number objects downloaded from HealthOmics (QDNAseq .rds, absolute fits per cell, failed_cells.tsv)."),
    ("results/pilot_cn/qc/", "no", "Per-cell QC tables (cn_qc_<bin>kb.tsv) and breakend calls (segments, breakends, cells) per bin."),
    ("results/pilot_cn/calibration/, fit_stability/", "no", "Calibration tables and absolute-fit stability tables."),
    ("results/pilot_cn/notebook/, review/", "no", "Exported profiles of example cells; review and example-breakend PDFs."),
    ("results/enrichment/pilot/", "no", "Event tables per cell set (main, nocuts, samecells), test results (test_a*, test_b, test_c), master tables (TSV + Excel), funnel tables, pipeline log."),
    ("s3://<PROJECT_BUCKET>/cnv_cfs/", "S3", "Staged SRA files, hg19 reference and bwa index, HealthOmics run outputs (BAMs, relCN, absCN) in eu-west-1, account <AWS_ACCOUNT_ID>."),
]

SKELETON = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">'
            '<title>{title}</title></head><body>\n{back}\n{body}\n</body></html>\n')
BACK = ('<p style="font:600 13px/1.4 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:12px 16px">'
        '<a href="../index.html" style="color:#0f6e73">← All documents</a></p>')


def copy_report(r, out):
    d = os.path.join(out, r["slug"])
    os.makedirs(d, exist_ok=True)
    txt = open(os.path.join(ROOT, r["src"]), encoding="utf-8").read()
    if txt.lstrip().lower().startswith("<!doctype"):
        txt = re.sub(r"(<body[^>]*>)", r"\1" + BACK.replace("\\", "\\\\"), txt, count=1, flags=re.I)
    else:
        m = re.search(r"<title>(.*?)</title>", txt, flags=re.S)
        txt = SKELETON.format(title=m.group(1) if m else r["title"], back=BACK, body=txt)
    open(os.path.join(d, "index.html"), "w", encoding="utf-8").write(txt)
    if r.get("assets"):
        import glob
        for f in glob.glob(os.path.join(ROOT, r["assets"])):
            shutil.copy(f, d)
    return datetime.fromtimestamp(os.path.getmtime(os.path.join(ROOT, r["src"]))).strftime("%Y-%m-%d")


def export_notebook(n, out):
    d = os.path.join(out, n["slug"])
    os.makedirs(d, exist_ok=True)
    subprocess.run(["python3", "-m", "nbconvert", "--to", "html", "--output-dir", d, "--output", "index",
                    os.path.join(ROOT, n["src"])], check=True, capture_output=True)
    f = os.path.join(d, "index.html")
    txt = open(f, encoding="utf-8").read()
    txt = re.sub(r"(<body[^>]*>)", r"\1" + BACK.replace("\\", "\\\\"), txt, count=1, flags=re.I)
    open(f, "w", encoding="utf-8").write(txt)
    return datetime.fromtimestamp(os.path.getmtime(os.path.join(ROOT, n["src"]))).strftime("%Y-%m-%d")


INDEX_CSS = """
:root{--bg:#f7f8f9;--surface:#fff;--fg:#1b232b;--muted:#56636f;--rule:#d6dde3;--accent:#0f6e73}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#11161b;--surface:#171e24;--fg:#e3e8ec;--muted:#9aa7b2;--rule:#2c3640;--accent:#5cc3c7;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#11161b;--surface:#171e24;--fg:#e3e8ec;--muted:#9aa7b2;--rule:#2c3640;--accent:#5cc3c7;color-scheme:dark}
body{background:var(--bg);color:var(--fg);font:15px/1.55 "Public Sans",-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;
padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px));padding-block:40px 80px}
main{max-width:880px;margin:0 auto}h1{font:600 2rem/1.2 "Newsreader",Georgia,serif;margin:0 0 .3rem}
h2{font:600 1.25rem/1.3 "Newsreader",Georgia,serif;margin:2.2rem 0 .8rem}
.meta{color:var(--muted);font-size:.9rem}.eyebrow{font:600 .72rem/1 "JetBrains Mono",ui-monospace,monospace;letter-spacing:.08em;text-transform:uppercase;color:var(--accent)}
.list{display:grid;gap:10px}.item{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:4px 16px;padding:12px 14px;background:var(--surface);border:1px solid var(--rule);border-radius:6px}
.item a{color:var(--accent);font-weight:600;text-decoration:none}.item a:hover,.item a:focus-visible{text-decoration:underline}
.item p{grid-column:1/-1;margin:0;color:var(--muted);font-size:.9rem}.date{font:.8rem "JetBrains Mono",ui-monospace,monospace;color:var(--muted)}
code{font:.85em "JetBrains Mono",ui-monospace,monospace}
.tbox{overflow-x:auto;border:1px solid var(--rule);border-radius:6px;background:var(--surface)}
.tbox table{border-collapse:collapse;width:100%;font-size:.86rem}.tbox th,.tbox td{padding:6px 10px;border-bottom:1px solid var(--rule);text-align:left;vertical-align:top}
.tbox th{font-weight:600}.tbox td:first-child{white-space:nowrap}.tbox td:nth-child(2){white-space:nowrap;color:var(--muted)}ul{padding-left:1.1rem}li+li{margin-top:.3rem}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = os.path.join(ROOT, a.out)
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    rep = [(r, copy_report(r, out)) for r in REPORTS]
    nbs = [(n, export_notebook(n, out)) for n in NOTEBOOKS]

    def item(x, date):
        return (f'<div class="item"><a href="{x["slug"]}/index.html">{html.escape(x["title"])}</a>'
                f'<span class="date">{date}</span><p>{html.escape(x["desc"])}</p></div>')
    rel = os.path.relpath(ROOT, out)
    struct = "".join(f"<tr><td><code>{html.escape(p)}</code></td><td>{html.escape(g)}</td><td>{html.escape(t)}</td></tr>"
                     for p, g, t in STRUCTURE)
    local = "".join(f'<li><a href="{rel}/{p}">{html.escape(t)}</a> · <code>{p}</code></li>' for t, p in LOCAL)
    page = f"""<title>CNV_cfs documents</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,600&family=Public+Sans:wght@400;600&family=JetBrains+Mono&display=swap">
<style>{INDEX_CSS}</style>
<main>
<span class="eyebrow">CNV_cfs · Lin et al. 2024 (PRJNA1003661) · hg19</span>
<h1>Project documents</h1>
<p class="meta">Copy-number breakends in normal breast epithelial cells and their enrichment at common fragile sites. Index generated on {datetime.now().strftime("%Y-%m-%d")} by <code>scripts/docs/build_hub.py</code>; each document's date is its last modification.</p>
<h2>Reports</h2><div class="list">{"".join(item(r, d) for r, d in rep)}</div>
<h2>Notebooks (exported to HTML)</h2><div class="list">{"".join(item(n, d) for n, d in nbs)}</div>
<h2>Project structure</h2>
<p class="meta">What each folder holds. “In git” says whether its content is versioned in the repository; the rest is regenerated by the scripts or downloaded.</p>
<div class="tbox"><table><thead><tr><th>Path</th><th>In git</th><th>Contents</th></tr></thead><tbody>{struct}</tbody></table></div>
<h2>Local-only material</h2>
<p class="meta">Files outside git; these links open only from a local copy of the repository (not in the published version).</p>
<ul>{local}</ul>
</main>"""
    open(os.path.join(out, "index.html"), "w", encoding="utf-8").write(page)
    n = sum(len(fs) for _, _, fs in os.walk(out))
    print(f"wrote {out}/index.html and {n - 1} files")


if __name__ == "__main__":
    main()

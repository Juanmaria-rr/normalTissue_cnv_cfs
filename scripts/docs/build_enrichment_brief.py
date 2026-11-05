#!/usr/bin/env python3
"""Build docs/enrichment/enrichment_brief.html: a concise version of the enrichment document,
focused on the results of tests A, B and C, keeping only the method needed to read them.

Same inputs as build_enrichment_html.py (results/enrichment/pilot/: main/test_a.tsv,
test_a_shuffle_controls.tsv, test_b.tsv, test_c.tsv, events_<bin>kb.tsv; robustness_*.tsv;
power_scaleup.tsv; funnel_cells.tsv) and its figures (fig_a_forest, fig_top_donors, fig_b_irr,
fig_c_donors, written to the same folder by enrichment_figures.py). Every number is read from
those tables; the full document (enrichment_pilot.html) keeps every table and variant. The
site-set section reuses CATALOG_METHODS of build_enrichment_html.py, with set sizes from
data/features/cfs/cfs_catalogs_summary.tsv, cfs_derived_hg19.tsv and data/features/genes/.

Usage: python3 scripts/docs/build_enrichment_brief.py --res results/enrichment/pilot --out docs/enrichment
"""
import argparse
import html
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))
from build_enrichment_html import CATALOG_METHODS, CSS, MAIN_SETS, fig, fmt, fp, method, table  # noqa: E402
from enrichment_figures import top_row  # noqa: E402

BINS = [200, 500, 1000]
DIRS = [("loss", "losses"), ("gain", "gains"), ("all", "all")]
R = "results/enrichment/pilot"

EXTRA_CSS = """
.mx td,.mx th{text-align:right;padding:4px 8px}.mx td:first-child,.mx th:first-child{text-align:left}
.mx thead tr:first-child th{text-align:center;border-bottom:1px solid var(--rule)}
.mx td small{display:block;color:var(--muted);font-size:.72rem;line-height:1.2}
.mx td.hit{background:color-mix(in srgb,var(--accent) 14%,transparent)}.mx td.strong{font-weight:600}
.mx td.dim{color:var(--muted);opacity:.55}.mx th.sep,.mx td.sep{border-left:1px solid var(--rule)}
.keys{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin:1rem 0}
.key{background:var(--surface);border:1px solid var(--rule);border-radius:6px;padding:.6rem .8rem;min-width:0}
.key b{display:block;font:600 .72rem/1.3 var(--mono);letter-spacing:.06em;text-transform:uppercase;color:var(--accent);margin-bottom:.25rem}
.key p{margin:0;font-size:.88rem}
.legend{font-size:.78rem;color:var(--muted);margin:.2rem 0 .6rem}
"""


def esc(x):
    return html.escape(str(x))


def matrix(rows, cell, sets):
    """Sets x (direction, bin) table. cell(set, direction, bin) -> (main, sub, classes) or None."""
    h1 = "<tr><th></th>" + "".join(f'<th colspan="3" class="sep">{lab}</th>' for _, lab in DIRS) + "</tr>"
    h2 = "<tr><th>set</th>" + "".join(f'<th class="{"sep" if k == 0 else ""}">{"1 Mb" if b == 1000 else f"{b} kb"}</th>'
                                       for _ in DIRS for k, b in enumerate(BINS)) + "</tr>"
    body = []
    for s in sets:
        tds = []
        for d, _ in DIRS:
            for k, b in enumerate(BINS):
                c = cell(s, d, b)
                cls = (["sep"] if k == 0 else [])
                if c is None:
                    tds.append(f'<td class="{" ".join(cls)}">–</td>')
                    continue
                main, sub, extra = c
                tds.append(f'<td class="{" ".join(cls + extra)}">{main}<small>{sub}</small></td>')
        body.append(f"<tr><td>{esc(s)}</td>{''.join(tds)}</tr>")
    return f'<div class="tablebox"><table class="mx"><thead>{h1}{h2}</thead><tbody>{"".join(body)}</tbody></table></div>'


def flow_line(inputs, script, outputs):
    c = lambda xs: ", ".join(f"<code>{x}</code>" for x in xs)
    return f'<div class="flow"><b>Data flow.</b> <span>Input:</span> {c(inputs)} → <span>Script:</span> {c(script)} → <span>Output:</span> {c(outputs)}.</div>'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = a.res
    ta = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    tb = pd.read_csv(f"{res}/main/test_b.tsv", sep="\t")
    tc = pd.read_csv(f"{res}/main/test_c.tsv", sep="\t")
    fc = pd.read_csv(f"{res}/funnel_cells.tsv", sep="\t")
    rv = pd.read_csv(f"{res}/robustness_variants.tsv", sep="\t")
    lc = pd.read_csv(f"{res}/robustness_loco.tsv", sep="\t", dtype={"chromosome_out": str})
    top = top_row(ta)
    tam = ta[ta.stratum == "main"].set_index(["feature_set", "direction", "bin_kb"])
    tbm = tb[(tb.stratum == "main") & (tb.model == "adjusted")].set_index(["feature_set", "direction", "bin_kb"])
    tcp = tc[tc.donor == "pooled"].set_index(["feature_set", "direction", "bin_kb"])
    n_glob = int(((ta.stratum == "main") & ~ta.feature_set.isin(["genes_300kb", "genes_1mb"])).sum())

    # cells and events
    tot = fc[fc.donor == "total"].set_index(["stage", "bin_kb"]).cells
    st = lambda s, b: int(tot.loc[(s, b)])
    ev_rows = []
    for b in BINS:
        e = pd.read_csv(f"{res}/main/events_{b}kb.tsv", sep="\t")
        e = e[(e.type == "interstitial") & e.main_def]
        ev_rows.append((b, int((e.direction == "loss").sum()), int((e.direction == "gain").sum()), len(e)))

    def cell_a(s, d, b):
        if (s, d, b) not in tam.index:
            return None
        r = tam.loc[(s, d, b)]
        if not r.n_events or not np.isfinite(r.fold):
            return None
        cls = []
        if r.expected < 2:
            cls.append("dim")
        if np.isfinite(r.p_perm) and r.p_perm < 0.05:
            cls.append("hit")
        if np.isfinite(r.get("q_perm_global", np.nan)) and r.q_perm_global < 0.1:
            cls.append("strong")
        return f"{fmt(r.fold)}", f"{int(r.observed)}/{fmt(r.expected, 1)} · p {fp(r.p_perm)}", cls

    def cell_b(s, d, b):
        if (s, d, b) not in tbm.index:
            return None
        r = tbm.loc[(s, d, b)]
        if not np.isfinite(r.irr):
            return "–", "separation", ["dim"]
        cls = (["hit"] if r.p < 0.05 else []) + (["strong"] if np.isfinite(r.q_global) and r.q_global < 0.1 else [])
        return f"{fmt(r.irr)}", f"{fmt(r.irr_lo)}–{fmt(r.irr_hi)} · p {fp(r.p)}", cls

    def cell_c(s, d, b):
        if (s, d, b) not in tcp.index:
            return None
        r = tcp.loc[(s, d, b)]
        if not np.isfinite(r.ratio):
            return None
        cls = (["hit"] if r.ratio_lo > 1 else []) + (["dim"] if r.expected < 2 else [])
        return f"{fmt(r.ratio)}", f"{fmt(r.ratio_lo)}–{fmt(r.ratio_hi)} · {int(r.observed)}/{fmt(r.expected, 1)}", cls

    sets_c = [s for s in MAIN_SETS if s in set(tc.feature_set)]
    sets_b = [s for s in MAIN_SETS if s in set(tb.feature_set)]
    tb_top = tb[(tb.feature_set == top.feature_set) & (tb.bin_kb == top.bin_kb) & (tb.direction == top.direction) & (tb.stratum == "main")].set_index("model")
    tc_top = tc[(tc.feature_set == top.feature_set) & (tc.bin_kb == top.bin_kb) & (tc.direction == top.direction)].set_index("donor")
    rv_top = rv[(rv.feature_set == top.feature_set) & (rv.direction == top.direction)]
    lc_top = lc[(lc.feature_set == top.feature_set) & (lc.direction == top.direction) & (lc.bin_kb == top.bin_kb)]
    lo_c = lc_top.loc[lc_top.fold.idxmin()]
    donors_top = tc_top.drop(index="pooled", errors="ignore")
    tdir = {"loss": "losses", "gain": "gains", "all": "events (all directions)"}[top.direction]

    parts = []
    P = parts.append
    P(f"""<title>CFS enrichment brief</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,600&family=Public+Sans:wght@400;600&family=JetBrains+Mono&display=swap">
<style>{CSS}{EXTRA_CSS}</style><main>
<span class="eyebrow">CNV_cfs · pilot (5 donors) · hg19 · concise version</span>
<h1>Breakends at common fragile sites: results of tests A, B and C</h1>
<p class="meta">Generated on {date.today().isoformat()} by <code>scripts/docs/build_enrichment_brief.py</code> from <code>{R}/</code>. Every number is read from the result tables.
The full document (<code>docs/enrichment/enrichment_pilot.html</code>) has every table, stratum and sensitivity variant; the tests explainer works through each test on one example.</p>""")

    P(f"""<h2 id="question">1 · Question and data</h2>
<p>Do interstitial copy-number breakends of normal breast epithelial cells (Lin et al. 2024, PRJNA1003661) overlap common fragile sites (CFS) more often than expected under a null that keeps chromosome and replication timing? Three tests answer it with different units and nulls; all are exploratory (no primary hypothesis was declared before the pilot).</p>
<div class="keys">
<div class="key"><b>Cells</b><p>Main set: {st('8 rt_rho <= 0.10', 200):,} / {st('8 rt_rho <= 0.10', 500):,} / {st('8 rt_rho <= 0.10', 1000):,} cells at 200 kb / 500 kb / 1 Mb, of {st('2 imaging', 500):,} imaging-selected cells. Cuts: ≥ 15 reads per usable bin per copy (diploid), noise_excess ≤ 1.5, rt_rho ≤ 0.10 (S phase).</p></div>
<div class="key"><b>Breakend</b><p>Boundary between two segments (≥ 1 Mb of usable bins) of different state in one cell; interval [lo, hi] = one bin on each side. Main definition E2–E4: flank ratio &gt; 1.3 or &lt; 0.7 (the paper's threshold).</p></div>
<div class="key"><b>Event</b><p>Breakends of one donor with the same chromosome, direction and type within 2 Mb, collapsed into one: a clone counts once. Interstitial = interval does not touch a pericentromeric block ± 1 Mb.</p></div>
<div class="key"><b>Events tested</b><p>Interstitial E2–E4 losses / gains: {" · ".join(f"{'1 Mb' if b == 1000 else f'{b} kb'} {l} / {g}" for b, l, g, _ in ev_rows)}.</p></div>
</div>
<p><b>Site sets.</b> Eight published CFS catalogues, derived sets and two gene-size controls; see section 2.</p>
<p><b>Multiple testing.</b> Benjamini–Hochberg over all {n_glob} main-stratum rows of the CFS sets (sets × bins × directions; controls excluded): <b>q global</b> on the exact p, <b>q perm. global</b> on the permutation p. A q computed within one set only (q family) is in the full tables but does not control the error across sets.</p>""")

    # site sets: sizes per catalogue (all cell types) and per controls
    cs = pd.read_csv(f"{ROOT}/data/features/cfs/cfs_catalogs_summary.tsv", sep="\t")
    dv = pd.read_csv(f"{ROOT}/data/features/cfs/cfs_derived_hg19.tsv", sep="\t")
    size = {c: (int(g.n.sum()), float(g.total_mb.sum())) for c, g in cs.groupby("catalog")}
    for c, g in dv.groupby("catalog"):
        size[c] = (len(g), float((g.end - g.start).sum() / 1e6))
    for lab in ("300kb", "1mb"):
        g = pd.read_csv(f"{ROOT}/data/features/genes/large_genes_{lab}.tsv", sep="\t")
        size[f"genes_{lab}"] = (len(g), float(g.span_kb.sum() / 1e3))
    kind = lambda c: "control" if c.startswith("genes_") else ("derived" if c in set(dv.catalog) else "measured")
    rows = []
    for sets_, src, assay, ct, site, proc, script in CATALOG_METHODS:
        c = sets_.split(":")[0].split(" ")[0]
        n, mbs = size.get(c, (None, None))
        rows.append({"set(s)": sets_, "kind": kind(c), "source": src, "what is measured": assay, "cell type": ct, "one site is": site,
                     "processing": proc, "sites (Mb)": f"{n:,} ({mbs:,.0f})" if n is not None else "–"})
    cm = pd.DataFrame(rows)
    P(f"""<h2 id="sets">2 · Site sets (CFS catalogues)</h2>
<p><b>Names.</b> Each set is written <code>&lt;catalogue&gt;:&lt;cell type&gt;</code>: the catalogue is the source study (first author and year) and the cell type the line in which the sites were measured; <code>:all</code> is the union of all the catalogue's cell types (for example <code>macheret2020:all</code> = MiDAS domains of HS68, HeLa and U2OS; <code>letallec2013:MCF10A</code> = sites mapped by FISH in MCF10A breast epithelial cells). Gene controls have no cell type.</p>
<p><b>Three kinds.</b> <i>Measured</i> catalogues each capture one facet of fragility in specific cell types (MiDAS = DNA synthesis in mitosis after mild aphidicolin; replication delay; FISH or cytogenetic breaks; CNV or cancer-SV hotspots). <i>Derived</i> sets combine them (consensus) or predict fragility in breast epithelium (breast_predicted). <i>Controls</i> are genes &gt; 300 kb and &gt; 1 Mb: most CFS lie in very large genes, so an enrichment at CFS should be read against an enrichment at large genes as such. All coordinates are hg19 (hg38 sources lifted with the UCSC chain).</p>""")
    P('<div class="tablebox wrap">' + table(cm, list(cm.columns))[len('<div class="tablebox">'):])
    P(method("<b>What is counted.</b> Sites per catalogue over all its cell types (overlapping sites of different cell types counted separately, so Mb can exceed the genome they cover). "
             "<b>Limitation.</b> Measured catalogues come from other cell types (mostly tumour lines, lymphocytes, fibroblasts); only Le Tallec's MCF10A sites (3) are breast epithelium, and breast_predicted is a prediction, not a measurement."))
    P(flow_line(["published supplements (data/features/cfs/raw/)", "UCSC refGene hg19", "GSE272504 scATAC"],
                ["scripts/05_features/build_cfs_catalogs.py", "build_gene_features.py", "build_breast_predicted_cfs.py", "build_cfs_consensus.py"],
                ["data/features/cfs/cfs_catalogs_hg19.tsv", "cfs_derived_hg19.tsv", "data/features/genes/large_genes_*.tsv"]))

    P(f"""<h2 id="testa">3 · Test A: events re-placed on the genome</h2>
<ul><li><b>Unit:</b> one event. <b>Observed O:</b> events whose [lo, hi] overlaps the set.</li>
<li><b>Null:</b> each event is moved, with its width, to every usable-bin start of its own chromosome where the placement has the event's replication-timing quintile and does not touch a pericentromere; p<sub>i</sub> = fraction of those placements that overlap the set (exact). <b>E</b> = Σ p<sub>i</sub>; <b>fold</b> = O / E.</li>
<li><b>p perm.</b> (reference): the set is shuffled 5,000 times (each site to a random usable bin of the same chromosome and RT quintile) and the fold recomputed; it keeps the clustering of events and sites. The exact Poisson-binomial p assumes independent events and is smaller.</li></ul>
<p class="legend">Cell: fold; below, observed / expected events and p perm. Tinted: p perm. &lt; 0.05. Bold: q perm. global &lt; 0.1. Faded: expected &lt; 2 events (fold dominated by counting noise).</p>""")
    P(matrix(None, cell_a, MAIN_SETS))
    P(method("<b>What is counted.</b> Main-set interstitial E2–E4 events per set, bin and direction (stratum main of <code>test_a.tsv</code>). "
             "<b>How.</b> As above; fold CIs (event and donor bootstrap) are in the table file and the forest plot. "
             "<b>Limitation.</b> Only chromosome and RT are matched; intervals (≥ 2 bins) are often wider than the sites, so overlap measures proximity at bin scale."))
    P(fig("fig_a_forest", "<b>What is shown.</b> Test A fold with event-bootstrap 95% CI per set and bin (colours); grey bar: 2.5–97.5% of the fold over the 5,000 shuffled sets; hollow: expected &lt; 2. "
          "<b>Limitation.</b> Axis truncated at 6; the bootstrap CI assumes independent events."))
    P(flow_line([f"{R}/main/events_&lt;bin&gt;kb.tsv", f"{R}/main/usable_bins_rt_&lt;bin&gt;kb.tsv", "data/features/cfs/*.tsv", "data/features/genes/*.tsv"],
                ["scripts/06_enrichment/test_a_permutation.py"], [f"{R}/main/test_a.tsv", f"{R}/main/test_a_shuffle_controls.tsv"]))

    sens = rv_top.pivot_table(index="variant", columns="bin_kb", values="fold")
    sens = sens.reindex([v for v in dict.fromkeys(rv.variant) if v in sens.index])
    srows = "".join(f"<tr><td>{esc(v)}</td>" + "".join(f'<td class="num">{fmt(sens.loc[v].get(b, np.nan))}</td>' for b in BINS) + "</tr>" for v in sens.index)
    P(f"""<h3>3.1 Strongest row: {esc(top.feature_set)}, {int(top.bin_kb)} kb, {tdir}</h3>
<p>Chosen by a fixed rule (smallest p perm. among main-stratum CFS rows), so it is the pilot's best lead, not a planned test.
Fold <b>{fmt(top.fold)}</b> ({int(top.observed)} observed / {fmt(top.expected)} expected), p exact {fp(top.p_enrich)}, p perm. {fp(top.p_perm)},
q global {fp(top.q_global)}, q perm. global {fp(top.q_perm_global)}. Leaving out one chromosome at a time the fold stays between {fmt(lc_top.fold.min())} and {fmt(lc_top.fold.max())} (lowest without chr{lo_c.chromosome_out}).</p>""")
    P(fig("fig_top_donors", "<b>What is shown.</b> Observed / expected events per donor and pooled for this row, with exact Poisson 95% CI; right: O / E and number of events. "
          "<b>Limitation.</b> Per-donor expectations are below 1–2 events, so only the direction across donors is informative."))
    P(f'<div class="tablebox"><table><thead><tr><th>variant</th><th class=num>200 kb</th><th class=num>500 kb</th><th class=num>1 Mb</th></tr></thead><tbody>{srows}</tbody></table></div>')
    P(method("<b>What is counted.</b> Test A fold of this set and direction in each robustness variant: main; nocuts (no noise / S-phase cut); samecells (cells passing in all bins); "
             "cells500 (500 kb main cells called at every bin: bins differ only in resolution); refined (500 kb breakends refined to 200 kb in the same cell; 200 kb column); "
             "linkage 1 / 5 Mb; collapse on E2–E4 only; leave one donor out. "
             "<b>Limitation.</b> Variants share most events; they test stability, not independent replication."))
    P(flow_line([f"{R}/&lt;variant&gt;/test_a.tsv", f"{R}/main/test_a_event_probs.tsv"], ["scripts/06_enrichment/robustness.py"],
                [f"{R}/robustness_variants.tsv", f"{R}/robustness_loco.tsv"]))

    gb = tb_top.loc["adjusted_genes300"] if "adjusted_genes300" in tb_top.index else None
    P(f"""<h2 id="testb">4 · Test B: regression on 1 Mb windows</h2>
<ul><li><b>Unit:</b> donor × usable 1 Mb window outside pericentromeres. Each event is assigned to the window it overlaps most (events overlapping no eligible window are excluded and counted).</li>
<li><b>Model:</b> Poisson (negative binomial if over-dispersed): events ~ overlaps set (yes/no) + RT + GC + mappability + chromosome + donor, offset log(usable fraction). <b>IRR</b> = event rate in windows of the set relative to the rest.</li>
<li>Standard error = the larger of the model-based and the chromosome-clustered one. Rows with no event in set windows have no estimate (separation).</li>
<li>Also fitted: crude (set + donor) and adjusted + fraction of the window in genes &gt; 300 kb. Large genes are probably part of the fragility mechanism, so that model is conservative.</li></ul>
<p class="legend">Cell: adjusted IRR; below, 95% CI and p. Tinted: p &lt; 0.05. Bold: q global &lt; 0.1.</p>""")
    P(matrix(None, cell_b, sets_b))
    P(method(f"<b>What is counted.</b> Adjusted model, stratum main, of <code>test_b.tsv</code>. For the strongest test A row: crude IRR {fmt(tb_top.loc['crude'].irr)} (p {fp(tb_top.loc['crude'].p)}), "
             f"adjusted {fmt(tb_top.loc['adjusted'].irr)} (p {fp(tb_top.loc['adjusted'].p)}, q global {fp(tb_top.loc['adjusted'].q_global)})"
             + (f", adjusted + genes &gt; 300 kb {fmt(gb.irr)} (p {fp(gb.p)})" if gb is not None else "") + ". "
             "<b>Limitation.</b> Overlap per window is binary and the 1 Mb grid can credit a site to the neighbouring window; counts are sparse (mostly 0–1 per donor and window)."))
    P(fig("fig_b_irr", "<b>What is shown.</b> IRR (log scale, 95% CI) per set, bin and direction: crude (grey), adjusted (red), adjusted + genes &gt; 300 kb (blue). "
          "<b>Limitation.</b> Rows with separation are not drawn."))
    P(flow_line([f"{R}/main/events_&lt;bin&gt;kb.tsv", f"{R}/main/usable_bins_rt_1000kb.tsv", "bins_with_replicationTime_50kb.tsv", "data/features/genes/large_genes_300kb.tsv"],
                ["scripts/06_enrichment/test_b_regression.py"], [f"{R}/main/test_b.tsv"]))

    P(f"""<h2 id="testc">5 · Test C: every cell's breakend, per donor</h2>
<ul><li><b>Unit:</b> one breakend of one cell, not collapsed: a clone weighs as many times as it has cells. o and p<sub>i</sub> as in test A.</li>
<li><b>Per donor and pooled:</b> O/E with a 95% CI from resampling cells (2,000 times; one-sided upper bound 3.0 / E when O = 0). Heterogeneity between donors: X² on donors with E ≥ 1.</li>
<li>It asks whether the excess is shared by donors and how clonal expansion changes it; it mixes break rate with clonal expansion.</li></ul>
<p class="legend">Cell: pooled O/E; below, 95% CI and observed / expected breakends. Tinted: CI above 1. Faded: expected &lt; 2.</p>""")
    P(matrix(None, cell_c, sets_c))
    don_txt = ", ".join(f"{d} {fmt(r.ratio)} ({int(r.observed)}/{fmt(r.expected, 1)})" for d, r in donors_top.iterrows())
    pooled = tc_top.loc["pooled"] if "pooled" in tc_top.index else None
    P(method(f"<b>What is counted.</b> Pooled rows of <code>test_c.tsv</code>. For the strongest test A row, per donor: {don_txt}"
             + (f"; heterogeneity p {fp(pooled.het_p)}" if pooled is not None and np.isfinite(pooled.het_p) else "; heterogeneity not computable (fewer than two donors with E ≥ 1)") + ". "
             "<b>Limitation.</b> One large clone (e.g. P08's gain clone) can dominate its donor and the pooled ratio; the cell bootstrap treats cells of a clone as independent."))
    P(fig("fig_c_donors", "<b>What is shown.</b> Breakend-level O/E per donor (colour) and pooled (black), 95% cell-bootstrap CI, for four sets; rows of panels = losses, gains, all. "
          "<b>Limitation.</b> Donors with few breakends have very wide CIs."))
    P(flow_line([f"{R}/main/breakends_&lt;bin&gt;kb.tsv", f"{R}/main/usable_bins_rt_&lt;bin&gt;kb.tsv", "data/features/cfs/*.tsv"],
                ["scripts/06_enrichment/test_c_per_donor.py"], [f"{R}/main/test_c.tsv"]))

    def row3(s, b, d):
        ra = tam.loc[(s, d, b)] if (s, d, b) in tam.index else None
        rb = tbm.loc[(s, d, b)] if (s, d, b) in tbm.index else None
        rc = tcp.loc[(s, d, b)] if (s, d, b) in tcp.index else None
        return (f"<tr><td>{esc(s)} · {'1 Mb' if b == 1000 else f'{b} kb'} · {d}</td>"
                f"<td class=num>{fmt(ra.fold) if ra is not None else '–'} ({fp(ra.p_perm) if ra is not None else '–'})</td>"
                f"<td class=num>{fp(ra.q_perm_global) if ra is not None else '–'}</td>"
                f"<td class=num>{(fmt(rb.irr) + ' (' + fp(rb.p) + ')') if rb is not None and np.isfinite(rb.irr) else '–'}</td>"
                f"<td class=num>{(fmt(rc.ratio) + ' [' + fmt(rc.ratio_lo) + '–' + fmt(rc.ratio_hi) + ']') if rc is not None else '–'}</td></tr>")
    cases = [(top.feature_set, int(top.bin_kb), top.direction), ("macheret2020:all", 500, "loss"), ("consensusA:core", 500, "loss"), ("union_all:all", 500, "loss")]
    P(f"""<h2 id="together">6 · The three tests side by side</h2>
<div class="tablebox"><table><thead><tr><th>set · bin · direction</th><th class=num>A: fold (p perm.)</th><th class=num>A: q perm. global</th><th class=num>B: adjusted IRR (p)</th><th class=num>C: pooled O/E [CI]</th></tr></thead>
<tbody>{"".join(row3(*c) for c in cases)}</tbody></table></div>""")
    P(method("<b>What is counted.</b> The strongest test A row and three reference sets at 500 kb (losses): the MiDAS union, the consensus core, and the union of all catalogues. "
             "<b>Reading.</b> A and B agree → the excess is not explained by RT, GC or mappability; A and C agree → it is not one donor or one clone. "
             "<b>Limitation.</b> The tests share the same events, so agreement is not independent replication."))

    pw = pd.read_csv(f"{res}/power_scaleup.tsv", sep="\t")
    pt = pw[(pw.feature_set == top.feature_set) & (pw.bin_kb == top.bin_kb) & (pw.direction == top.direction)]
    pl = pw[(pw.bin_kb == 500) & (pw.direction == "loss") & pw.feature_set.isin(["macheret2020:all", "consensusA:core"])]
    P(f"""<h2 id="limits">7 · What limits these results</h2><ul>
<li>Five donors and few events: the strongest row rests on {int(top.observed)} events; per-donor expectations are below 1–2.</li>
<li>No pre-declared primary test: the strongest row was picked from {n_glob} tests; only the global q values account for that.</li>
<li>Events cluster (clones, recurrent regions, donor-specific losses): the permutation p is larger than the exact p.</li>
<li>Resolution: breakend intervals (≥ 2 bins) are about as wide as the sites; the bins are not independent replications (different cells and widths), hence the cells500 and refined variants.</li>
<li>Catalogues and replication timing come from other cell types (tumour lines, lymphocytes, fibroblasts).</li>
<li>Scale-up to 49 donors (expected events × 9.05): {"; ".join(f"{r.feature_set} at {int(r.bin_kb)} kb {r.direction}: smallest fold detectable with 80% power at a global α ≈ {fmt(r.min_fold_80_global)}" for r in pd.concat([pt, pl]).drop_duplicates(["feature_set"]).itertuples())} (full document, section 10).</li></ul>
<p class="meta">Reproduce: <code>bash scripts/06_enrichment/run_all.sh</code>; this document: <code>python3 scripts/docs/build_enrichment_brief.py --res {R} --out docs/enrichment</code>.</p></main>""")

    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/enrichment_brief.html", "w") as f:
        f.write("\n".join(parts))
    print("wrote", f"{a.out}/enrichment_brief.html")


if __name__ == "__main__":
    main()

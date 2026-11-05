#!/usr/bin/env python3
"""Build docs/enrichment/tests_explainer.html: a teaching walk-through of enrichment tests A, B
and C, each worked step by step on real pilot data (what the test computes, what question it
answers, what it does not, and its limitations).

Worked examples (fixed): test A and B on interstitial E2-E4 losses at 500 kb against the HS68
MiDAS domains (macheret2020:HS68), the strongest test A result of the pilot (selected post hoc;
enrichment document, section 6.5); test C on the same subset plus the gains of the same set
(clonal weighting). The single event used to draw the null is chosen from the data (the
overlapping event with most null placements). Every number is read from the result tables or
recomputed with the functions of the test scripts themselves (test_a_permutation.py,
test_b_regression.py), so they match the published results; the null distribution of test A
is the exact Poisson-binomial of the real p_i.
Updated 2026-11-05 for the revised tests (grid-aligned null, exact and permutation p-values,
global q, test B window assignment by overlap, three models and conservative SE).

Inputs: results/enrichment/pilot/main/{events,breakends,usable_bins_rt}_<bin>kb.tsv,
test_a.tsv, test_a_shuffle_controls.tsv, test_b.tsv, test_c.tsv; data/features/cfs/
cfs_catalogs_hg19.tsv (+ cfs_derived_hg19.tsv); data/features/genes/; the 50 kb RT/GC/
mappability annotation of sWGS_relativeCN.

Usage: python3 scripts/docs/build_tests_explainer.py --res results/enrichment/pilot --out docs/enrichment
"""
import argparse
import html
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, f"{ROOT}/scripts/06_enrichment")
sys.path.insert(0, os.path.dirname(__file__))
from build_enrichment_html import CSS, flow, fmt, fp, table  # noqa: E402
from genome_annot import PERICEN  # noqa: E402
from test_a_permutation import Null, feature_sets, overlaps, pb_pmf  # noqa: E402
from test_b_regression import assign_windows, windows  # noqa: E402

SET, BIN, DIRN = "macheret2020:HS68", 500, "loss"
RT_FILE = f"{ROOT}/external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv"
CFS, GENES = f"{ROOT}/data/features/cfs/cfs_catalogs_hg19.tsv", f"{ROOT}/data/features/genes"
R = "results/enrichment/pilot"

EXTRA_CSS = """
:root{--obs:#c2571a;--hit:#0f6e73;--miss:#c3ccd4;--q:#56636f}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--obs:#f08a4b;--hit:#5cc3c7;--miss:#3a4651;--q:#9aa7b2}}
:root[data-theme="dark"]{--obs:#f08a4b;--hit:#5cc3c7;--miss:#3a4651;--q:#9aa7b2}
svg.chart{width:100%;height:auto;background:var(--surface);border:1px solid var(--rule);border-radius:6px}
svg.chart text{fill:var(--muted);font:11px var(--font)}svg.chart text.ink{fill:var(--fg)}
svg.chart .axis{stroke:var(--rule)}svg.chart .ref{stroke:var(--muted);stroke-dasharray:4 3}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px;margin:1rem 0}
.card{background:var(--surface);border:1px solid var(--rule);border-radius:8px;padding:.7rem .95rem;font-size:.9rem}
.card h4{margin:0 0 .3rem;font:600 .8rem/1.2 var(--mono);letter-spacing:.04em;text-transform:uppercase;color:var(--accent)}
.card p{margin:.2rem 0}.yes{border-left:3px solid var(--hit)}.no{border-left:3px solid var(--obs)}
.formula{font:.9rem/1.7 var(--mono);background:var(--code);border-radius:6px;padding:.6rem .9rem;margin:.6rem 0;max-width:90ch;overflow-x:auto}
ol.steps>li{margin:.45rem 0}.key{font-weight:600;color:var(--fg)}
.warn{border:1px solid var(--obs);border-radius:8px;padding:.6rem .9rem;margin:1rem 0;background:var(--surface);max-width:90ch}
.warn b{color:var(--obs)}
"""


def esc(x):
    return html.escape(str(x))


def mb(x, d=2):
    return f"{x / 1e6:.{d}f}"


def svg(w, h, body, label):
    return f'<svg class="chart" viewBox="0 0 {w} {h}" role="img" aria-label="{esc(label)}">{body}</svg>'


def caption(counted, how, limits):
    return (f'<div class="method"><b>What is counted.</b> {counted} <b>How.</b> {how} '
            f'<b>Limitations.</b> {limits}</div>')


def cards(answers, not_answers):
    a = "".join(f"<p>{x}</p>" for x in answers)
    n = "".join(f"<p>{x}</p>" for x in not_answers)
    return (f'<div class="cards"><div class="card yes"><h4>What it answers</h4>{a}</div>'
            f'<div class="card no"><h4>What it does not answer</h4>{n}</div></div>')


# ---------------------------------------------------------------- data


def test_a_example(res, fs):
    """Per-event o_i and p_i for the worked subset, plus the candidate placements of EXAMPLE."""
    ev = pd.read_csv(f"{res}/main/events_{BIN}kb.tsv", sep="\t", dtype={"chromosome": str})
    ev = ev[(ev.type == "interstitial") & ev.main_def & (ev.direction == DIRN) & ev.rt_quintile.notna()].copy()
    bins = pd.read_csv(f"{res}/main/usable_bins_rt_{BIN}kb.tsv", sep="\t", dtype={"chromosome": str})
    null = Null(bins)
    feat = fs[SET]
    rows, cand = [], {}
    for i, r in enumerate(ev.itertuples()):
        a, b = null.candidates(r.chromosome, int(r.rt_quintile), int(r.hi - r.lo))
        hit = overlaps(feat, r.chromosome, a, b)
        rows.append(dict(o=bool(overlaps(feat, r.chromosome, r.lo, r.hi)[0]), n_cand=len(a), n_hit=int(hit.sum()),
                         p=hit.mean() if len(a) else np.nan))
        cand[i] = (a, b, hit)
    ev = pd.concat([ev.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    i = ev[ev.o & (ev.n_hit > 0)].sort_values("n_cand").index[-1]
    sites = pd.read_csv(CFS, sep="\t")
    sites = sites[(sites.catalog == SET.split(":")[0]) & (sites.cell_type == SET.split(":")[1])]
    r = ev.loc[i]
    names = sites[(sites.chrom.str.replace("chr", "") == r.chromosome) & (sites.start <= r.hi) & (sites.end >= r.lo)]["name"]
    ex = dict(event=next(ev.loc[[i]].itertuples()), a=cand[i][0], b=cand[i][1], hit=cand[i][2],
              site=", ".join(sorted(set(names.astype(str)))))
    return ev, null, ex


def p_zero_hits(res, fs):
    """Events of SET that overlap the set although their null p_i is 0, per bin and direction."""
    out = []
    for b in (200, 500, 1000):
        ev = pd.read_csv(f"{res}/main/events_{b}kb.tsv", sep="\t", dtype={"chromosome": str})
        ev = ev[(ev.type == "interstitial") & ev.main_def & ev.rt_quintile.notna()]
        null = Null(pd.read_csv(f"{res}/main/usable_bins_rt_{b}kb.tsv", sep="\t", dtype={"chromosome": str}))
        for r in ev.itertuples():
            a, bb = null.candidates(r.chromosome, int(r.rt_quintile), int(r.hi - r.lo))
            o = overlaps(fs[SET], r.chromosome, r.lo, r.hi)[0]
            p = overlaps(fs[SET], r.chromosome, a, bb).mean() if len(a) else np.nan
            out.append(dict(bin_kb=b, direction=r.direction, o=o, p=p))
    return pd.DataFrame(out)


def test_b_example(res, fs):
    w = windows(f"{res}/main", RT_FILE, GENES)
    w["feature"] = [int(overlaps(fs[SET], c, s, e)[0]) for c, s, e in zip(w.chromosome, w.start, w.end)]
    ev = pd.read_csv(f"{res}/main/events_{BIN}kb.tsv", sep="\t", dtype={"chromosome": str})
    ev = ev[ev.type == "interstitial"]
    donors = sorted(ev.donor.unique())
    x = ev[ev.main_def & (ev.direction == DIRN)].copy()
    x["widx"] = assign_windows(x, w)
    x["in_window"] = x.widx >= 0
    cnt = x[x.in_window].groupby(["donor", "widx"]).size()
    df = pd.concat([w.assign(donor=d, widx=w.index) for d in donors], ignore_index=True)
    df["events"] = cnt.reindex(list(zip(df.donor, df.widx))).fillna(0).astype(int).values
    return w, df, x, donors


# ---------------------------------------------------------------- figures


def fig_concept():
    """Three cells of one clone -> three breakends (test C) -> one event (tests A, B)."""
    W, H, L, Rm = 900, 210, 150, 30
    X = lambda m: L + (m - 0) / 10 * (W - L - Rm)
    s = [f'<line class="axis" x1="{L}" y1="40" x2="{W - Rm}" y2="40" stroke-width="6" stroke-linecap="round"/>',
         f'<text x="{L - 10}" y="44" text-anchor="end" class="ink">chromosome</text>',
         f'<rect x="{X(4.2)}" y="28" width="{X(5.0) - X(4.2)}" height="24" rx="3" fill="var(--hit)" opacity=".85"><title>a fragile site of the set</title></rect>',
         f'<text x="{(X(4.2) + X(5.0)) / 2}" y="20" text-anchor="middle" class="ink">site</text>']
    cells = [("cell 1", 3.5, 5.0), ("cell 2", 3.5, 5.0), ("cell 3", 4.0, 5.5)]
    for i, (lab, lo, hi) in enumerate(cells):
        y = 78 + i * 22
        s.append(f'<text x="{L - 10}" y="{y + 4}" text-anchor="end">{lab} · breakend [lo, hi]</text>')
        s.append(f'<line x1="{X(lo)}" y1="{y}" x2="{X(hi)}" y2="{y}" stroke="var(--obs)" stroke-width="5" stroke-linecap="round"><title>{lab}: breakend interval {lo}–{hi} Mb</title></line>')
    y = 78 + 3 * 22 + 18
    s.append(f'<text x="{L - 10}" y="{y + 4}" text-anchor="end" class="ink">donor event (median)</text>')
    s.append(f'<line x1="{X(3.5)}" y1="{y}" x2="{X(5.0)}" y2="{y}" stroke="var(--obs)" stroke-width="9" stroke-linecap="round"><title>collapsed event: median lo / hi of its breakends</title></line>')
    s.append(f'<text x="{X(6.0)}" y="{88}" class="ink">test C: 3 breakends, each counts</text>')
    s.append(f'<text x="{X(6.0)}" y="{y + 4}" class="ink">tests A and B: 1 event</text>')
    for m in range(0, 11, 2):
        s.append(f'<text x="{X(m)}" y="{H - 6}" text-anchor="middle">{m} Mb</text>')
    return svg(W, H, "".join(s), "Breakends of three cells of one clone collapse into one donor event")


def fig_null(ex, null, fs):
    r, a, b, hit = ex["event"], ex["a"], ex["b"], ex["hit"]
    c = r.chromosome
    bins = null.bins[null.bins.chromosome == c]
    clen = float(bins.end.max())
    W, L, Rm = 960, 190, 20
    X = lambda v: L + v / clen * (W - L - Rm)
    q = int(r.rt_quintile)
    allowed = set(a.tolist())
    s = []
    rows = [("usable bins (allowed start?)", 26), (f"{SET} sites", 62), ("null placements (same RT quintile)", 98),
            ("observed event", 168)]
    for lab, y in rows:
        s.append(f'<text x="{L - 10}" y="{y + 4}" text-anchor="end" class="ink">{esc(lab)}</text>')
    for t in bins.itertuples():
        ok = int(t.start) in allowed
        f = "var(--q)" if ok else "var(--miss)"
        s.append(f'<rect x="{X(t.start):.1f}" y="18" width="{max(X(t.end) - X(t.start), 1):.1f}" height="16" fill="{f}">'
                 f'<title>chr{c}:{mb(t.start)}–{mb(t.end)} Mb · bin RT {t.rt:.1f} · '
                 f'{"a placement starting here has the event width and RT quintile Q" + str(q) if ok else "not an allowed start"}</title></rect>')
    lo, hi = PERICEN[int(c)]
    s.append(f'<rect x="{X((lo - 1) * 1e6):.1f}" y="12" width="{X((hi + 1) * 1e6) - X((lo - 1) * 1e6):.1f}" height="170" '
             f'fill="none" stroke="var(--muted)" stroke-dasharray="3 3"><title>pericentromeric block ± 1 Mb: no placement may overlap it</title></rect>')
    s.append(f'<text x="{X((lo + hi) / 2 * 1e6):.1f}" y="9" text-anchor="middle">pericentromere ± 1 Mb (excluded)</text>')
    S, E = fs[SET].get(c, (np.array([]), np.array([])))
    for s0, e0 in zip(S, E):
        s.append(f'<rect x="{X(s0):.1f}" y="52" width="{max(X(e0) - X(s0), 2):.1f}" height="20" rx="2" fill="var(--fg)">'
                 f'<title>{SET} site chr{c}:{mb(s0)}–{mb(e0)} Mb</title></rect>')
    for i, (x0, x1, h) in enumerate(zip(a, b, hit)):
        y = 88 + (i % 4) * 11
        col = "var(--hit)" if h else "var(--miss)"
        s.append(f'<line x1="{X(x0):.1f}" y1="{y}" x2="{X(x1):.1f}" y2="{y}" stroke="{col}" stroke-width="{7 if h else 6}" '
                 f'stroke-linecap="round"><title>candidate {mb(x0)}–{mb(x1)} Mb · {"overlaps a site" if h else "no overlap"}</title></line>')
    s.append(f'<line x1="{X(r.lo):.1f}" y1="168" x2="{X(r.hi):.1f}" y2="168" stroke="var(--obs)" stroke-width="9" stroke-linecap="round">'
             f'<title>observed: {r.donor} chr{c}:{mb(r.lo)}–{mb(r.hi)} Mb loss, RT quintile Q{q}</title></line>')
    s.append(f'<text x="{X(r.lo):.1f}" y="192" class="ink">{r.donor} loss, chr{c}:{mb(r.lo, 1)}–{mb(r.hi, 1)} Mb (Q{q}), p_i = {hit.mean():.3f}</text>')
    for m in range(0, int(clen / 1e6) + 1, 20):
        s.append(f'<text x="{X(m * 1e6):.1f}" y="214" text-anchor="middle">{m} Mb</text>')
    leg = [("var(--q)", f"allowed start (placement RT in Q{q})"), ("var(--miss)", "not allowed / miss"),
           ("var(--hit)", "placement that overlaps a site"), ("var(--obs)", "observed event")]
    for i, (col, lab) in enumerate(leg):
        x = L + i * 190
        s.append(f'<rect x="{x}" y="226" width="12" height="12" rx="2" fill="{col}"/><text x="{x + 18}" y="236">{esc(lab)}</text>')
    return svg(W, 248, "".join(s), f"Null placements of one event on chromosome {c}")


def fig_hist(p, O, E):
    pmf = pb_pmf(p)
    kmax = max(int(np.nonzero(pmf > 1e-4)[0].max()), O) + 1
    cnt = pmf[:kmax + 1]
    W, H, L, B, T = 900, 260, 60, 40, 20
    bw = (W - L - 20) / (kmax + 1)
    top = cnt.max()
    Y = lambda f: H - B - f / top * (H - B - T)
    s = [f'<line class="axis" x1="{L}" y1="{H - B}" x2="{W - 20}" y2="{H - B}"/>']
    for k, f in enumerate(cnt):
        col = "var(--obs)" if k >= O else "var(--miss)"
        x = L + k * bw + 2
        h = H - B - Y(f)
        if f > 0:
            s.append(f'<path d="M{x:.1f},{H - B} v{-max(h - 4, 0):.1f} q0,-4 4,-4 h{bw - 12:.1f} q4,0 4,4 v{max(h - 4, 0):.1f} z" fill="{col}">'
                     f'<title>{k} overlaps: probability {f:.2e}</title></path>')
        s.append(f'<text x="{x + (bw - 4) / 2:.1f}" y="{H - B + 15}" text-anchor="middle">{k}</text>')
    xe = L + E * bw + bw / 2
    s.append(f'<line class="ref" x1="{xe:.1f}" y1="{T}" x2="{xe:.1f}" y2="{H - B}"/>'
             f'<text x="{xe + 5:.1f}" y="{T + 10}" class="ink">E = Σp<tspan baseline-shift="sub" font-size="8">i</tspan> = {E:.2f}</text>')
    xo = L + O * bw + bw / 2
    tail = float(pmf[O:].sum())
    s.append(f'<text x="{xo:.1f}" y="{H - B - 30}" text-anchor="middle" class="ink">O = {O}</text>'
             f'<text x="{xo:.1f}" y="{H - B - 16}" text-anchor="middle">P(S ≥ O) = {tail:.1e}</text>')
    s.append(f'<text x="{(W + L) / 2}" y="{H - 6}" text-anchor="middle">overlapping events under the null (S)</text>')
    s.append(f'<text x="14" y="{(H - B) / 2}" transform="rotate(-90 14 {(H - B) / 2})" text-anchor="middle">probability</text>')
    return svg(W, H, "".join(s), "Exact null distribution of the number of overlapping events"), tail


def fig_donors(tc):
    W, L, Rm, rowh = 900, 120, 220, 30
    t = tc.copy()
    t["order"] = t.donor.map(lambda d: 99 if d == "pooled" else 0)
    t = t.sort_values(["order", "donor"])
    xmax = max(np.nanmax(t.ratio_hi), np.nanmax(t.ratio)) * 1.08
    X = lambda v: L + v / xmax * (W - L - Rm)
    H = 30 + rowh * len(t) + 30
    s = [f'<line class="ref" x1="{X(1):.1f}" y1="14" x2="{X(1):.1f}" y2="{H - 30}"/>',
         f'<text x="{X(1):.1f}" y="10" text-anchor="middle">O/E = 1 (no enrichment)</text>']
    for i, r in enumerate(t.itertuples()):
        y = 34 + i * rowh
        pooled = r.donor == "pooled"
        col = "var(--fg)" if pooled else "var(--obs)"
        s.append(f'<text x="{L - 10}" y="{y + 4}" text-anchor="end" class="ink">{"pooled" if pooled else r.donor}</text>')
        if np.isfinite(r.ratio_lo) and np.isfinite(r.ratio_hi):
            s.append(f'<line x1="{X(r.ratio_lo):.1f}" y1="{y}" x2="{X(r.ratio_hi):.1f}" y2="{y}" stroke="{col}" stroke-width="2"/>')
        if np.isfinite(r.ratio):
            s.append(f'<circle cx="{X(r.ratio):.1f}" cy="{y}" r="{6 if pooled else 5}" fill="{col}" stroke="var(--surface)" stroke-width="2">'
                     f'<title>{r.donor}: O/E {r.ratio:.2f} [{r.ratio_lo:.2f}–{r.ratio_hi:.2f}], {r.observed} observed / {r.expected:.2f} expected, {r.breakends} breakends in {r.cells} cells</title></circle>')
        s.append(f'<text x="{W - Rm + 12}" y="{y + 4}">{r.observed} obs / {r.expected:.2f} exp · {r.breakends} breakends</text>')
    for v in range(0, int(xmax) + 1, 2):
        s.append(f'<text x="{X(v):.1f}" y="{H - 10}" text-anchor="middle">{v}</text>')
    return svg(W, H, "".join(s), "Observed / expected per donor with cell-bootstrap 95% CI")


# ---------------------------------------------------------------- sections


def section_shared():
    return f"""<h2 id="shared">1 · The shared question and its vocabulary</h2>
<p>All three tests ask the same question: <span class="key">do copy-number breakends of normal breast cells fall inside a set of
fragile sites more often than they would by chance?</span> They differ in what they count, what “by chance” means and which
confounders they neutralise. Running three tests is a way of checking that an answer does not depend on one arbitrary choice.</p>
<ul>
<li><b>Breakend</b>: the boundary between two segments of different copy-number state in <i>one cell</i>. Its position is uncertain
by about one bin, so it is stored as an interval <code>[lo, hi]</code> (from the start of the last bin of the left segment to the end
of the first bin of the right one).</li>
<li><b>Event</b>: breakends of the <i>same donor</i>, chromosome, direction and type that lie within 2 Mb of each other, collapsed
into one. A clone of 100 cells that inherited the same deletion is 100 breakends but 1 event: the break happened once.</li>
<li><b>Set</b>: a list of genomic intervals (a CFS catalogue, a consensus, or a control such as genes &gt; 1 Mb). “Overlap” means
the breakend or event interval shares at least one base with any site of the set.</li>
<li><b>O, E and fold</b>: O is how many events overlap the set; E is how many would overlap if the events had been placed at
random under a defined null; fold = O / E. A fold of 1 means “as expected”, 2 means “twice as often as expected”.</li>
<li><b>Subsets</b>: every test is run per bin size (200 kb, 500 kb, 1 Mb), direction (loss, gain, all) and deviation stratum
(E1–E4; “main” = E2–E4, the paper's 1.3 / 0.7 ratio threshold). Only interstitial events (not touching a centromere) enter.</li>
</ul>
<figure>{fig_concept()}{caption(
        "A schematic (not data): three cells of one donor carrying the same loss, and the site set.",
        "Breakends of the same donor, chromosome, direction and type within 2 Mb are linked; the event takes the medians of their lo and hi.",
        "Linking at 2 Mb can merge two independent nearby breaks, or split one event whose cells have different intervals.")}</figure>
<p><b>Why a null at all.</b> Simply counting overlaps says nothing: a set that covers 10% of the genome will catch about 10% of
random events. Worse, fragile sites are not random places: they are late-replicating, AT-rich, often inside very large
genes, and they are unevenly spread across chromosomes. Copy-number breaks (and copy-number calling artefacts) follow some of the
same features. A useful null must therefore place random events in genomic context comparable to the real ones. The three tests
do this in different ways.</p>"""


def section_a(ev, ex, null, fs, ta, sh, pz):
    r = ex["event"]
    n_c, n_h = len(ex["a"]), int(ex["hit"].sum())
    O, E = int(ev.o.sum()), float(ev.p.sum())
    row = ta[(ta.feature_set == SET) & (ta.bin_kb == BIN) & (ta.direction == DIRN) & (ta.stratum == "main")].iloc[0]
    hist, p_ex = fig_hist(ev.p.values, O, E)
    hits = ev[ev.o].sort_values(["donor", "chromosome", "pos"])
    tab = pd.DataFrame({
        "donor": hits.donor, "event interval": [f"chr{c}:{mb(lo)}–{mb(hi)} Mb" for c, lo, hi in zip(hits.chromosome, hits.lo, hits.hi)],
        "RT quintile": hits.rt_quintile.astype(int).map(lambda q: f"Q{q}"), "cells": hits.n_cells,
        "null placements": hits.n_cand, "placements overlapping": hits.n_hit, "p_i": hits.p.map(lambda v: f"{v:.3f}")})
    s = sh[(sh.feature_set == SET) & (sh.bin_kb == BIN) & (sh.direction == DIRN)].fold
    pzs = pz.groupby(["bin_kb", "direction"]).apply(lambda g: pd.Series({"O": int(g.o.sum()), "of which p_i = 0": int((g.o & (g.p == 0)).sum())})).reset_index()
    pzs = pzs[pzs.direction != "mixed"]
    pz_tab = table(pzs.assign(bin_kb=pzs.bin_kb.map(lambda b: f"{b} kb")), ["bin_kb", "direction", "O", "of which p_i = 0"], num=("O", "of which p_i = 0"))
    return f"""<h2 id="testa">2 · Test A: positional permutation of events</h2>
<p class="meta">Script <code>scripts/06_enrichment/test_a_permutation.py</code> · unit: one collapsed event · worked example:
interstitial E2–E4 losses at {BIN} kb against <code>{SET}</code> (MiDAS domains of HS68 normal fibroblasts, Macheret 2020).</p>
<p><b>The idea in one sentence.</b> Take each real event and ask: <i>if this exact event (same chromosome, same width, same
replication timing) had happened somewhere else, how likely would it have been to hit a fragile site?</i> Sum those
probabilities over all events and compare with how many actually hit one.</p>

<h3>Step 1 · One event, its null placements and its p<sub>i</sub></h3>
<ol class="steps">
<li>The event: donor {r.donor}, a loss with interval <code>chr{r.chromosome}:{mb(r.lo)}–{mb(r.hi)} Mb</code> (width {mb(r.hi - r.lo + 1)} Mb),
overlapping the HS68 site(s) annotated {esc(ex["site"]) or "–"} in Macheret et al. 2020. Its replication timing (mean RT of the usable bins
of the interval) falls in <b>quintile Q{int(r.rt_quintile)}</b> of the usable {BIN} kb bins (Q1 = latest replicating fifth of the genome).</li>
<li>Allowed alternative positions: the event, with its width kept, is started at every usable {BIN} kb bin of the same chromosome
(the observed intervals also start at a bin start, so placements sit on the same grid). For each placement the RT is computed
exactly as for the event (mean RT of the usable bins it covers), and only placements in Q{int(r.rt_quintile)} are kept; placements touching
the pericentromeric block ± 1 Mb are removed. That leaves <b>{n_c} placements</b>.</li>
<li>Of these, <b>{n_h}</b> overlap an HS68 site. All placements are equally likely, so the event's null probability of overlap is
exact, no simulation needed: <b>p<sub>i</sub> = {n_h} / {n_c} = {n_h / n_c:.3f}</b>.</li>
<li>Observed: the real interval does overlap a site, so <b>o<sub>i</sub> = 1</b>.</li>
</ol>
<figure>{fig_null(ex, null, fs)}{caption(
        f"Chromosome {r.chromosome}: usable {BIN} kb bins coloured by whether their RT quintile equals the event's (Q{int(r.rt_quintile)}); the HS68 MiDAS sites; the {n_c} null placements of the event (staggered in four rows only for legibility), coloured when they overlap a site; and the observed event. Hover any mark for its coordinates.",
        f"Placements = intervals of the event's width starting at each usable bin of chromosome {r.chromosome}, kept when their interval RT is in quintile Q{int(r.rt_quintile)}, excluding the pericentromere ± 1 Mb (test_a_permutation.Null.candidates).",
        "The null only moves an event along its own chromosome and within its RT quintile; anything else that makes a region break-prone (gene size, transcription, GC) is not held constant.")}</figure>

<h3>Step 2 · All events: observed O and expected E</h3>
<p>The same is done for each of the <b>{len(ev)}</b> interstitial E2–E4 losses at {BIN} kb. Most have p<sub>i</sub> close to 0
(median {np.median(ev.p):.3f}) because HS68 sites cover little of the genome. Then:</p>
<div class="formula">O = Σ o<sub>i</sub> = {O} events overlap a site<br>E = Σ p<sub>i</sub> = {E:.2f} events expected to overlap by chance<br>
fold = O / E = {O} / {E:.2f} = <b>{O / E:.2f}</b></div>
<p>These are the {O} overlapping events (the other {len(ev) - O} contribute only their p<sub>i</sub> to E):</p>
{table(tab, list(tab.columns), num=("cells", "null placements", "placements overlapping", "p_i"))}
<div class="method"><b>What is counted.</b> Main-set interstitial E2–E4 loss events at {BIN} kb whose interval overlaps an HS68 MiDAS domain; “cells” = cells carrying the event (it still counts once). <b>How.</b> As in step 1, for every event. <b>Limitations.</b> Rows with p<sub>i</sub> = 0 are discussed below (“a null that cannot reproduce the observation”).</div>

<h3>Step 3 · Is O = {O} surprising? The p-value</h3>
<p>Under the null each event overlaps independently with its own probability p<sub>i</sub>, so the number of overlaps S is a sum of
{len(ev)} unequal coin flips (a Poisson-binomial distribution). Its distribution can be computed exactly, by convolving the events
one at a time, so no simulation is needed: <b>p exact = P(S ≥ {O})</b>.</p>
<figure>{hist}{caption(
        f"Exact distribution of S = Σ Bernoulli(p<sub>i</sub>) for the {len(ev)} real p<sub>i</sub>; bars at or above O = {O} (orange) form the p-value tail.",
        f"Recomputed here: p = {p_ex:.1e}; in <code>test_a.tsv</code> p exact = {fp(row.p_enrich)}.",
        "Assumes events are independent of each other; see the permutation p below for a null that keeps their clustering.")}</figure>
<p><b>Confidence interval of the fold.</b> The {len(ev)} pairs (o<sub>i</sub>, p<sub>i</sub>) are resampled with replacement 2,000
times and O/E recomputed: fold {row.fold:.2f} [{row.fold_lo:.2f}–{row.fold_hi:.2f}]. It reflects how much the fold depends on
which events happened to be sampled, not uncertainty in the site set.</p>
<p><b>The permutation p (the reference p-value).</b> Events are not independent: a donor's recurrent losses, nearby events and the
clustering of sites along the genome all make the exact p too optimistic. So the HS68 set itself is shuffled {len(s):,} times (every
site moved to a random usable bin of the same chromosome and RT quintile, width kept) and test A is re-run on each shuffled set. The
shuffled folds have median {s.median():.2f} and 2.5–97.5% range {s.quantile(.025):.2f}–{s.quantile(.975):.2f}; the real set gives
{row.fold:.2f}. <b>p perm.</b> = (1 + shuffles with fold ≥ {row.fold:.2f}) / ({len(s):,} + 1) = <b>{fp(row.p_perm)}</b>, against p exact {fp(row.p_enrich)}.
The shuffles centring near 1 also checks that the null is not biased.</p>
<p><b>Multiple testing.</b> This is one of about 200 main-stratum tests (sets × bins × directions). Benjamini–Hochberg gives
q family = {fp(row.q_family)} (only within the HS68 set), q global = {fp(row.q_global)} (all sets) and q perm. global = {fp(row.q_perm_global)}
(all sets, on permutation p). With no set chosen in advance, the global values are the ones that control the error.</p>
{cards([
        "Do events (each break counted once, however many cells inherited it) overlap the set more than events of the same width and replication timing placed elsewhere on the same chromosome?",
        "How large is the excess (fold) and how consistent across bins, directions and strata?"],
        ["Whether the excess is explained by large genes, GC or mappability (only RT and chromosome are matched): that is test B.",
         "Whether one donor or one expanded clone drives it: that is test C (and the leave-one-donor-out runs).",
         f"Which individual event is due to fragility: a fold of {row.fold:.1f} with {O} overlaps still means some of the {O} would be there by chance (E ≈ {E:.1f})."])}
<h3>Limitations of test A</h3>
<ul>
<li><b>RT is the only genomic covariate.</b> Fragile sites concentrate in very large, late-replicating, transcribed genes. Matching
RT removes part of that, not gene size or transcription.</li>
<li><b>Independence.</b> The exact p and the event bootstrap assume each event's overlap is independent. Events of the same donor or
region may not be; clonal collapse removes the most obvious dependence, the permutation p keeps the clustering of events, and a
bootstrap over donors (in <code>test_a.tsv</code>) shows the spread between donors.</li>
<li><b>Coarse positions.</b> Event intervals are 0.5 to several Mb wide, comparable to the sites. A wide interval overlaps more easily;
the null keeps the width so this is fair, but it limits how precisely a break can be assigned to a site.</li>
<li><b>Many sets, many subsets.</b> q family corrects within one set only; q global and q perm. global correct across all sets.
With no pre-registered hypothesis, one isolated small p is a lead, not a finding.</li>
</ul>
<div class="warn"><b>A null that could not reproduce the observation (fixed 5 November 2026).</b> In the first version, null
placements were centred on single bins that had to be in the event's RT quintile, while the event's quintile came from the mean RT
of its whole interval. An event spanning a Q5 and a Q3 bin averaged to Q4, none of its own bins was an allowed centre, and it got
p<sub>i</sub> = 0 although it overlapped a site: +1 to O and 0 to E. The null now places intervals on the same bin grid and computes
their RT exactly as for events, so every event can land where it was observed. Overlapping events with p<sub>i</sub> = 0 for this set
with the current null:
{pz_tab}</div>
{flow([f"{R}/main/events_{BIN}kb.tsv", f"{R}/main/usable_bins_rt_{BIN}kb.tsv", "data/features/cfs/cfs_catalogs_hg19.tsv", f"{R}/main/test_a.tsv", f"{R}/main/test_a_shuffle_controls.tsv"],
           ["scripts/docs/build_tests_explainer.py (calls test_a_permutation.Null, overlaps, feature_sets)"], ["this section (not saved separately)"],
           f"o<sub>i</sub>, p<sub>i</sub> recomputed with the test's own functions for the {SET}, {BIN} kb, loss, E2–E4 subset; fold, CI, p and q read from test_a.tsv; shuffle quantiles from test_a_shuffle_controls.tsv.")}"""


def section_b(w, df, x, donors, tb, ta, ex):
    f1, f0 = w[w.feature == 1], w[w.feature == 0]
    e1, e0 = int(df[df.feature == 1].events.sum()), int(df[df.feature == 0].events.sum())
    n1, n0 = len(f1) * len(donors), len(f0) * len(donors)
    r1, r0 = e1 / n1, e0 / n0
    lost = int((~x.in_window).sum())
    c0, p0 = ex["event"].chromosome, ex["event"].pos
    ex_w = df[(df.chromosome == c0) & (df.start > p0 - 3e6) & (df.start < p0 + 2e6)]
    piv = ex_w.pivot_table(index=["chromosome", "start", "end", "feature", "rt", "gc", "map", "genes300", "usable"], columns="donor", values="events").reset_index()
    tab = pd.DataFrame({"window": [f"chr{c}:{mb(s - 1, 0)}–{mb(e, 0)} Mb" for c, s, e in zip(piv.chromosome, piv.start, piv.end)],
                        "overlaps HS68": piv.feature.map({1: "yes", 0: "no"}), "RT (z)": piv.rt.round(2), "GC (z)": piv.gc.round(2),
                        "mappability (z)": piv["map"].round(2), "fraction in genes > 300 kb": piv.genes300.round(2),
                        "usable fraction": piv.usable.round(2), **{d: piv[d].astype(int) for d in donors}})
    cov = w.groupby("feature")[["rt", "gc", "map", "genes300"]].mean().round(2)
    m = tb[(tb.feature_set == SET) & (tb.bin_kb == BIN) & (tb.direction == DIRN) & (tb.stratum == "main")].set_index("model")
    u, a, g = m.loc["crude"], m.loc["adjusted"], m.loc["adjusted_genes300"]
    ra = ta[(ta.feature_set == SET) & (ta.bin_kb == BIN) & (ta.direction == DIRN) & (ta.stratum == "main")].iloc[0]
    fitted = tb[tb.irr.notna()]
    return f"""<h2 id="testb">3 · Test B: regression on genomic windows</h2>
<p class="meta">Script <code>scripts/06_enrichment/test_b_regression.py</code> · unit: one 1 Mb window in one donor · same
worked example (E2–E4 losses at {BIN} kb, <code>{SET}</code>).</p>
<p><b>The idea in one sentence.</b> Cut the genome into windows, count events in each, and ask whether windows that contain a
fragile site have a higher event rate <i>after accounting for</i> the other things that make a window break-prone (RT, GC,
mappability, chromosome, donor, and optionally large genes). It is the regression counterpart of test A: instead of matching on RT
it adjusts for several covariates at once.</p>

<h3>Step 1 · Build the table: one row per donor × window</h3>
<p>The genome is divided into the usable 1 Mb windows outside pericentromeric blocks ± 1 Mb: <b>{len(w):,} windows</b>, the same for
the three event bin sizes. With {len(donors)} donors that is {len(w) * len(donors):,} rows. Each event is assigned to the eligible window
with which its interval [lo, hi] overlaps most. {lost} of the {len(x)} events overlap no eligible window (non-usable or too close to a
centromere) and drop out, leaving {len(x) - lost}. For each window the script records:</p>
{table(tab, list(tab.columns), num=tuple(c for c in tab.columns if c not in ("window", "overlaps HS68")))}
<div class="method"><b>What is counted.</b> Windows of chromosome {c0} around the event drawn in test A ({ex["event"].donor}, chr{c0}:{mb(ex["event"].lo, 1)}–{mb(ex["event"].hi, 1)} Mb); the donor columns are the number of E2–E4 loss events assigned to the window. <b>How.</b> RT, GC and mappability are the means of the 50 kb annotation of sWGS_relativeCN, z-scored over all windows; “fraction in genes &gt; 300 kb” from refGene loci; usable fraction = share of the window's 50 kb bins that pass the QDNAseq filters (used as an offset). <b>Limitations.</b> Overlap with the set is yes / no, whatever fraction of the window the site covers.</div>

<h3>Step 2 · The crude comparison (what the crude model estimates)</h3>
{table(pd.DataFrame({"windows": ["overlap HS68", "do not overlap"], "windows × donors": [f"{n1:,}", f"{n0:,}"], "events": [e1, e0],
                     "events per window and donor": [f"{r1:.4f}", f"{r0:.4f}"]}), ["windows", "windows × donors", "events", "events per window and donor"], num=("windows × donors", "events", "events per window and donor"))}
<div class="formula">rate ratio = {r1:.4f} / {r0:.4f} = <b>{r1 / r0:.2f}</b></div>
<p>This is the incidence rate ratio (IRR) of the set. The crude model (<code>events ~ feature + donor</code>, with the
usable fraction as offset) gives practically the same number, IRR {u.irr:.2f} [{u.irr_lo:.2f}–{u.irr_hi:.2f}], because donor only
shifts the baseline rate and almost every window is fully usable.</p>
<p>Note that <b>{e1}</b> events fall in HS68 windows here, against <b>{int(ra.observed)}</b> overlapping events in test A. The two tests
define “hitting the set” differently: test A asks whether the event interval touches a site, test B whether the window the event is
assigned to overlaps a site.</p>

<h3>Step 3 · Adjusting for what else makes a window break-prone</h3>
<p>HS68 windows are not average windows. Their mean covariates (z-scores; genes as a fraction):</p>
{table(cov.reset_index().assign(feature=lambda d: d.feature.map({1: "overlap HS68", 0: "do not overlap"})).rename(columns={"feature": "windows", "rt": "RT", "gc": "GC", "map": "mappability", "genes300": "fraction in genes > 300 kb"}),
           ["windows", "RT", "GC", "mappability", "fraction in genes > 300 kb"], num=("RT", "GC", "mappability", "fraction in genes > 300 kb"))}
<p>Two adjusted models are fitted:</p>
<div class="formula">adjusted: log E[events] = β·feature + γ<sub>1</sub>·RT + γ<sub>2</sub>·GC + γ<sub>3</sub>·mappability + chromosome + donor + log(usable fraction)<br>
adjusted_genes300: the same + γ<sub>4</sub>·genes300<br>IRR = exp(β)</div>
<p>Result: adjusted IRR <b>{a.irr:.2f} [{a.irr_lo:.2f}–{a.irr_hi:.2f}]</b>, p = {fp(a.p)}, q family = {fp(a.q_family)}, q global = {fp(a.q_global)};
adjusted_genes300 IRR <b>{g.irr:.2f} [{g.irr_lo:.2f}–{g.irr_hi:.2f}]</b>, p = {fp(g.p)} (IRR of genes300 {g.irr_genes300:.2f}, of RT {g.irr_rt:.2f} per SD).
Large genes are probably part of how fragile sites break (a mediator, not only a confounder), so adjusting for them is conservative:
it can remove part of the very effect being tested.</p>
<p><b>Which distribution and which errors.</b> Counts are modelled as negative binomial if they are more variable than Poisson. Here
the estimated over-dispersion is negative (α = {a.alpha:.2f}): almost every donor × window has 0 or 1 events, so the model is Poisson
(<code>family = {a.family}</code>, as in {int((fitted.family == 'poisson').sum())} of the {len(fitted)} fitted models). Windows of one
chromosome are not independent, so the standard error of the set is computed both from the model and cluster-robust by chromosome,
and the larger of the two is used ({a.se_source} here). When no event falls in windows of a set the model has no finite estimate
(separation) and the row is left empty instead of reporting p = 0.</p>
{cards([
        "Is the event rate higher in windows that contain the set, after adjusting for RT, GC, mappability, chromosome and donor?",
        "How much of it remains when large genes are adjusted for too? (Compare adjusted and adjusted_genes300.)"],
        ["Anything about the exact position of a break inside the window: an event counts in one window only.",
         "Whether one clone drives the result (events are collapsed, as in test A)."])}
<h3>Limitations of test B</h3>
<ul>
<li><b>Fixed 1 Mb grid.</b> A site near a window border may be credited to the neighbouring window; {lost} of {len(x)} events in this
example overlap no eligible window and are excluded.</li>
<li><b>Binary overlap.</b> A window touched by 20 kb of a site counts the same as one fully covered.</li>
<li><b>Sparse counts.</b> With few events and many fixed effects (22 chromosomes, 5 donors), estimates are unstable for small sets;
CIs are wide.</li>
<li><b>Linear covariates.</b> RT, GC and mappability enter linearly on the log scale; a non-linear effect of RT is only partly
removed (test A's quintile matching is more flexible for RT).</li>
</ul>
{flow([f"{R}/main/events_{BIN}kb.tsv", f"{R}/main/usable_bins_rt_1000kb.tsv", "bins_with_replicationTime_50kb.tsv", "data/features/genes/large_genes_300kb.tsv", f"{R}/main/test_b.tsv"],
           ["scripts/docs/build_tests_explainer.py (calls test_b_regression.windows and assign_windows)"], ["this section (not saved separately)"],
           "window table and crude rates recomputed with the test's own functions; IRRs, p, q, α, family and SE source read from test_b.tsv.")}"""


def section_c(res, tc, ta):
    sub = tc[(tc.feature_set == SET) & (tc.bin_kb == BIN) & (tc.direction == DIRN)]
    pooled = sub[sub.donor == "pooled"].iloc[0]
    don = sub[sub.donor != "pooled"]
    Rr = pooled.observed / pooled.expected
    use = don[don.expected >= 1]
    terms = [(r.donor, r.observed, r.expected, (r.observed - Rr * r.expected) ** 2 / (Rr * r.expected)) for r in use.itertuples()]
    X2 = sum(t[3] for t in terms)
    dfree = len(terms) - 1
    x2_lines = "<br>".join(f"{d}: ({o} − {Rr:.2f} × {e:.2f})² / ({Rr:.2f} × {e:.2f}) = {t:.2f}" for d, o, e, t in terms)
    be = pd.read_csv(f"{res}/main/breakends_{BIN}kb.tsv", sep="\t", dtype={"chromosome": str})
    be = be[(be.type == "interstitial") & be.stratum.isin(["E2", "E3", "E4"])]
    per = be.groupby(["donor", "direction"]).agg(breakends=("cell_id", "size"), cells=("cell_id", "nunique"), events=("event_id", "nunique")).reset_index()
    per = per[per.direction.isin(["loss", "gain"])]
    per["breakends per event"] = (per.breakends / per.events).round(1)
    big = be.groupby(["donor", "event_id", "direction", "chromosome"]).agg(cells=("cell_id", "nunique"), lo=("lo", "min"), hi=("hi", "max")).reset_index().sort_values("cells").iloc[-1]
    g = tc[(tc.feature_set == SET) & (tc.bin_kb == BIN) & (tc.direction == "gain")].set_index("donor")
    ga = ta[(ta.feature_set == SET) & (ta.bin_kb == BIN) & (ta.direction == "gain") & (ta.stratum == "main")].iloc[0]
    al = tc[(tc.feature_set == SET) & (tc.bin_kb == BIN) & (tc.direction == "all") & (tc.donor == "pooled")].iloc[0]
    return f"""<h2 id="testc">4 · Test C: per breakend and per donor</h2>
<p class="meta">Script <code>scripts/06_enrichment/test_c_per_donor.py</code> · unit: one breakend of one cell (no collapsing) ·
same worked example, then the gains of the same set.</p>
<p><b>The idea in one sentence.</b> Repeat test A's observed / expected comparison, but count every cell's breakend (so an expanded
clone weighs as many times as it has cells) and compute it separately for each donor, to see whether the excess is shared or
driven by one donor.</p>

<h3>Step 1 · Same p<sub>i</sub>, different unit</h3>
<p>Each breakend gets o (does its interval overlap the set?) and p (its exact null probability, computed exactly as in test A from its
own chromosome, RT quintile and width). Then per donor O<sub>d</sub> = Σo, E<sub>d</sub> = Σp, ratio O<sub>d</sub>/E<sub>d</sub>.
The pooled ratio is ΣO / ΣE over all donors.</p>
<p><b>Why this matters: clones.</b> The breast epithelium contains expanded clones. Breakends per event by donor and direction
({BIN} kb, E2–E4):</p>
{table(per, ["donor", "direction", "breakends", "cells", "events", "breakends per event"], num=("breakends", "cells", "events", "breakends per event"))}
<p>The largest single event is a {big.direction} in {big.donor}, chr{big.chromosome}:{mb(big.lo, 1)}–{mb(big.hi, 1)} Mb, present in
<b>{big.cells} cells</b>. In tests A and B it counts once; in test C it counts {big.cells} times. Test C asks “are the
copy-number changes we <i>observe in cells</i> enriched?” (which mixes break rate with clonal expansion), while A asks “are the
<i>break events</i> enriched?”.</p>

<h3>Step 2 · Per-donor ratios and their uncertainty</h3>
<figure>{fig_donors(sub)}{caption(
        f"Interstitial E2–E4 loss breakends at {BIN} kb against {SET}: observed overlaps, expected overlaps (Σp) and O/E per donor and pooled.",
        "95% CI by bootstrap of cells (2,000 resamples): the cells of a donor are resampled with replacement, carrying all their breakends; the pooled CI resamples cells within each donor.",
        "Donors with E &lt; 1 have very wide CIs (when O = 0 the upper limit is the one-sided exact Poisson bound 3.0 / E); cells of one clone are not independent, so the cell bootstrap understates the uncertainty when clones dominate.")}</figure>
<p>Pooled: O/E = {pooled.observed} / {pooled.expected:.2f} = <b>{pooled.ratio:.2f} [{pooled.ratio_lo:.2f}–{pooled.ratio_hi:.2f}]</b>.
{int((don.ratio > 1).sum())} of the {len(don)} donors point above 1 ({", ".join(f"{r.donor} {r.ratio:.2f}" for r in don.itertuples())}). With {int(don.observed.min())}–{int(don.observed.max())} overlaps per donor, every per-donor estimate is very imprecise.</p>

<h3>Step 3 · Do donors disagree? The heterogeneity test</h3>
<p>If all donors shared the pooled ratio R, donor d would be expected to show R·E<sub>d</sub> overlaps. The test sums the squared
deviations from that, Pearson style, over donors with E<sub>d</sub> ≥ 1 (the others are too small to say anything):</p>
<div class="formula">R = {pooled.observed} / {pooled.expected:.2f} = {Rr:.2f}<br>{x2_lines}<br>{f"X² = {X2:.2f}, df = {dfree} → p = {chi2.sf(X2, dfree):.2f}" if dfree > 0 else f"only {len(terms)} donor with E ≥ 1: df = 0, the test cannot be computed"}</div>
<p>{"With only " + str(len(terms)) + " informative donors the test has almost no power: absence of heterogeneity is not evidence of homogeneity." if dfree > 0 else "Here the per-donor expected counts are too small for the test; the per-donor ratios above are the only information on consistency between donors."}</p>

<h3>Step 4 · When A and C disagree: the gains</h3>
<p>For gains against the same set, test A gives fold {ga.fold:.2f} ({int(ga.observed)} / {ga.expected:.2f} events). Test C gives pooled
{g.loc['pooled'].ratio:.2f} ({g.loc['pooled'].observed} / {g.loc['pooled'].expected:.2f} breakends), and P08 alone {g.loc['P08'].ratio:.2f}
({g.loc['P08'].observed} / {g.loc['P08'].expected:.2f}, {g.loc['P08'].breakends} breakends in {g.loc['P08'].cells} cells). P08's gains are
dominated by a few large clones (the {big.cells}-cell event above among them) whose breakends carry a sizeable null probability but do not
overlap a site; counted once per cell, they add many times to E and nothing to O, swamp the other donors and pull the pooled ratio below 1. Pooling losses and gains, test C's heterogeneity p is {fp(al.het_p)}. This is the kind of discrepancy test C exists to
expose: a ratio driven by one donor's clonal structure, not by break location.</p>
<p><b>Change from the original plan.</b> The plan was a mixed logistic model (breakend overlaps ~ set, with a donor random effect).
With five donors the random-effect variance cannot be estimated, so the test reports per-donor ratios plus this heterogeneity test.</p>
{cards([
        "Is the enrichment present in each donor, or carried by one?",
        "Does the result change when every cell's breakend counts (i.e. weighting by clonal expansion)?"],
        ["Whether the break rate itself is higher at sites: clonal expansion and selection are mixed in.",
         "A formal population-level estimate: five donors are too few for a random-effects model."])}
<h3>Limitations of test C</h3>
<ul>
<li><b>Clones dominate.</b> One expanded clone can outweigh all other donors; the pooled ratio then describes that clone.</li>
<li><b>Dependence within clones.</b> The bootstrap resamples cells as if independent; cells of one clone share one event, so the CI is
too narrow when clones are large.</li>
<li><b>Few donors.</b> Five donors, of which often only one or two have E ≥ 1: the heterogeneity test has very little power.</li>
<li><b>Same null as test A.</b> It inherits the RT-only matching.</li>
</ul>
{flow([f"{R}/main/breakends_{BIN}kb.tsv", f"{R}/main/test_c.tsv", f"{R}/main/test_a.tsv"], ["scripts/docs/build_tests_explainer.py"],
           ["this section (not saved separately)"], "per-donor and pooled ratios, CIs and heterogeneity p read from test_c.tsv; the X² terms recomputed from those rows; clone sizes from breakends_&lt;bin&gt;kb.tsv.")}"""


def section_compare(ta, tb, tc):
    cases = [("macheret2020:HS68", 500, "loss"), ("genes_1mb", 200, "loss"), ("union_all:all", 500, "loss")]
    rows = []
    for name, b, d in cases:
        a = ta[(ta.feature_set == name) & (ta.bin_kb == b) & (ta.direction == d) & (ta.stratum == "main")].iloc[0]
        m = tb[(tb.feature_set == name) & (tb.bin_kb == b) & (tb.direction == d) & (tb.stratum == "main") & (tb.model == "adjusted")].iloc[0]
        c = tc[(tc.feature_set == name) & (tc.bin_kb == b) & (tc.direction == d) & (tc.donor == "pooled")].iloc[0]
        rows.append({"set · bin · direction": f"{name} · {b} kb · {d}",
                     "A: fold [CI] (O/E)": f"{a.fold:.2f} [{a.fold_lo:.2f}–{a.fold_hi:.2f}] ({int(a.observed)}/{a.expected:.1f})",
                     "A: p perm.": fp(a.p_perm), "A: q global": fp(a.q_global),
                     "B: adjusted IRR [CI]": f"{m.irr:.2f} [{m.irr_lo:.2f}–{m.irr_hi:.2f}]" if np.isfinite(m.irr) else "– (separation)", "B: q global": fp(m.q_global),
                     "C: pooled O/E [CI]": f"{c.ratio:.2f} [{c.ratio_lo:.2f}–{c.ratio_hi:.2f}]", "C: heterogeneity p": fp(c.het_p)})
    t = pd.DataFrame(rows)
    return f"""<h2 id="compare">5 · Reading the three tests together</h2>
<div class="tablebox"><table><thead><tr><th></th><th>Test A</th><th>Test B</th><th>Test C</th></tr></thead><tbody>
<tr><td>unit</td><td>event (clone counts once)</td><td>donor × 1 Mb window</td><td>breakend of a cell (clone counts per cell)</td></tr>
<tr><td>“hits the set”</td><td>event interval overlaps a site</td><td>event assigned (largest overlap) to a window that overlaps a site</td><td>breakend interval overlaps a site</td></tr>
<tr><td>chance defined by</td><td>re-placing the event on its chromosome, same RT quintile and width</td><td>regression on RT, GC, mappability, large genes, chromosome, donor</td><td>as test A, per breakend</td></tr>
<tr><td>effect size</td><td>fold = O / E</td><td>IRR</td><td>O / E per donor and pooled</td></tr>
<tr><td>main strength</td><td>exact, flexible RT matching; keeps event width</td><td>adjusts several confounders at once, incl. large genes</td><td>shows donor-to-donor consistency and clonal weight</td></tr>
<tr><td>main weakness</td><td>RT is the only covariate</td><td>1 Mb grid and binary overlap; sparse counts</td><td>clones dominate; 5 donors</td></tr>
</tbody></table></div>
<p>Three real cases (E2–E4):</p>
{table(t, list(t.columns))}
<div class="method"><b>What is counted.</b> One row per set, bin and direction; A and C on interstitial E2–E4 events / breakends; B the adjusted model (without genes300). <b>How.</b> Read from test_a.tsv (stratum main), test_b.tsv (stratum main, model adjusted) and test_c.tsv (donor pooled); q global = BH across all CFS sets. <b>Limitations.</b> The three tests share most of the data, so agreement between them is not three independent replications; genes_1mb is a control and has no q global.</div>
<p>{"<br>".join(f"<b>{r['set · bin · direction']}</b>: test A fold {r['A: fold [CI] (O/E)']} (p perm. {r['A: p perm.']}); test B adjusted IRR {r['B: adjusted IRR [CI]']}; test C pooled O/E {r['C: pooled O/E [CI]']}." for r in rows)}</p>
<p><b>Rule of thumb.</b> A agrees with B → the excess is not explained by large genes or GC. A agrees with C → it is not an artefact of
one donor or one clone. A disagrees with C → look at clone sizes. A disagrees with B → look at event widths and large genes.</p>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = a.res
    fs = feature_sets(CFS, GENES)
    ta = pd.read_csv(f"{res}/main/test_a.tsv", sep="\t")
    sh = pd.read_csv(f"{res}/main/test_a_shuffle_controls.tsv", sep="\t")
    tb = pd.read_csv(f"{res}/main/test_b.tsv", sep="\t")
    tc = pd.read_csv(f"{res}/main/test_c.tsv", sep="\t")
    ev, null, ex = test_a_example(res, fs)
    pz = p_zero_hits(res, fs)
    w, df, x, donors = test_b_example(res, fs)

    parts = [f"""<title>Enrichment tests explained</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,600&family=Public+Sans:wght@400;600&family=JetBrains+Mono&display=swap">
<style>{CSS}{EXTRA_CSS}</style><main>
<span class="eyebrow">CNV_cfs · pilot · hg19 · teaching document</span>
<h1>Tests A, B and C, step by step on real data</h1>
<p class="meta">Generated on {date.today().isoformat()} by <code>scripts/docs/build_tests_explainer.py</code>. Numbers are read from
<code>results/enrichment/pilot/main/</code> or recomputed with the test scripts' own functions. Results of all sets are in the
enrichment document (<code>docs/enrichment/enrichment_pilot.html</code>).</p>
<nav><ol><li><a href="#shared">The shared question</a></li><li><a href="#testa">Test A: positional permutation</a></li>
<li><a href="#testb">Test B: window regression</a></li><li><a href="#testc">Test C: per breakend and donor</a></li>
<li><a href="#compare">Reading them together</a></li></ol></nav>""",
             section_shared(), section_a(ev, ex, null, fs, ta, sh, pz), section_b(w, df, x, donors, tb, ta, ex),
             section_c(res, tc, ta), section_compare(ta, tb, tc), "</main>"]
    os.makedirs(a.out, exist_ok=True)
    out = f"{a.out}/tests_explainer.html"
    with open(out, "w") as fh:
        fh.write("\n".join(parts))
    print(out)


if __name__ == "__main__":
    main()

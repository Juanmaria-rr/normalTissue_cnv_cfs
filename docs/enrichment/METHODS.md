# Methods: copy-number breakends at common fragile sites (pilot)

This document describes how the pilot results were produced, from reads to the three
enrichment tests. The results themselves are in the concise report
(`docs/enrichment/enrichment_brief.html`) and the full report (`enrichment_pilot.html`).
Every step is a script in this repository and documents its inputs and outputs in its header.
`bash scripts/06_enrichment/run_all.sh` rebuilds the analysis from the per-cell copy-number
tables.

## 1. Data

- **Dataset.** Single-cell DNA-seq of normal breast epithelium, Lin et al. 2024 (*Nature*
  636:663), ENA project PRJNA1003661.
- **Genome build.** hg19.
- **Pilot.** Five donors (P07, P08, P18, P20, P36).
  - 9,181 cells: wells the authors' imaging scored as one good-quality cell
    (`qc_flag = Good`, `cell_flag = yes`), control wells excluded.
  - The full dataset has 49 donors and about 91,000 cells.

## 2. Reads to copy number

Run on AWS HealthOmics with the workflows in `workflows/` and the containers in `containers/`.

1. **Trimming.** R1 only; cutadapt removes the Tn5 mosaic end; reads under 20 bp are discarded.
2. **Alignment.** bwa mem to hg19. Duplicates are marked with samtools markdup, after
   concatenating the runs of each cell.
3. **Relative copy number.** QDNAseq at three bin sizes in parallel (200 kb, 500 kb, 1 Mb).
   - Reads with MAPQ ≥ 37 and not duplicated are counted.
   - Bins are filtered by `applyFilters`, then GC and mappability are corrected.
   - Segmentation is CBS with `undo.SD = 1`.
   - Pipeline: `external/sWGS_relativeCN`.
4. **Absolute fits.** `fitting_absCN` is run on cells with ≥ 20,000 usable reads. It is used
   only as a flag, because breakends come from relative copy number. 97.5% of the pilot's
   altered cells have a stable diploid fit.

## 3. Cell selection (per bin size)

A cell enters the main set at a given bin size when it passes all three cuts.
Script: `scripts/03_cn_qc/cn_qc_metrics.R`.

| Cut | Rule | Meaning |
|---|---|---|
| Depth (`depth_ok`) | mean corrected reads per usable bin / 2 ≥ 15 | enough reads to tell one copy from counting noise: 15 reads per bin per chromosome copy in a diploid cell |
| Noise (`noise_excess`) | ≤ 1.5 | robust SD of neighbour-bin differences divided by its Poisson expectation; 1 = pure counting noise |
| S phase (`rt_rho`) | ≤ 0.10 | Spearman correlation between the cell's relative copy number and replication timing; a high value means the cell was caught replicating |

Main-set cells: 2,052, 5,032 and 6,297 at 200 kb, 500 kb and 1 Mb.

The depth rule was recalibrated on 5 November 2026 (`scripts/03_cn_qc/recalibrate_depth.py`).
The former rule assumed that 91% of reads land in usable bins; the observed fraction is about
77% and depends on the donor.

## 4. Breakends and events

Script: `scripts/06_enrichment/build_event_table.py`.

- **Segment state.** Segment ratio = segment value / median over the cell's usable bins of their
  segmented values.
  - Gain if the ratio is > 1.2, loss if < 0.8, otherwise neutral.
  - Consecutive segments with the same state are merged.
  - Merged segments with less than 1 Mb of usable bins are dropped.
- **Breakend.** The boundary between consecutive retained segments of a chromosome with
  different state. Its interval `[lo, hi]` runs from the start of the last bin of the left
  segment to the end of the first bin of the right one, so it spans at least two bins.
- **Deviation and strata.** Deviation = largest |ratio − 1| of the two flanks.
  - E1: 0.2–0.3; E2: 0.3–0.4; E3: 0.4–0.6; E4: ≥ 0.6.
  - **Main definition: E2–E4**, which corresponds to the paper's ratio threshold of 1.3 / 0.7.
- **Direction.** Loss or gain, taken from the flank with the larger deviation.
- **Type.**
  - Arm-level if the interval overlaps a pericentromeric block ± 1 Mb (UCSC acen, gvar and
    stalk bands; the whole p arm of acrocentric chromosomes).
  - Otherwise interstitial. Only interstitial breakends are tested.
- **Terminal.** A breakend is terminal when its altered segment runs to the end of the
  chromosome. This flag is used only in the robustness split.
- **Event.** Within one donor, breakends with the same chromosome, direction and type are
  chained when each lies within 2 Mb of the previous one.
  - Each chain becomes one event, with the medians of its breakends' lo, hi and deviation.
  - A clone therefore counts once. The same locus in two donors counts as two events.

## 5. Site sets

Scripts: `scripts/05_features/`. All sets are in hg19; hg38 sources are lifted over with the UCSC
chain, dropping intervals that change chromosome or whose length changes by more than 20%.

Names follow `<catalogue>:<cell type>`, where `:all` is the union of the catalogue's cell types.

| Kind | Sets |
|---|---|
| Measured | Macheret 2020 MiDAS (HS68, HeLa, U2OS); Ji 2020 MiDAS (HBE, HeLa, U2OS); HumCFS cytogenetic bands; Le Tallec 2013 FISH (all, epithelial lines, MCF10A); Brison 2019 aphidicolin replication delay; Sarni 2020 fragility signature; Wilson 2015 aphidicolin CNV hotspots; Simpson 2021 cancer-SV hotspots |
| Derived | consensusA core (≥ 1 replication-class and ≥ 1 breakage-class source) / reproducible (≥ 2 cell types) / all; breast_predicted (genes > 300 kb, open in breast scATAC GSE272504 and late replicating); A_plus_C; union_all |
| Controls | genes > 300 kb and genes > 1 Mb (refGene loci) |

## 6. Enrichment tests

All tests are run per bin size and direction (loss, gain, all) on interstitial E2–E4 events.
Test A is also run per deviation stratum.

### Test A: positional permutation

Script: `test_a_permutation.py`. Unit: one event.

- **Observed.** O = number of events whose interval overlaps the set.
- **Null.** Each event is re-placed, with its width, at every usable-bin start of its own
  chromosome. A placement is kept when its replication timing falls in the event's
  replication-timing quintile and it does not touch a pericentromeric block ± 1 Mb.
  - Replication timing is computed the same way for placements and events: the mean over the
    usable bins the interval covers.
  - The event's null overlap probability p_i is the fraction of its placements that overlap
    the set. It is exact, because all placements are enumerated.
- **Expected and fold.** E = Σ p_i; fold = O / E.
- **p exact.** The Poisson-binomial tail P(S ≥ O), computed exactly. It assumes events are
  independent.
- **p perm (reference).** The site set is shuffled 5,000 times: each site moves to a random
  usable bin of the same chromosome and replication-timing quintile, keeping its width. The
  fold is recomputed each time, and p = (1 + number of shuffled folds ≥ observed) / 5,001.
  This keeps the clustering of events and sites.
- **Confidence intervals.** Bootstrap of events (2,000 resamples), and a bootstrap of donors as
  a rough guide.

### Test B: window regression

Script: `test_b_regression.py`.

- **Unit.** Donor × usable 1 Mb window outside pericentromeric blocks: 1,917 windows × 5 donors.
- **Assignment.** Each event goes to the eligible window with which its interval overlaps most.
  Events that overlap no eligible window are excluded and counted.
- **Model.** Poisson, or negative binomial when over-dispersed:
  `events ~ set overlap (yes/no) + RT + GC + mappability + chromosome + donor`, with offset
  log(usable fraction of the window).
  - A second model adds the fraction of the window in genes > 300 kb. Large genes are probably
    part of the fragility mechanism, so this model is conservative.
  - A crude model (set + donor) is also reported.
- **Effect.** IRR = exp(β) of the set term.
- **Standard error.** The larger of the model-based SE and the SE clustered by chromosome.
- **Separation.** When no event falls in windows of the set, no estimate is reported.

### Test C: per breakend and donor

Script: `test_c_per_donor.py`. Unit: one breakend of one cell, not collapsed, so a clone weighs
by its number of cells.

- o and p_i are defined as in test A.
- O/E is reported per donor and pooled. The 95% CI comes from a bootstrap of cells; when
  O = 0, the upper limit is the one-sided Poisson bound 3.0 / E.
- Heterogeneity between donors: X² over donors with E ≥ 1.

### Multiple testing

There is no pre-declared primary test, so all results are exploratory. Benjamini–Hochberg is
applied over all 198 main-stratum rows of the CFS sets (sets × bins × directions; controls
excluded):

- **q global** on the exact p;
- **q perm. global** on the permutation p;
- **q family** (BH within one set) is also given, but it does not control the error across sets.

## 7. Robustness and power

- **Variants of test A.** Script: `robustness.py`. Each variant re-runs test A:
  - no noise or S-phase cuts;
  - the same cells in all three bins;
  - the 500 kb main cells called at every bin (`cells500`), so bins differ only in resolution;
  - 500 kb breakends refined to 200 kb in the same cell (`refined`);
  - linkage at 1 Mb and 5 Mb instead of 2 Mb;
  - collapsing on E2–E4 breakends only;
  - leaving out each donor;
  - leaving out each chromosome;
  - terminal versus internal events.
- **Power.** Script: `power_scaleup.py`. Expected events are scaled to 49 donors (×9.05, the
  ratio of cells). The smallest fold detectable with 80% power is computed with a one-sided
  exact Poisson test, at α = 0.05 and at 0.05 / 198.

## 8. Known limitations

- Five donors and few events: the strongest row (HS68 MiDAS, 500 kb losses) rests on 10 events.
- Events cluster (clones, recurrent regions, donor-specific losses), so the exact p is
  optimistic and the permutation p is the reference.
- Breakend intervals (≥ 2 bins) are about as wide as the sites, so overlap measures proximity
  at the scale of the bin.
- CFS catalogues and replication timing come from other cell types. Only Le Tallec's three
  MCF10A sites are breast epithelium.
- Only breaks that leave a copy-number change are detected.

## 9. Code and outputs

| Step | Script | Output |
|---|---|---|
| QC metrics | `scripts/03_cn_qc/cn_qc_metrics.R`, `recalibrate_depth.py` | `results/pilot_cn/qc/cn_qc_<bin>kb.tsv` |
| Segments | `scripts/03_cn_qc/call_breakends.R` | `results/pilot_cn/qc/breakends_<bin>kb_segments.tsv` |
| Events | `scripts/06_enrichment/build_event_table.py` | `results/enrichment/pilot/<set>/events_<bin>kb.tsv` |
| Tests | `test_a_permutation.py`, `test_b_regression.py`, `test_c_per_donor.py` | `results/enrichment/pilot/main/test_{a,b,c}.tsv` |
| Robustness, power | `robustness.py`, `power_scaleup.py` | `robustness_*.tsv`, `power_scaleup.tsv` |
| Reports | `scripts/docs/build_enrichment_brief.py`, `build_enrichment_html.py`, `summary_figures.py` | `docs/enrichment/`, `docs/slides/figs/` |

Public copy of the code: https://github.com/Juanmaria-rr/normalTissue_cnv_cfs

# Manual curation (audit) of single-cell CN fits

Breakends are taken from the **relative** CN. Manual curation is an **audit** of ~100 pilot
cells: it checks that the mathematical fit matches the biological one, and that
ratio-based breakends agree with those from the curated integer solution (decision
2026-10-29, see `WORKLOG.md`).

The app is the group's curation app (Shiny for Python), vendored from the
`<curation-app delivery>` delivery. Only the code is copied, and the first commit is unmodified
(`6b69998`). Changes for single cells are listed at the top of `app/curation_app/app.py`:

- On sample change the app starts from the **unfiltered** optimum. With single cells the
  constraint filter can force a wrong low ploidy, because the powered threshold falls with
  ploidy. The filtered optimum is shown with a "⚠ differs" warning when it disagrees.
- Star rating and notes reset on sample change. Saved samples are marked ✓.
- A per-cell QC panel shows reads, nrpbpcc, noise_excess, breakpoints, frac_nondiploid and
  rt_rho.
- There is no downsampled panel for single cells. Saved rows carry a timestamp.

## Steps

```bash
# 0. one-time: virtualenv in curation/.venv (Python >= 3.11)
bash curation/setup_mac.sh

# 1. export relative CN of the selected cells (cnv-cfs/copynumber image)
Rscript scripts/04_curation/export_relcn.R --relcn <cn out>/relCN --bin 200 \
    --cells cells.txt --out curation/staging/<batch>

# 2. build the app batch: HDF5 + full ploidy x purity grid + QC panel + blinding key
curation/.venv/bin/python scripts/04_curation/build_curation_batch.py \
    --export curation/staging/<batch> --bin 200 --batch <batch> \
    --out-dir curation/data --key-dir curation/keys --qc <cn_qc_200kb.tsv>

# 3. curate (http://127.0.0.1:8501); one batch per app run
bash curation/run_curation.sh <batch>

# 4. unblind -> <batch>_curated.tsv + bestfitting.txt (input of bestsolution_absCN.R)
curation/.venv/bin/python scripts/04_curation/unblind_curated.py --batch <batch> \
    --data-dir curation/data --key-dir curation/keys --out results/curation
```

## How to judge a fit

Integer CN states should sit on the horizontal gridlines. Start from the unfiltered optimum.

- **2 vs 4:** a flat profile fits equally well at ploidy 2, 3 or 4. Prefer 2 unless a
  segment sits clearly at an odd state relative to the rest.
- **Star rating** is confidence in the fit, not sample quality: 3 = unambiguous,
  0 = no fit looks right.
- **Notes:** write anything odd (wave-like profile, many short segments, half-integer
  segments).
- **Blinding:** sample codes are random, and the key lives in `curation/keys/`, outside the
  folder the app opens. Do not look up the key while curating.

## Grid

`build_curation_batch.py` recomputes the grid with the exact formulas of
`pureSamples_absCN.R`: ploidy 1.4–8, purity 0.90–1.00, the same clonality, `powered` and
`LOHok`. On the 20 test cells it matches upstream `resfitting.txt` to 6e-15 in clonality,
with 0 `powered` mismatches over 1,340 cell × ploidy pairs. It keeps every purity, whereas
upstream keeps the best purity per ploidy. As a result the app's "filtered optimum" can
differ from upstream's choice: for P36_chip1_co_r51c33 the app gives 1.4, upstream 1.6 and
the unfiltered optimum 2.0.

## Layout

```
curation/app/       vendored app code (curation_app/ + the pylib modules it imports)
curation/data/      batches (git-ignored; rebuilt from the scripts)
curation/keys/      blinding keys (versioned)
curation/staging/   long-format exports (git-ignored)
```

Known limitation of the vendored plot: the "σ̂Δ*" noise in the profile title is computed on
segmented values and reads ~0 for single cells. Use `noise_excess` in the QC panel instead.

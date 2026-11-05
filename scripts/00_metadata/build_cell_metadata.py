#!/usr/bin/env python3
"""Build per-cell and per-donor metadata tables for PRJNA1003661.

Sources:
  - ENA run report (one SRA run = one single cell).
  - Paper code repo (navinlabcode/normalbreastDNA): co-assay cell -> ATAC barcode map.

Outputs (in --out-dir):
  cells.tsv   one row per WGS run (cell_id groups runs of the same cell)
  donors.tsv  per-donor summary (cell counts, QC flags, depth vs nrpbpcc thresholds)

Usage:
  python3 scripts/00_metadata/build_cell_metadata.py --out-dir data/metadata
"""
import argparse
import csv
import io
import os
import re
import statistics
import subprocess
from collections import defaultdict

PROJECT = "PRJNA1003661"
ENA_FIELDS = [
    "run_accession", "sample_alias", "sample_accession", "experiment_accession",
    "library_name", "library_strategy", "library_layout", "instrument_model",
    "read_count", "base_count", "tissue_type", "fastq_ftp",
]
ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}"
           "&result=read_run&fields={fields}&format=tsv&limit=0")
COASSAY_URL = ("https://raw.githubusercontent.com/navinlabcode/normalbreastDNA/main/"
               "ATAC_analysis/HBCA_aneuploidy_coassay_barcode_cellname_dict.tsv")

# Usable autosomal bins per bin size after the sWGS_relativeCN filters; computed by
# scripts/00_metadata/count_usable_bins.R (fitting_absCN's nbins_ref_genome).
USABLE_BINS = {200: 12017, 500: 4407, 1000: 1942}
NRPBPCC = 15
PLOIDY = 2
USED_FRACTION = 0.91   # same constant as fitting_absCN

# library_name is not uniformly formatted across sequencing batches, e.g.
#   BCMHBCA25L_coassay_r47_c41_yes_B1_Good
#   Sample_BCMHBCA09L_r0c9_yes_A2_Good
#   chip2_BCMHBCA..._r3_c5_yes_C1_Good_S12
#   BCMHBCA22L_r31c22_yes_Good_1
#   BCMHBCA75L_r49_c38_yes_B1_ManualUse_MultipleCells_S3203_R1_001
#   BCMHBCA10R_coassay_r7_c2_yes_P24_NoCells      (A24/P24 = control wells)
# so it is parsed token by token rather than with one anchored regex.
TISSUE_RE = re.compile(r"BCMHBCA\d+[LR]")
WELL_RE = re.compile(r"_r(\d+)_?c(\d+)_")
CHIP_RE = re.compile(r"(?:^|_)(chip\d+)(?:_|$)")
FLAG_RE = re.compile(r"_(yes|no)_(?:([A-Z]\d{1,2}|W\d|NA)_)?")
SUFFIX_RE = re.compile(r"(_S\d+)?(_R\d_001)?$")
QC_TOKENS = ("Good", "ManualUse", "NoCells", "CellsInControl", "LowConfidence",
             "MultipleCells", "Cluster")
CONTROL_WELLS = {"A24", "P24"}


def parse_library_name(name):
    """Return dict of parsed fields; missing pieces are empty strings."""
    out = {"tissue": "", "coassay": 0, "chip": "chip1", "row": "", "col": "",
           "cell_flag": "", "quadrant": "", "qc": ""}
    if m := TISSUE_RE.search(name):
        out["tissue"] = m.group(0)
    out["coassay"] = int("_coassay_" in name)
    if m := CHIP_RE.search(name):
        out["chip"] = m.group(1)
    if m := WELL_RE.search(name):
        out["row"], out["col"] = m.groups()
    if m := FLAG_RE.search(name):
        out["cell_flag"], out["quadrant"] = m.group(1), m.group(2) or ""
    tail = SUFFIX_RE.sub("", name)
    toks = [t for t in tail.split("_") if any(t.startswith(q) for q in QC_TOKENS)]
    out["qc"] = parse_qc("_".join(toks)) if toks else "unparsed"
    return out


def fetch(url):
    # curl instead of urllib: python.org macOS builds often lack a CA bundle
    return subprocess.run(["curl", "-sSfL", "--retry", "3", url],
                          check=True, capture_output=True, text=True).stdout


def threshold_reads(bin_kb):
    return NRPBPCC * PLOIDY * USABLE_BINS[bin_kb] / USED_FRACTION


def parse_qc(raw):
    # normalise "ManualUseLowConfidence" -> "ManualUse_LowConfidence"
    return re.sub(r"^ManualUse(?=[A-Z])", "ManualUse_", raw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    runs = list(csv.DictReader(io.StringIO(
        fetch(ENA_URL.format(acc=PROJECT, fields=",".join(ENA_FIELDS)))), delimiter="\t"))
    coassay = {r["Sample_ID"]: r["barcode"] for r in csv.DictReader(
        io.StringIO(fetch(COASSAY_URL)), delimiter="\t")}

    cells, unparsed = [], 0
    for r in runs:
        if r["library_strategy"] != "WGS":
            continue          # drops the Visium RNA-seq runs
        g = parse_library_name(r["library_name"])
        if g["qc"] == "unparsed" or not g["row"]:
            unparsed += 1
        n_reads = int(r["read_count"] or 0)
        n_files = max(1, len(r["fastq_ftp"].split(";")))
        rec = {
            "run_accession": r["run_accession"],
            "donor": r["sample_alias"],
            "tissue_id": g["tissue"],
            "library_name": r["library_name"],
            "coassay": g["coassay"],
            "atac_barcode": coassay.get(r["library_name"], ""),
            "chip": g["chip"],
            "row": g["row"], "col": g["col"],
            "cell_flag": g["cell_flag"],
            "quadrant": g["quadrant"],
            "control_well": int(g["quadrant"] in CONTROL_WELLS),
            "qc_flag": g["qc"],
            "instrument": r["instrument_model"],
            "layout": r["library_layout"],
            "read_count": n_reads,                         # spots = fragments
            "read_length": round(int(r["base_count"] or 0) / max(1, n_reads) / n_files),
            "tissue_type": r["tissue_type"],
            "fastq_ftp": r["fastq_ftp"],
        }
        cells.append(rec)

    # Some libraries were sequenced in two runs (top-up on another instrument; same
    # library_name and well). One cell = one well; runs of the same cell are merged
    # downstream (concatenated FASTQs), so depth flags use the per-cell total.
    for c in cells:
        c["cell_id"] = (f"{c['donor']}_{c['chip']}_{'co' if c['coassay'] else 'dna'}"
                        f"_r{c['row']}c{c['col']}")
    per_cell = defaultdict(list)
    for c in cells:
        per_cell[c["cell_id"]].append(c)
    for c in cells:
        runs_c = per_cell[c["cell_id"]]
        c["runs_per_cell"] = len(runs_c)
        c["cell_read_count"] = sum(x["read_count"] for x in runs_c)
        for kb in USABLE_BINS:
            c[f"pass_depth_{kb}kb"] = int(c["cell_read_count"] >= threshold_reads(kb))

    out = args.out_dir.rstrip("/")
    os.makedirs(out, exist_ok=True)
    with open(f"{out}/cells.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cells[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(cells)

    by_donor = defaultdict(list)
    for c in cells:
        by_donor[c["donor"]].append(c)
    rows = []
    for d, cs in sorted(by_donor.items()):
        good_runs = [c for c in cs if c["qc_flag"] == "Good" and c["cell_flag"] == "yes"
                     and not c["control_well"]]
        good = list({c["cell_id"]: c for c in good_runs}.values())   # one row per cell
        row = {
            "donor": d,
            "tissue_id": cs[0]["tissue_id"],
            "n_runs": len(cs),
            "n_cells": len({c["cell_id"] for c in cs}),
            "n_good_single": len(good),
            "n_good_multi_run": sum(c["runs_per_cell"] > 1 for c in good),
            "n_coassay": len({c["cell_id"] for c in cs if c["coassay"]}),
            "median_reads_good": int(statistics.median([c["cell_read_count"] for c in good])) if good else 0,
            "instruments": ",".join(sorted({c["instrument"] for c in cs})),
            "read_lengths": ",".join(str(x) for x in sorted({c["read_length"] for c in cs})),
        }
        for kb in USABLE_BINS:
            row[f"n_good_pass_{kb}kb"] = sum(c[f"pass_depth_{kb}kb"] for c in good)
        rows.append(row)
    with open(f"{out}/donors.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"runs: {len(cells)}  cells: {len(per_cell)}  unparsed library_name: {unparsed}"
          f"  donors: {len(rows)}")
    print("depth thresholds (fragments): " + ", ".join(
        f"{kb}kb={int(threshold_reads(kb))}" for kb in USABLE_BINS))


if __name__ == "__main__":
    main()

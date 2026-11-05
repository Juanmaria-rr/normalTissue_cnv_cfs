#!/usr/bin/env python3
"""Inventory of all sequencing data released with Lin et al. 2024 (normal breast scDNA).

Sources (per the paper's Data availability statement):
  - PRJNA1003661 (SRA): single-cell DNA (WGS strategy; DNA-only and DNA-ATAC co-assay
    wells) and Visium spatial transcriptomics (RNA-Seq strategy).
  - GSE272504 / PRJNA1136981 (GEO/SRA): scATAC part of the DNA-ATAC co-assay, one pooled
    run per donor; cells are linked to their DNA well by the ATAC barcode in cells.tsv.
No single-cell RNA-seq was generated in the paper.

Per-cell scDNA counts come from data/metadata/cells.tsv (build_cell_metadata.py), so the
cell definitions (cell_id, Good/yes/non-control) are the same as everywhere else.

Outputs (in --out-dir):
  data_inventory_datasets.tsv   one row per data type
  data_inventory_donors.tsv     one row per donor, data types as columns

Usage:
  python3 scripts/00_metadata/build_data_inventory.py --cells data/metadata/cells.tsv \
      --pilot data/pilot/pilot_summary.tsv --out-dir data/metadata
"""
import argparse
import csv
import io
import os
import re
import subprocess
from collections import defaultdict

ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}"
           "&result=read_run&fields={fields}&format=tsv&limit=0")
FIELDS = ["run_accession", "library_strategy", "library_layout", "instrument_model",
          "read_count", "base_count", "sample_title", "library_name"]
DNA_PROJECT, ATAC_PROJECT = "PRJNA1003661", "PRJNA1136981"


def fetch_runs(acc):
    # curl instead of urllib: python.org macOS builds often lack a CA bundle
    txt = subprocess.run(["curl", "-sSfL", "--retry", "3",
                          ENA_URL.format(acc=acc, fields=",".join(FIELDS))],
                         check=True, capture_output=True, text=True).stdout
    return list(csv.DictReader(io.StringIO(txt), delimiter="\t"))


def donor_of(sample_title):
    m = re.match(r"(P\d+)", sample_title)
    return m.group(1) if m else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--pilot", help="pilot_summary.tsv, to flag pilot donors")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with open(args.cells) as f:
        runs = list(csv.DictReader(f, delimiter="\t"))
    pilot = set()
    if args.pilot:
        with open(args.pilot) as f:
            pilot = {r["donor"] for r in csv.DictReader(f, delimiter="\t")}

    # scDNA, per cell (first run of each cell carries the cell-level fields)
    cells = {}
    for r in runs:
        cells.setdefault(r["cell_id"], r)
    good = lambda c: c["qc_flag"] == "Good" and c["cell_flag"] == "yes" and c["control_well"] == "0"

    donors = defaultdict(lambda: defaultdict(int))
    for c in cells.values():
        kind = "coassay" if c["coassay"] == "1" else "dnaonly"
        d = donors[c["donor"]]
        d[f"scdna_{kind}_cells"] += 1
        d[f"scdna_{kind}_cells_good"] += good(c)
        if kind == "coassay" and c["atac_barcode"]:
            d["coassay_cells_with_atac_barcode"] += 1
    for r in runs:
        donors[r["donor"]]["scdna_runs"] += 1

    dna_runs = fetch_runs(DNA_PROJECT)
    atac_runs = fetch_runs(ATAC_PROJECT)
    visium = [r for r in dna_runs if r["library_strategy"] == "RNA-Seq"]
    other = [r for r in dna_runs if r["library_strategy"] not in ("WGS", "RNA-Seq")]
    assert not other, f"unexpected strategies: {set(r['library_strategy'] for r in other)}"
    assert sum(r["library_strategy"] == "WGS" for r in dna_runs) == len(runs), \
        "cells.tsv is out of date with ENA"

    def attach(rs, prefix):
        for r in rs:
            d = donors[donor_of(r["sample_title"])]
            d[f"{prefix}_runs"] += 1
            d[f"{prefix}_reads"] += int(r["read_count"] or 0)
            d[f"{prefix}_accessions"] = ",".join(filter(None, [d.get(f"{prefix}_accessions") or "",
                                                               r["run_accession"]]))
    attach(atac_runs, "atac")
    attach(visium, "visium")

    os.makedirs(args.out_dir, exist_ok=True)
    cols = ["donor", "pilot", "scdna_runs",
            "scdna_dnaonly_cells", "scdna_dnaonly_cells_good",
            "scdna_coassay_cells", "scdna_coassay_cells_good", "coassay_cells_with_atac_barcode",
            "atac_runs", "atac_reads", "atac_accessions",
            "visium_runs", "visium_reads", "visium_accessions"]
    with open(os.path.join(args.out_dir, "data_inventory_donors.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for donor in sorted(donors, key=lambda x: int(x[1:])):
            d = donors[donor]
            w.writerow([donor, int(donor in pilot)] + [d.get(c, 0) for c in cols[2:]])

    def summ(rs):
        return (len(rs), len({donor_of(r["sample_title"]) for r in rs}),
                sum(int(r["read_count"] or 0) for r in rs),
                ";".join(sorted({r["library_layout"] for r in rs})),
                ";".join(sorted({r["instrument_model"] for r in rs})))

    dna_by_kind = defaultdict(list)
    for r in runs:
        dna_by_kind["coassay" if r["coassay"] == "1" else "dnaonly"].append(r)
    rows = []
    for kind, label, note in [
        ("dnaonly", "scDNA-seq (DNA-only wells)", "one SRA run per well; some wells have 2 runs (top-up)"),
        ("coassay", "scDNA-seq (DNA part of DNA-ATAC co-assay)", "one SRA run per well; ATAC barcode links the cell to its ATAC profile"),
    ]:
        rs = dna_by_kind[kind]
        ncell = len({r["cell_id"] for r in rs})
        ngood = len({r["cell_id"] for r in rs if good(r)})
        rows.append([label, "SRA", DNA_PROJECT, "WGS", len(rs), len({r["donor"] for r in rs}),
                     ncell, ngood, sum(int(r["read_count"] or 0) for r in rs),
                     ";".join(sorted({r["layout"] for r in rs})),
                     ";".join(sorted({r["instrument"] for r in rs})), note])
    n, nd, reads, lay, ins = summ(atac_runs)
    rows.append(["scATAC-seq (ATAC part of DNA-ATAC co-assay)", "GEO/SRA", f"GSE272504 / {ATAC_PROJECT}",
                 "ATAC-seq", n, nd, "", "", reads, lay, ins,
                 "one pooled run per donor; demultiplex by cell barcode; paper: 31,176 epithelial cells, 19 donors"])
    n, nd, reads, lay, ins = summ(visium)
    rows.append(["Visium spatial transcriptomics", "SRA", DNA_PROJECT, "RNA-Seq", n, nd, "", "", reads,
                 lay, ins, "one run per section; H&E/tissue images are not in SRA"])
    with open(os.path.join(args.out_dir, "data_inventory_datasets.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["data_type", "repository", "accession", "library_strategy", "runs", "donors",
                    "cells", "cells_good_single_noncontrol", "reads", "layouts", "instruments", "note"])
        w.writerows(rows)
    for r in rows:
        print("\t".join(map(str, r[:9])))


if __name__ == "__main__":
    main()

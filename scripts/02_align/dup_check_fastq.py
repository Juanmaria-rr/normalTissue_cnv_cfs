#!/usr/bin/env python3
"""Sequence-level duplicate rates from raw ENA FASTQs, to explain per-donor duplicate gaps.

For each selected cell, downloads its runs' FASTQs (ENA fastq_ftp) into --raw and computes
the fraction of reads whose key was already seen:
  r1        first --k bases of R1 (approximates single-end, R1-only marking)
  pair      first --k bases of R1 + first --k bases of R2 (approximates paired marking)
per run and for the runs concatenated. Sequence identity is a proxy for position identity:
it misses duplicates with sequencing errors in the first k bases and does not see
adapter-trimmed or unmapped reads, so absolute values differ from samtools markdup; the
comparison between keys, runs and donors is what this check is for.

Usage: python3 scripts/02_align/dup_check_fastq.py --cells data/metadata/cells.tsv \
    --donor P20 --assay 0 --n 4 --raw data/pilot/dup_check/raw --out data/pilot/dup_check/dup_check_P20.tsv
"""
import argparse
import gzip
import os
import subprocess

import pandas as pd


def fetch(url, path):
    if not os.path.exists(path):
        subprocess.run(["curl", "-sSLf", "-o", path, "https://" + url if not url.startswith("http") else url], check=True)
    return path


def reads(path, k):
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().strip(); f.readline(); f.readline()
            yield s[:k]


def dup_rates(pairs, k):
    seen_r1, seen_pair = set(), set()
    n = d1 = dp = 0
    for r1, r2 in pairs:
        n += 1
        if r1 in seen_r1:
            d1 += 1
        else:
            seen_r1.add(r1)
        key = r1 + "|" + (r2 or "")
        if key in seen_pair:
            dp += 1
        else:
            seen_pair.add(key)
    return n, d1 / max(n, 1), dp / max(n, 1)


def run_pairs(files, k):
    r1 = reads(files[0], k)
    r2 = reads(files[1], k) if len(files) > 1 else None
    for a in r1:
        yield a, (next(r2) if r2 is not None else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--assay", type=int, required=True)
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.raw, exist_ok=True); os.makedirs(os.path.dirname(a.out), exist_ok=True)
    c = pd.read_csv(a.cells, sep="\t", low_memory=False)
    c = c[(c.donor == a.donor) & (c.coassay == a.assay) & (c.qc_flag == "Good")]
    ids = c.cell_id.drop_duplicates().sample(a.n, random_state=a.seed)
    rows = []
    for cid in ids:
        runs = c[c.cell_id == cid].sort_values("run_accession")
        all_pairs = []
        for r in runs.itertuples():
            urls = [u for u in str(r.fastq_ftp).split(";") if u]
            urls = [u for u in urls if u.endswith("_1.fastq.gz") or u.endswith("_2.fastq.gz")] or urls[:1]
            files = [fetch(u, os.path.join(a.raw, os.path.basename(u))) for u in sorted(urls)]
            pr = list(run_pairs(files, a.k))
            n, d1, dp = dup_rates(pr, a.k)
            rows.append(dict(cell_id=cid, run=r.run_accession, instrument=r.instrument, layout=r.layout,
                             reads=n, dup_r1=round(d1, 4), dup_pair=round(dp, 4) if len(files) > 1 else None))
            all_pairs += pr
        if len(runs) > 1:
            n, d1, dp = dup_rates(all_pairs, a.k)
            rows.append(dict(cell_id=cid, run="combined", instrument="", layout="", reads=n,
                             dup_r1=round(d1, 4), dup_pair=round(dp, 4)))
        print(cid, "done")
    t = pd.DataFrame(rows)
    t.to_csv(a.out, sep="\t", index=False)
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()

"""Shared genome annotations for the enrichment analysis (hg19, autosomes).

PERICEN        approximate pericentromeric blocks (centromere + adjacent heterochromatin;
               whole p arm for acrocentrics), Mb, from UCSC cytoBand acen/gvar/stalk bands
usable_bins()  bins kept by applyFilters (use == TRUE), from any export_relcn.R output
rt_per_bin()   mean replication timing (50 kb annotation of sWGS_relativeCN) per bin
"""
import numpy as np
import pandas as pd

PERICEN = {1: (121.5, 142.6), 2: (90.5, 96.0), 3: (87.9, 93.9), 4: (48.2, 52.7), 5: (46.1, 50.7),
           6: (58.5, 62.6), 7: (57.9, 61.9), 8: (43.1, 48.1), 9: (38.0, 71.0), 10: (38.0, 42.9),
           11: (51.0, 55.8), 12: (33.3, 38.2), 13: (0, 19.5), 14: (0, 19.1), 15: (0, 20.7),
           16: (34.6, 46.5), 17: (22.2, 25.8), 18: (15.0, 18.9), 19: (24.2, 28.6), 20: (25.6, 29.8),
           21: (0, 14.3), 22: (0, 16.3)}
AUTOSOMES = [str(i) for i in range(1, 23)]


def overlaps_pericen(chrom, lo, hi, pad_mb=1.0):
    bl, bh = PERICEN[int(chrom)]
    return (lo / 1e6 <= bh + pad_mb) and (hi / 1e6 >= bl - pad_mb)


def usable_bins(profiles_tsv_gz):
    # use flags are bin-level (applyFilters), identical for every cell: take the first cell
    first = None
    rows = []
    for chunk in pd.read_csv(profiles_tsv_gz, sep="\t", dtype={"chromosome": str}, chunksize=200000):
        if first is None:
            first = chunk.cell_id.iloc[0]
        c = chunk[chunk.cell_id == first]
        rows.append(c)
        if (chunk.cell_id != first).any():
            break
    b = pd.concat(rows)
    b = b[b.use & b.chromosome.isin(AUTOSOMES)][["chromosome", "start", "end"]].reset_index(drop=True)
    return b


def rt_per_bin(rt_file, bins, bin_bp):
    rt = pd.read_csv(rt_file, sep="\t", usecols=["seqnames", "start", "replicationTime"])
    rt["seqnames"] = rt.seqnames.astype(str)
    rt["bstart"] = ((rt.start - 1) // bin_bp) * bin_bp + 1
    m = rt.groupby(["seqnames", "bstart"]).replicationTime.mean()
    return m.reindex(list(zip(bins.chromosome, bins.start))).values


def interval_rt(bins_rt, chrom, lo, hi):
    # mean RT of usable bins overlapping [lo, hi]; NaN if none
    b = bins_rt[(bins_rt.chromosome == chrom) & (bins_rt.end >= lo) & (bins_rt.start <= hi)]
    return float(np.nanmean(b.rt)) if len(b) and np.isfinite(b.rt).any() else np.nan

#!/usr/bin/env Rscript
# Per-cell copy-number QC metrics for one bin size, from the copynumber workflow outputs.
#
# Inputs (layout as published by workflows/copynumber):
#   <relcn>/<donor>_<bin>kb/copyNumbersSegmented_<bin>kb.rds   relative CN, read-count space
#   <abscn>/<donor>_<bin>kb/absoluteCN.rds                     chosen fit (pData: ploidy, purity)
#   <abscn>/<donor>_<bin>kb/<cell>_resfitting.txt              candidate fits (powered, clonality)
#   --rt: bins_with_replicationTime_50kb.tsv (sWGS_relativeCN), high value = early replication
#
# Metrics (autosomes, bins with use == TRUE only):
#   reads            QDNAseq total.reads (MAPQ >= 37, non-duplicate reads counted in bins)
#   nrpbpcc          mean corrected reads per used bin / ploidy (reads per bin per chromosome copy)
#   depth_ok         reads >= 15 x 2 x used_bins / 0.91: the 15-nrpbpcc requirement at the DIPLOID
#                    threshold (breakends come from relative CN, so depth is not tied to the fit)
#   ploidy, clonality  fit chosen by pureSamples_absCN.R --filters powered (kept for reference)
#   ploidy_best, clonality_best  minimum-clonality fit with no filter, from <cell>_resfitting.txt;
#                    clonality = mean squared distance to integer CN (lower = cleaner)
#   fit_disagree     1 if ploidy_best differs from the filtered choice (cell near the depth
#                    threshold; the filter can force a wrong low ploidy)
#   ploidy_flag      1 if ploidy_best >= 2.6 (not near-diploid: review or exclude from ratio calls)
#   noise_excess     robust SD (MAD) of neighbour differences (x[i+1] - x[i]) / sqrt(seg[i+1] + seg[i]),
#                    within chromosome. x = corrected counts, seg = segment mean (read-count space).
#                    Pure Poisson counting noise gives 1; values > 1 are noise beyond depth.
#   breakpoints      changes of integer CN (rounded absolute segments) between adjacent used bins,
#                    within chromosomes
#   frac_nondiploid  fraction of used bins with integer CN != 2 (descriptive, not an aneuploidy call)
#   rt_rho           Spearman rho between relative CN (x / median x) and replication timing of the
#                    bin (50 kb RT averaged into the bin). S-phase cells: positive (early > late)
#
# Usage:
#   Rscript scripts/03_cn_qc/cn_qc_metrics.R --relcn <dir> --abscn <dir> --bin 200 \
#       --rt external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv \
#       [--cells data/metadata/cells.tsv] --out <table.tsv>
# Runs in the cnv-cfs/copynumber image (Biobase, QDNAseqmod).

suppressPackageStartupMessages({
    library(argparse)
    library(Biobase)
    library(QDNAseqmod)
})

p <- ArgumentParser()
p$add_argument("--relcn", required = TRUE)
p$add_argument("--abscn", required = TRUE)
p$add_argument("--bin", type = "integer", required = TRUE)
p$add_argument("--rt", required = TRUE)
p$add_argument("--cells", help = "cells.tsv, to add donor/assay/instrument columns")
p$add_argument("--out", required = TRUE)
a <- p$parse_args()

# Replication timing averaged into bins of the requested size
rt <- read.delim(a$rt)
bin_bp <- a$bin * 1000
rt$key <- paste0(rt$seqnames, ":", floor((rt$start - 1) / bin_bp) * bin_bp + 1)
rt_bin <- tapply(rt$replicationTime, rt$key, mean, na.rm = TRUE)

robust_sd <- function(v) mad(v, na.rm = TRUE)

cell_metrics <- function(rel, abscn, resdir) {
    fd <- fData(rel)
    keep <- fd$use & fd$chromosome %in% as.character(1:22)
    chr <- fd$chromosome[keep]
    key <- paste0(fd$chromosome, ":", fd$start)[keep]
    x_all <- assayDataElement(rel, "copynumber")[keep, , drop = FALSE]
    s_all <- assayDataElement(rel, "segmented")[keep, , drop = FALSE]
    a_all <- assayDataElement(abscn, "segmented")
    a_all <- a_all[match(rownames(x_all), rownames(a_all)), , drop = FALSE]
    rtv <- rt_bin[key]
    pd <- pData(abscn)

    do.call(rbind, lapply(colnames(abscn), function(cell) {
        x <- x_all[, cell]; s <- s_all[, cell]; icn <- round(a_all[, cell])
        same <- c(chr[-1] == chr[-length(chr)])
        d <- (x[-1] - x[-length(x)]) / sqrt(pmax(s[-1] + s[-length(s)], 1e-9))
        brk <- sum(same & (icn[-1] != icn[-length(icn)]), na.rm = TRUE)
        res <- read.delim(file.path(resdir, paste0(cell, "_resfitting.txt")))
        fit <- res[which.min(abs(res$ploidy - pd[cell, "ploidy"]) + abs(res$purity - pd[cell, "purity"])), ]
        data.frame(
            cell_id = cell,
            bin_kb = a$bin,
            reads = rel$total.reads[colnames(rel) == cell],
            used_bins = length(x),
            ploidy = pd[cell, "ploidy"],
            nrpbpcc = round(mean(x, na.rm = TRUE) / pd[cell, "ploidy"], 2),
            depth_ok = as.integer(rel$total.reads[colnames(rel) == cell] >= 15 * 2 * length(x) / 0.91),
            clonality = round(fit$clonality, 4),
            ploidy_best = res$ploidy[which.min(res$clonality)],
            clonality_best = round(min(res$clonality), 4),
            fit_disagree = as.integer(abs(res$ploidy[which.min(res$clonality)] - pd[cell, "ploidy"]) > 1e-6),
            ploidy_flag = as.integer(res$ploidy[which.min(res$clonality)] >= 2.6),
            noise_excess = round(robust_sd(d[same]), 3),
            breakpoints = brk,
            frac_nondiploid = round(mean(icn != 2, na.rm = TRUE), 4),
            rt_rho = round(suppressWarnings(cor(x / median(x, na.rm = TRUE), rtv,
                                                method = "spearman", use = "complete.obs")), 3)
        )
    }))
}

dirs <- list.dirs(a$abscn, recursive = FALSE)
dirs <- dirs[grepl(paste0("_", a$bin, "kb$"), basename(dirs))]
out <- do.call(rbind, lapply(dirs, function(d) {
    tag <- basename(d)
    rel <- readRDS(file.path(a$relcn, tag, paste0("copyNumbersSegmented_", a$bin, "kb.rds")))
    abscn <- readRDS(file.path(d, "absoluteCN.rds"))
    # pureSamples_absCN.R replaces "(" and ")" in sample names; align relative names likewise
    colnames(rel) <- gsub("[()]", "_", colnames(rel))
    cell_metrics(rel[, colnames(abscn)], abscn, d)
}))

if (!is.null(a$cells)) {
    meta <- read.delim(a$cells, colClasses = "character")
    meta <- meta[!duplicated(meta$cell_id), c("cell_id", "donor", "coassay", "instrument", "layout", "read_length")]
    out <- merge(meta, out, by = "cell_id", all.y = TRUE)
}
write.table(out, a$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(out), "cells ->", a$out, "\n")

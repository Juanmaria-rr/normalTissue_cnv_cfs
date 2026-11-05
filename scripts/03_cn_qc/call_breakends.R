#!/usr/bin/env Rscript
# Ratio-based CN states, segments and breakends per cell (decision 2026-10-29: breakends come
# from the relative CN, not from the absolute fit).
#
# For each cell (autosomes, bins with use == TRUE):
#   rel     = segment mean / median of the cell's segmented values (robust to aneuploid segments)
#   state   = gain if rel > --gain, loss if rel < --loss, else neutral
#             (defaults 1.3 / 0.7, the ratio thresholds of Lin et al. 2024; a one-copy step
#             in a diploid cell is 1.5 / 0.5)
#   segments shorter than --min-seg-mb are dropped before calling breakends, so short
#             isolated segments do not add pairs of breakends
#   breakend = boundary between consecutive retained segments on the same chromosome with a
#             different state; pos = midpoint between them; [lo, hi] = from the start of the last
#             bin of the left segment to the end of the first bin of the right one (one bin
#             either side, plus any gap of filtered or dropped bins)
#
# Inputs:  <relcn>/<donor>_<bin>kb/copyNumbersSegmented_<bin>kb.rds
# Outputs: <out>_segments.tsv   cell_id, chromosome, start, end, n_bins, rel, state
#          <out>_breakends.tsv  cell_id, chromosome, pos, lo, hi, left_state, right_state, left_rel, right_rel
#          <out>_cells.tsv      cell_id, n_segments, n_breakends, n_gain_seg, n_loss_seg,
#                               frac_altered (fraction of retained bins not neutral), altered_mb
#
# Usage (cnv-cfs/copynumber image):
#   Rscript scripts/03_cn_qc/call_breakends.R --relcn <cn out>/relCN --bin 200 --out results/breakends/pilot_200kb
#   [--gain 1.3 --loss 0.7 --min-seg-mb 1]

suppressPackageStartupMessages({
    library(argparse)
    library(Biobase)
    library(QDNAseqmod)
})

p <- ArgumentParser()
p$add_argument("--relcn", required = TRUE)
p$add_argument("--bin", type = "integer", required = TRUE)
p$add_argument("--gain", type = "double", default = 1.3)
p$add_argument("--loss", type = "double", default = 0.7)
p$add_argument("--min-seg-mb", type = "double", default = 1)
p$add_argument("--out", required = TRUE)
a <- p$parse_args()
options(scipen = 999)
bin_bp <- a$bin * 1000

call_cell <- function(cell, fd, seg) {
    k <- fd$use & fd$chromosome %in% as.character(1:22) & is.finite(seg)
    f <- fd[k, ]; s <- seg[k]
    rel_all <- s / median(s)
    # run-length segments: consecutive bins, same chromosome, same segment value
    brk <- c(TRUE, f$chromosome[-1] != f$chromosome[-nrow(f)] | abs(diff(s)) > 1e-9)
    id <- cumsum(brk)
    sg <- data.frame(
        cell_id = cell,
        chromosome = tapply(f$chromosome, id, `[`, 1),
        start = tapply(f$start, id, min),
        end = tapply(f$end, id, max),
        n_bins = as.integer(tapply(id, id, length)),
        rel = round(as.numeric(tapply(rel_all, id, `[`, 1)), 4),
        stringsAsFactors = FALSE)
    sg$state <- ifelse(sg$rel > a$gain, "gain", ifelse(sg$rel < a$loss, "loss", "neutral"))
    keep <- (sg$end - sg$start + 1) / 1e6 >= a$min_seg_mb
    kept <- sg[keep, ]
    be <- NULL
    if (nrow(kept) > 1) {
        i <- which(kept$chromosome[-1] == kept$chromosome[-nrow(kept)] &
                   kept$state[-1] != kept$state[-nrow(kept)])
        if (length(i)) {
            lo <- kept$end[i] - bin_bp + 1; hi <- kept$start[i + 1] + bin_bp - 1
            be <- data.frame(cell_id = cell, chromosome = kept$chromosome[i],
                             pos = round((lo + hi) / 2), lo = lo, hi = hi,
                             left_state = kept$state[i], right_state = kept$state[i + 1],
                             left_rel = kept$rel[i], right_rel = kept$rel[i + 1])
        }
    }
    alt <- kept$state != "neutral"
    summ <- data.frame(cell_id = cell, n_segments = nrow(kept), n_breakends = if (is.null(be)) 0L else nrow(be),
                       n_gain_seg = sum(kept$state == "gain"), n_loss_seg = sum(kept$state == "loss"),
                       frac_altered = round(sum(kept$n_bins[alt]) / sum(kept$n_bins), 4),
                       altered_mb = round(sum((kept$end - kept$start + 1)[alt]) / 1e6, 2))
    list(seg = sg, be = be, cell = summ)
}

dirs <- list.dirs(a$relcn, recursive = FALSE)
dirs <- dirs[grepl(paste0("_", a$bin, "kb$"), basename(dirs))]
res <- list()
for (d in dirs) {
    obj <- readRDS(file.path(d, paste0("copyNumbersSegmented_", a$bin, "kb.rds")))
    fd <- fData(obj); sgm <- assayDataElement(obj, "segmented")
    for (c in colnames(obj)) res[[c]] <- call_cell(c, fd, sgm[, c])
}
dir.create(dirname(a$out), showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write.table(x, paste0(a$out, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(do.call(rbind, lapply(res, `[[`, "seg")), "_segments.tsv")
be <- do.call(rbind, lapply(res, `[[`, "be"))
w(if (is.null(be)) data.frame(cell_id = character()) else be, "_breakends.tsv")
w(do.call(rbind, lapply(res, `[[`, "cell")), "_cells.tsv")
cat(length(res), "cells,", if (is.null(be)) 0 else nrow(be), "breakends ->", a$out, "\n")

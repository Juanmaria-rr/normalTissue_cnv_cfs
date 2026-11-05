#!/usr/bin/env Rscript
# Count usable bins per bin size, replicating the filter that sWGS_relativeCN applies
# (QDNAseqmod applyFilters(residual=TRUE, blacklist=TRUE), default chromosomes
# excluded = X, Y, MT). fitting_absCN uses this count (nbins_ref_genome) in the
# 15-nrpbpcc depth threshold.
#
# 200 kb: annotation shipped in external/sWGS_relativeCN/bin.annotations;
# other sizes: QDNAseq.hg19 (same source sWGS_relativeCN uses).
#
# Usage: Rscript scripts/00_metadata/count_usable_bins.R data/metadata/usable_bins.tsv

suppressPackageStartupMessages({ library(QDNAseq); library(Biobase) })

args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args) >= 1) args[1] else "data/metadata/usable_bins.tsv"
ann_200 <- "external/sWGS_relativeCN/bin.annotations/BinAnnotations_QDNAseq_200kb.rds"

usable <- function(b) {
  d <- if (is(b, "AnnotatedDataFrame")) pData(b) else as.data.frame(b)
  keep <- !d$chromosome %in% c("X", "Y", "MT") & !is.na(d$gc) &
          !is.na(d$residual) & d$blacklist == 0
  c(total = nrow(d), autosomal = sum(!d$chromosome %in% c("X", "Y", "MT")),
    usable = sum(keep))
}

sizes <- c(200, 500, 1000)
res <- do.call(rbind, lapply(sizes, function(s) {
  b <- if (s == 200) readRDS(ann_200)
       else suppressMessages(getBinAnnotations(binSize = s, genome = "hg19"))
  u <- usable(b)
  data.frame(bin_kb = s, total_bins = u[["total"]], autosomal_bins = u[["autosomal"]],
             usable_bins = u[["usable"]],
             # 15 nrpbpcc x ploidy 2 x usable bins / 0.91 (fitting_absCN, purity 1)
             reads_needed_diploid = round(15 * 2 * u[["usable"]] / 0.91))
}))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res)

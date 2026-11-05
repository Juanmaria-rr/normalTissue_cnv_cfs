#!/usr/bin/env Rscript
# Export relative CN of selected cells from the copynumber workflow outputs, in the long
# format that build_curation_batch.py turns into a curation-app batch.
#
# Input:  <relcn>/<donor>_<bin>kb/copyNumbersSegmented_<bin>kb.rds (QDNAseqCopyNumbers)
# Output: <out>_bins.tsv.gz   cell_id, chromosome, start, end, use, copynumber, segmented
#         <out>_cells.tsv     cell_id, total_reads, used_bins
# copynumber/segmented are in read-count space (as sWGS_relativeCN leaves them); bins
# filtered out by applyFilters have use = FALSE and NA values.
#
# Usage (cnv-cfs/copynumber image):
#   Rscript scripts/04_curation/export_relcn.R --relcn <cn out>/relCN --bin 200 \
#       [--cells cells.txt] --out curation/staging/test20_200kb
#   --cells: one cell_id per line; default = every cell found.

suppressPackageStartupMessages({
    library(argparse)
    library(Biobase)
    library(QDNAseqmod)
})

p <- ArgumentParser()
p$add_argument("--relcn", required = TRUE)
p$add_argument("--bin", type = "integer", required = TRUE)
p$add_argument("--cells")
p$add_argument("--out", required = TRUE)
a <- p$parse_args()

dirs <- list.dirs(a$relcn, recursive = FALSE)
dirs <- dirs[grepl(paste0("_", a$bin, "kb$"), basename(dirs))]
want <- if (!is.null(a$cells)) readLines(a$cells) else NULL

bins <- list(); cells <- list()
for (d in dirs) {
    obj <- readRDS(file.path(d, paste0("copyNumbersSegmented_", a$bin, "kb.rds")))
    keep <- if (is.null(want)) colnames(obj) else intersect(colnames(obj), want)
    if (!length(keep)) next
    fd <- fData(obj)
    cn <- assayDataElement(obj, "copynumber"); sg <- assayDataElement(obj, "segmented")
    for (c in keep) {
        bins[[c]] <- data.frame(cell_id = c, chromosome = fd$chromosome, start = fd$start,
                                end = fd$end, use = fd$use, copynumber = cn[, c], segmented = sg[, c])
        cells[[c]] <- data.frame(cell_id = c, total_reads = pData(obj)[c, "total.reads"],
                                 used_bins = sum(fd$use))
    }
}
if (!is.null(want) && length(setdiff(want, names(cells))))
    stop("cells not found: ", paste(setdiff(want, names(cells)), collapse = ", "))

dir.create(dirname(a$out), showWarnings = FALSE, recursive = TRUE)
gz <- gzfile(paste0(a$out, "_bins.tsv.gz"), "w")
write.table(do.call(rbind, bins), gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
write.table(do.call(rbind, cells), paste0(a$out, "_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(length(cells), "cells ->", a$out, "\n")

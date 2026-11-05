#!/usr/bin/env Rscript
# Data and figures for the QC levels 1-2 method document, from the 20-cell test.
#
# Writes to --out (docs/qc_levels_1_2):
#   read_fate.tsv        per cell: reads split into unmapped / low MAPQ (<37) / duplicate / usable
#   dup_vs_paper.tsv     per cell duplicate %, with the donor's duplicate rate in Suppl. Table 2
#   noise_vs_depth.tsv   per cell x bin: nrpbpcc, MAPD (classic) and noise_excess (Poisson-normalised)
#   clonality_curves.tsv ploidy vs clonality (best purity per ploidy) for two example cells, 200 kb
#   rt_deciles.tsv       per cell: mean relative CN by replication-timing decile, 200 kb
#   fig_*.png            genome-wide CN profiles of example cells
#
# Usage (cnv-cfs/copynumber image, which has samtools 1.9):
#   Rscript scripts/docs/qc_doc_data.R --bam-dir <align pubdir>/bam --relcn <cn out>/relCN \
#       --abscn <cn out>/absCN --align-qc data/pilot/align_qc_test20.tsv \
#       --paper-qc data/metadata/paper_scdna_qc.tsv \
#       --rt external/sWGS_relativeCN/bin.annotations/bins_with_replicationTime_50kb.tsv \
#       --out docs/qc_levels_1_2

suppressPackageStartupMessages({
    library(argparse)
    library(Biobase)
    library(QDNAseqmod)
    library(ggplot2)
})

p <- ArgumentParser()
for (a in c("--bam-dir", "--relcn", "--abscn", "--align-qc", "--paper-qc", "--rt", "--out"))
    p$add_argument(a, required = TRUE)
a <- p$parse_args()
dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
wr <- function(x, f) write.table(x, file.path(a$out, f), sep = "\t", quote = FALSE, row.names = FALSE)

aqc <- read.delim(a$align_qc)

# ---- read fate (level 1): samtools counts per BAM ------------------------------------------
cnt <- function(bam, args) as.numeric(system(paste("samtools view -c", args, shQuote(bam)), intern = TRUE))
fate <- do.call(rbind, lapply(aqc$cell_id, function(c) {
    bam <- file.path(a$bam_dir, paste0(c, ".bam"))
    total   <- cnt(bam, "-F 0x900")           # primary records
    mapped  <- cnt(bam, "-F 0x904")           # primary mapped
    dup     <- cnt(bam, "-F 0x904 -f 0x400")  # primary mapped duplicates
    usable  <- cnt(bam, "-q 37 -F 0xD04")     # primary mapped, MAPQ >= 37, not duplicate
    data.frame(cell_id = c, donor = sub("_.*", "", c), total = total,
               unmapped = total - mapped, duplicate = dup,
               low_mapq = mapped - dup - usable, usable = usable)
}))
wr(fate, "read_fate.tsv")

# ---- duplicates vs paper -----------------------------------------------------------------
pq <- read.delim(a$paper_qc)
dup <- merge(aqc[, c("cell_id", "donor", "coassay", "dup_pct")],
             data.frame(donor = pq$donor, paper_dup_pct = pq$duplication_rate), by = "donor")
wr(dup[order(dup$donor, dup$cell_id), ], "dup_vs_paper.tsv")

# ---- per-bin helpers ---------------------------------------------------------------------
rt <- read.delim(a$rt)
rt_for <- function(bin_kb) {
    bp <- bin_kb * 1000
    key <- paste0(rt$seqnames, ":", floor((rt$start - 1) / bp) * bp + 1)
    tapply(rt$replicationTime, key, mean, na.rm = TRUE)
}
load_bin <- function(bin_kb) {
    dirs <- list.dirs(a$abscn, recursive = FALSE)
    dirs <- dirs[grepl(paste0("_", bin_kb, "kb$"), basename(dirs))]
    lapply(dirs, function(d) list(
        rel = readRDS(file.path(a$relcn, basename(d), paste0("copyNumbersSegmented_", bin_kb, "kb.rds"))),
        abs = readRDS(file.path(d, "absoluteCN.rds")), dir = d))
}
auto_use <- function(obj) { fd <- fData(obj); fd$use & fd$chromosome %in% as.character(1:22) }

# ---- noise vs depth ----------------------------------------------------------------------
noise <- do.call(rbind, lapply(c(200, 500, 1000), function(b) do.call(rbind, lapply(load_bin(b), function(o) {
    k <- auto_use(o$rel); chr <- fData(o$rel)$chromosome[k]
    same <- chr[-1] == chr[-length(chr)]
    do.call(rbind, lapply(colnames(o$abs), function(c) {
        x <- assayDataElement(o$rel, "copynumber")[k, c]; s <- assayDataElement(o$rel, "segmented")[k, c]
        lr <- log2(pmax(x / median(x, na.rm = TRUE), 1e-3))
        z <- (x[-1] - x[-length(x)]) / sqrt(pmax(s[-1] + s[-length(s)], 1e-9))
        data.frame(cell_id = c, donor = sub("_.*", "", c), bin_kb = b,
                   nrpbpcc = round(mean(x, na.rm = TRUE) / pData(o$abs)[c, "ploidy"], 2),
                   mapd = round(median(abs(diff(lr)[same]), na.rm = TRUE), 4),
                   noise_excess = round(mad(z[same], na.rm = TRUE), 3))
    }))
}))))
wr(noise, "noise_vs_depth.tsv")

# ---- clonality curves --------------------------------------------------------------------
curves <- do.call(rbind, lapply(c("P36_chip1_co_r51c33", "P36_chip1_co_r65c23"), function(c) {
    r <- read.delim(file.path(a$abscn, "P36_200kb", paste0(c, "_resfitting.txt")))
    r <- r[order(r$ploidy), c("ploidy", "clonality", "powered")]
    data.frame(cell_id = c, r)
}))
wr(curves, "clonality_curves.tsv")

# ---- replication-timing deciles (200 kb) -------------------------------------------------
rtb <- rt_for(200)
rtd <- do.call(rbind, lapply(load_bin(200), function(o) {
    k <- auto_use(o$rel); fd <- fData(o$rel)[k, ]
    v <- rtb[paste0(fd$chromosome, ":", fd$start)]
    dec <- cut(v, quantile(v, 0:10 / 10, na.rm = TRUE), include.lowest = TRUE, labels = 1:10)
    do.call(rbind, lapply(colnames(o$abs), function(c) {
        x <- assayDataElement(o$rel, "copynumber")[k, c]; rel <- x / median(x, na.rm = TRUE)
        m <- tapply(rel, dec, mean, na.rm = TRUE)
        data.frame(cell_id = c, rt_decile = as.integer(names(m)), mean_rel_cn = round(as.numeric(m), 4))
    }))
}))
wr(rtd, "rt_deciles.tsv")

# ---- profile figures ---------------------------------------------------------------------
profile_df <- function(bin_kb, donor, cell) {
    d <- file.path(a$abscn, paste0(donor, "_", bin_kb, "kb"))
    o <- readRDS(file.path(d, "absoluteCN.rds")); k <- auto_use(o); fd <- fData(o)[k, ]
    chr <- as.integer(fd$chromosome)
    len <- tapply(fd$end, chr, max); off <- c(0, cumsum(as.numeric(len)))[seq_along(len)]; names(off) <- names(len)
    data.frame(chr = chr, start_mb = fd$start / 1e6, pos = (fd$start + off[as.character(chr)]) / 1e6,
               cn = assayDataElement(o, "copynumber")[k, cell],
               seg = assayDataElement(o, "segmented")[k, cell],
               bin = paste(bin_kb, "kb"), ploidy = pData(o)[cell, "ploidy"])
}
theme_doc <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"))
chr_lines <- function(df) {
    b <- aggregate(pos ~ chr, df, max)
    list(geom_vline(xintercept = b$pos, colour = "grey85", linewidth = 0.3))
}
plot_profiles <- function(dfs, title, subtitle, file, h = 6.5, ymax = 6) {
    df <- do.call(rbind, dfs)
    df$panel <- factor(df$panel, levels = unique(df$panel))
    g <- ggplot(df, aes(pos, pmin(cn, ymax))) +
        geom_vline(data = aggregate(pos ~ chr + panel, df, max), aes(xintercept = pos),
                   colour = "grey88", size = 0.3) +
        geom_point(size = 0.25, colour = "grey55", alpha = 0.6) +
        geom_point(aes(y = pmin(seg, ymax)), size = 0.35, colour = "#2563eb") +
        facet_wrap(~panel, ncol = 1) +
        scale_y_continuous(breaks = 0:ymax, limits = c(0, ymax)) +
        labs(x = "genomic position, autosomes 1-22 (Mb)", y = "absolute copy number",
             title = title, subtitle = subtitle) + theme_doc
    ggsave(file.path(a$out, file), g, width = 10, height = h, dpi = 150, bg = "white")
}
lab <- function(df, txt) { df$panel <- txt; df }

# 1. same cell at three bin sizes (focal chr7 gain)
c1 <- "P08_chip1_co_r71c48"
plot_profiles(lapply(c(200, 500, 1000), function(b) lab(profile_df(b, "P08", c1),
              sprintf("%s kb", b))), paste(c1, "- 1.2 M usable reads"),
              "grey = bin CN, blue = segment CN; the chr6 q-arm gain (68-95 Mb, CN 3) is called at all three bin sizes",
              "fig_profile_bins_P08_r71c48.png", h = 7)
# 2. depth: deep vs shallow diploid at 200 kb
plot_profiles(list(lab(profile_df(200, "P18", "P18_chip1_dna_r55c17"), "P18_chip1_dna_r55c17 - 602 k usable reads, powered"),
                   lab(profile_df(200, "P20", "P20_chip1_co_r40c13"), "P20_chip1_co_r40c13 - 33 k usable reads, not powered")),
              "Depth sets the scatter: same 200 kb bins, 18x fewer reads",
              "grey = bin CN, blue = segment CN at the chosen ploidy", "fig_profile_depth_200kb.png", h = 5)
# 3. bad fit: P36_r51c33 at 200 kb vs 500 kb
c3 <- "P36_chip1_co_r51c33"
plot_profiles(list(lab(profile_df(200, "P36", c3), "200 kb - fitted ploidy 1.6, clonality 0.18"),
                   lab(profile_df(500, "P36", c3), "500 kb - fitted ploidy 2.0")),
              paste(c3, "- the 200 kb fit is wrong"),
              "A ploidy that leaves segments between integers is flagged by high clonality",
              "fig_profile_badfit_P36_r51c33.png", h = 5)
# 4. single-bin blip: P07_co_r41c34 at 200 kb, zoom on the altered bin's chromosome
c4 <- "P07_chip1_co_r41c34"; d4 <- profile_df(200, "P07", c4)
chg <- which(round(d4$seg) != 2)
if (length(chg)) {
    zc <- d4$chr[chg[1]]; z <- d4[d4$chr == zc, ]; alt <- d4[chg, ]
    g <- ggplot(z, aes(start_mb, pmin(cn, 10))) + geom_point(size = 0.6, colour = "grey55") +
        geom_point(aes(y = pmin(seg, 10)), size = 0.8, colour = "#2563eb") +
        scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
        labs(x = sprintf("chr%s position (Mb)", zc), y = "absolute copy number",
             title = sprintf("%s, chr%s at 200 kb: a %d-bin segment at CN %d = 2 breakpoints", c4, zc,
                             nrow(alt), round(alt$seg[1])),
             subtitle = sprintf("chr%s:%.1f-%.1f Mb; short isolated segments inflate the breakpoint count",
                                zc, min(alt$start_mb), max(alt$start_mb) + 0.2)) + theme_doc
    ggsave(file.path(a$out, "fig_blip_P07_co_r41c34.png"), g, width = 10, height = 3.6, dpi = 150, bg = "white")
}
cat("done ->", a$out, "\n")

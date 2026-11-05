# Copy-number image for CNV_cfs: relative CN (sWGS_relativeCN, QDNAseqmod) and candidate
# absolute-CN fits (fitting_absCN). Built from the group's own conda environment file
# (R 4.0.2, QDNAseq 1.26.0, QDNAseq.hg19 1.20.0) plus the packages their
# Check_R_packages.R installs separately (argparse, this.path, Ckmeans.1d.dp) and QDNAseqmod
# from the bundled tarball. The pipelines are copied from the pinned submodules.
# Build context: repository root.
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="cnv-cfs/copynumber" \
      org.opencontainers.image.source="https://github.com/Juanmaria-rr/CNV_cfs"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps fonts-dejavu-core fontconfig \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER external/sWGS_relativeCN/config_env.yml /tmp/swgs_env.yml
RUN micromamba create -y -n swgs -f /tmp/swgs_env.yml && micromamba clean -a -y

ENV PATH=/opt/conda/envs/swgs/bin:$PATH

# Check_R_packages.R installs these from CRAN unpinned. Installing them with conda could
# re-solve the pinned R 4.0.2 environment, so take them from a CRAN snapshot dated like the
# environment (Feb 2022) and compile with the env's own toolchain.
RUN Rscript -e 'install.packages(c("argparse", "this.path", "Ckmeans.1d.dp"), repos = "https://packagemanager.posit.co/cran/2022-02-15")'

COPY --chown=$MAMBA_USER:$MAMBA_USER external/sWGS_relativeCN /opt/sWGS_relativeCN
COPY --chown=$MAMBA_USER:$MAMBA_USER external/fitting_absCN /opt/fitting_absCN

# QDNAseqmod is not on conda; install it from the tarball shipped with sWGS_relativeCN.
RUN Rscript -e 'install.packages("/opt/sWGS_relativeCN/QDNAseqmod_1.21.0.tar.gz", repos = NULL, type = "source")' \
 && Rscript -e 'for (p in c("QDNAseqmod","QDNAseq","QDNAseq.hg19","Biobase","argparse","this.path","Ckmeans.1d.dp","splitstackshape","doMC","foreach","tidyverse","matrixStats","GenomicRanges")) suppressPackageStartupMessages(library(p, character.only = TRUE)); cat("R packages OK\n")'

# Upstream bug: pureSamples_absCN.R and tumorSamples_absCN.R source
# "<script dir>/functions_absCN.R", but the file lives in utils/. Link it instead of
# editing the pinned submodule.
RUN ln -s /opt/fitting_absCN/utils/functions_absCN.R /opt/fitting_absCN/functions_absCN.R \
 && micromamba env export -n swgs > /opt/conda/env_export.yml

USER root

# CNV_cfs patch (pureSamples_absCN.safe.patch): per-cell error handling in the fitting loop.
# A cell that errors is written to failed_cells.tsv (with the message) instead of aborting
# the whole donor; cells skipped by --total.reads are listed there too. The pinned upstream
# script is left untouched; the patched copy is pureSamples_absCN_safe.R (upstream ships
# CRLF line endings; the copy is normalised to LF, which R does not care about).
COPY containers/copynumber/pureSamples_absCN.safe.patch /tmp/pureSamples_absCN.safe.patch
RUN apt-get update && apt-get install -y --no-install-recommends patch && rm -rf /var/lib/apt/lists/* \
 && sed 's/\r$//' /opt/fitting_absCN/pureSamples_absCN.R > /opt/fitting_absCN/pureSamples_absCN_safe.R \
 && patch /opt/fitting_absCN/pureSamples_absCN_safe.R /tmp/pureSamples_absCN.safe.patch

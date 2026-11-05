# Alignment image for CNV_cfs: SRA -> FASTQ (R1) -> trimmed -> BAM (hg19).
# Tool versions match macintyrelab/scDNAseq config_env.yml (bwa 0.7.17, samtools 1.15,
# cutadapt 3.7, trim-galore 0.6.7, fastqc 0.11.9). sra-tools is not in that environment and
# is added here. HealthOmics runs linux/amd64.
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="cnv-cfs/align" \
      org.opencontainers.image.source="https://github.com/Juanmaria-rr/CNV_cfs"

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.9 \
        sra-tools=3.1.1 \
        bwa=0.7.17 \
        samtools=1.15 \
        cutadapt=3.7 \
        trim-galore=0.6.7 \
        fastqc=0.11.9 \
        pigz \
        procps-ng \
    && micromamba env export -n base > /opt/conda/env_export.yml \
    && micromamba clean -a -y

# Nextflow runs task scripts with bash and needs the env on PATH without activation.
ENV PATH=/opt/conda/bin:$PATH
USER root

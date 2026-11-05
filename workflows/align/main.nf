// Per-cell alignment of PRJNA1003661 scDNA-seq runs staged in S3 (.sra files).
// Mirrors macintyrelab/scDNAseq steps 1-6 with the settings chosen for this dataset:
//   1. fasterq-dump each run; keep read 1 only (single-end runs: the only read), so
//      paired and single-end runs give comparable counts for QDNAseq.
//   2. concatenate all runs of a cell (top-up sequencing of the same library).
//   3. no pre-amplification hard trim (no constant 5' prefix in ACT/Tn5 reads);
//      cutadapt removes Tn5 mosaic-end read-through, dropping reads < 20 bp.
//   4. bwa mem -> samtools sort -> samtools markdup (duplicates flagged, not removed;
//      QDNAseq excludes them when counting) -> index, flagstat.
// Cells are processed in chunks so that each task stages the bwa index once.
nextflow.enable.dsl = 2

params.samplesheet = null   // TSV with header: cell_id <TAB> runs (comma-separated SRR)
params.sra_prefix  = null   // S3 prefix holding <SRR>/<SRR>
params.ref         = null   // S3 URI of genome.fa; .fai and bwa index files alongside
params.chunk_size  = 50
params.min_length  = 20
params.adapter     = 'CTGTCTCTTATACACATCT'   // Tn5 mosaic end
params.outdir      = '/mnt/workflow/pubdir'   // HealthOmics output location

process ALIGN_CHUNK {
    container "${params.container_align}"
    cpus 8
    memory '32 GB'
    errorStrategy 'retry'
    maxRetries 1
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(chunk_id), val(cells_tsv), path(sra_files)
    path ref_files

    output:
    path "bam/*.bam*", emit: bam
    path "qc/*",       emit: qc

    script:
    """
    mkdir -p bam qc tmp
    # Tool versions (outside pipefail: bwa without arguments exits 1)
    { bwa 2>&1 | grep -i '^Version'; samtools --version | head -1;
      echo "cutadapt \$(cutadapt --version)"; echo "fasterq-dump \$(fasterq-dump --version | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -1)"; } > qc/versions_chunk${chunk_id}.txt || true
    set -euo pipefail
    printf '%s\\n' "${cells_tsv}" > cells.tsv

    while IFS=\$'\\t' read -r cell runs; do
        : > tmp/\${cell}.fq
        for r in \${runs//,/ }; do
            fasterq-dump --split-3 -e ${task.cpus} -t tmp -O tmp/\$r ./\$r
            if [[ -f tmp/\$r/\${r}_1.fastq ]]; then cat tmp/\$r/\${r}_1.fastq; else cat tmp/\$r/\$r.fastq; fi >> tmp/\${cell}.fq
            rm -rf tmp/\$r
        done
        cutadapt -j ${task.cpus} -a ${params.adapter} -m ${params.min_length} \\
            -o tmp/\${cell}.trim.fq tmp/\${cell}.fq > qc/\${cell}.cutadapt.txt
        bwa mem -t ${task.cpus} -R "@RG\\tID:\${cell}\\tSM:\${cell}" genome.fa tmp/\${cell}.trim.fq 2> qc/\${cell}.bwa.log \\
            | samtools sort -@ 2 -o tmp/\${cell}.sorted.bam -
        samtools markdup -s tmp/\${cell}.sorted.bam bam/\${cell}.bam 2> qc/\${cell}.markdup.txt
        samtools index bam/\${cell}.bam
        samtools flagstat bam/\${cell}.bam > qc/\${cell}.flagstat.txt
        rm -f tmp/\${cell}.*
    done < cells.tsv
    """
}

workflow {
    ['samplesheet', 'sra_prefix', 'ref', 'container_align'].each {
        if (!params[it]) error "params.${it} is required"
    }
    // HealthOmics only accepts S3 prefixes (not objects) with a trailing '/'
    def sra_prefix = params.sra_prefix.toString().replaceAll('/+$', '')
    ref_files = ['', '.fai', '.amb', '.ann', '.bwt', '.pac', '.sa'].collect { file("${params.ref}${it}") }

    chunks = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .buffer(size: params.chunk_size as int, remainder: true)
        .toList()
        .flatMap { cs -> cs.withIndex().collect { c, i -> [i, c] } }
        .map { i, c ->
            def runs = c.collectMany { it.runs.tokenize(',') }
            tuple(i,
                  c.collect { "${it.cell_id}\t${it.runs}" }.join('\n'),
                  runs.collect { file("${sra_prefix}/${it}/${it}") })
        }

    ALIGN_CHUNK(chunks, ref_files)
}

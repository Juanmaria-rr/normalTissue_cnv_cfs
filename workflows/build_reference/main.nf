// Build the bwa index (and samtools faidx) for the hg19 FASTA used by CNV_cfs.
// One-off HealthOmics workflow; outputs land in the run's output URI (pubdir/hg19/).
nextflow.enable.dsl = 2

params.fasta = null   // S3 URI of genome.fa

process BWA_INDEX {
    container "${params.container_align}"
    cpus 2
    memory '16 GB'
    publishDir '/mnt/workflow/pubdir/hg19', mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}*", emit: index
    path "versions.txt"

    script:
    """
    samtools faidx ${fasta}
    bwa index -a bwtsw ${fasta}
    { bwa 2>&1 | grep -i '^Version'; samtools --version | head -1; } > versions.txt
    """
}

workflow {
    if (!params.fasta) error "params.fasta (S3 URI of genome.fa) is required"
    BWA_INDEX(file(params.fasta))
}

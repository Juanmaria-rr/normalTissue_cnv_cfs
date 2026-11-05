// Relative CN (sWGS_relativeCN) and candidate absolute-CN fits (fitting_absCN) for per-cell
// BAMs produced by workflows/align. Pinned upstream scripts are run unmodified from the
// cnv-cfs/copynumber image (/opt/sWGS_relativeCN, /opt/fitting_absCN).
//
// One task per (donor, bin size). Grouping by donor is only sharding: QDNAseq corrections,
// filters and segmentation are per sample, so results do not depend on which cells share
// a task.
//
// Note on --filters powered: pureSamples_absCN.R does not drop cells that fail the depth
// filter; it falls back to the best unfiltered fit. Depth filtering is done downstream
// (cn_qc_metrics.R, diploid threshold).
//
// Fitting runs pureSamples_absCN_safe.R (image patch): a cell that errors is listed in
// failed_cells.tsv instead of aborting the donor, as are cells below --total.reads.
// In single-cell mode the script only copies each BAM to a scratch folder that is deleted,
// so ABS_FITS gets empty placeholder BAMs instead of staging the real ones.
//
// Resume: --reuse_relcn lists donor_bin tags (e.g. P07_200) whose relative CN is taken from
// --relcn_prefix (a previous run's pubdir/relCN/) instead of being recomputed.
nextflow.enable.dsl = 2

params.samplesheet  = null      // TSV with header cell_id <TAB> ...; donor = cell_id prefix before '_'
params.bam_prefix   = null      // S3 prefix holding <cell_id>.bam and .bam.bai (align output bam/)
params.bin_sizes    = '200,500,1000'
params.undo_sd      = 1
params.cnaware      = 'NO'
params.mappability  = 'YES'
params.reptiming    = 'NO'      // only applied upstream at 30/50 kb
params.downsampling = 15        // nrpbpcc for the `powered` flag
params.filters      = 'powered'
params.total_reads  = 20000     // cells at or below are not fitted (listed in failed_cells.tsv); far
                                // below the 1 Mb diploid threshold (64k), so no usable cell is lost
params.relcn_prefix = null      // S3 prefix of a previous run's relCN/ (trailing /)
params.reuse_relcn  = ''        // comma-separated donor_bin tags to take from relcn_prefix
params.outdir       = '/mnt/workflow/pubdir'

process RELATIVE_CN {
    tag "${donor}_${bin}kb"
    container "${params.container_copynumber}"
    cpus 4
    memory '32 GB'
    errorStrategy 'retry'       // one retry for infrastructure failures, then the run fails
    maxRetries 1
    publishDir "${params.outdir}/relCN/${donor}_${bin}kb", mode: 'copy', saveAs: { it.replaceFirst(/^out\//, '') }

    input:
    tuple val(donor), path(bams), val(bin)

    output:
    tuple val(donor), val(bin), path("out/*.rds"), emit: rds
    path "out/*"

    script:
    """
    set -euo pipefail
    mkdir -p bams && mv *.bam *.bai bams/
    # swgs-relativeCN.R sources plotting/visualize_cn.R relative to the working directory
    WD=\$PWD
    cd /opt/sWGS_relativeCN
    Rscript swgs-relativeCN.R --bin.size ${bin} --genome hg19 --undo.SD ${params.undo_sd} \\
        --cnaware ${params.cnaware} --mappability ${params.mappability} --reptiming ${params.reptiming} \\
        --bams.dir \$WD/bams --out.dir \$WD/out
    cd \$WD
    Rscript -e 'cat(R.version.string, "\\n"); for (p in c("QDNAseqmod","QDNAseq","QDNAseq.hg19")) cat(p, as.character(packageVersion(p)), "\\n")' > out/versions.txt
    """
}

process ABS_FITS {
    tag "${donor}_${bin}kb"
    container "${params.container_copynumber}"
    cpus 10                     // pureSamples_absCN.R calls registerDoMC(10)
    memory '32 GB'
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.outdir}/absCN/${donor}_${bin}kb", mode: 'copy', saveAs: { it.replaceFirst(/^out\//, '') }

    input:
    tuple val(donor), val(bin), path(rds), val(cells)

    output:
    path "out/*"

    script:
    """
    set -euo pipefail
    mkdir -p bams reads && mv *.rds reads/
    for c in ${cells.join(' ')}; do : > bams/\$c.bam; done
    Rscript /opt/fitting_absCN/pureSamples_absCN_safe.R --bin.size ${bin} --genome hg19 \\
        --single.cells YES --downsampling ${params.downsampling} --downs.method original \\
        --filters ${params.filters} --total.reads ${params.total_reads} \\
        --bams.dir \$PWD/bams --reads.dir \$PWD/reads --out.dir \$PWD/out
    # single-cell mode copies every BAM into bams/downsampled_bams; drop the copies
    rm -rf bams/downsampled_bams
    """
}

workflow {
    ['samplesheet', 'bam_prefix', 'container_copynumber'].each {
        if (!params[it]) error "params.${it} is required"
    }
    // HealthOmics only accepts S3 prefixes (not objects) with a trailing '/'
    def bam_prefix = params.bam_prefix.toString().replaceAll('/+$', '')
    bins = Channel.fromList(params.bin_sizes.toString().tokenize(',').collect { it.trim() as int })

    by_donor = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.cell_id.tokenize('_')[0],
                          [file("${bam_prefix}/${r.cell_id}.bam"),
                           file("${bam_prefix}/${r.cell_id}.bam.bai")]) }
        .groupTuple()
        .map { donor, files -> tuple(donor, files.flatten()) }

    def reuse = params.reuse_relcn.toString().tokenize(',')*.trim()
    if (reuse && !params.relcn_prefix) error "params.relcn_prefix is required with reuse_relcn"
    def relcn_prefix = params.relcn_prefix ? params.relcn_prefix.toString().replaceAll('/+$', '') : ''
    combos = by_donor.combine(bins)     // (donor, files, bin)
    reused = combos.filter { d, f, b -> "${d}_${b}".toString() in reuse }
        .map { d, f, b -> tuple(d, b, [file("${relcn_prefix}/${d}_${b}kb/copyNumbersSegmented_${b}kb.rds")]) }
    RELATIVE_CN(combos.filter { d, f, b -> !("${d}_${b}".toString() in reuse) })

    cells = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.cell_id.tokenize('_')[0], r.cell_id) }
        .groupTuple()
    ABS_FITS(RELATIVE_CN.out.rds.mix(reused).combine(cells, by: 0))
}

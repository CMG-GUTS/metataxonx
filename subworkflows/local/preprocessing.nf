/*

    READ TRIM + QC + MERGE

*/
include { CUTADAPT } from '../../modules/nf-core/cutadapt.nf'
include { SEQKIT_SAMPLE } from '../../modules/nf-core/seqkit/sample.nf'
include { PEAR } from '../../modules/nf-core/pear.nf'

include { FASTQC as FASTQC_reads } from '../../modules/nf-core/fastqc.nf'
include { FASTQC as FASTQC_trim } from '../../modules/nf-core/fastqc.nf'
include { FASTQC as FASTQC_pear } from '../../modules/nf-core/fastqc.nf'

include { MULTIQC as MULTIQC_reads } from '../../modules/nf-core/multiqc.nf'
include { MULTIQC as MULTIQC_trim } from '../../modules/nf-core/multiqc.nf'
include { MULTIQC as MULTIQC_pear } from '../../modules/nf-core/multiqc.nf'

include { MERGE_MULTIQC_STATS } from '../../modules/local/merge_multiqc_stats.nf'

workflow PREPROCESSING {
    take:
    reads

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions = Channel.empty()    

    FASTQC_reads(reads)
    MULTIQC_reads(
        FASTQC_reads.out.zip.collect{ it[1] },
        "raw",
        [], [], [], [], []
    )
    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_reads.out.report)
    ch_versions = ch_versions.mix(FASTQC_reads.out.versions)
    ch_versions = ch_versions.mix(MULTIQC_reads.out.versions)

    if (!params.bypass_trim) {
        CUTADAPT(
            reads
        ).trimmed.set { ch_trimmed_reads }

        FASTQC_trim(ch_trimmed_reads)
        MULTIQC_trim(
            FASTQC_trim.out.zip.collect{ it[1] },
            "trimmed",
            [], [], [], [], []
        )
        ch_multiqc_stats_trim = MULTIQC_trim.out.multiqc_stats
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_trim.out.report)
    } else {
        ch_trimmed_reads = reads
        ch_multiqc_stats_trim = Channel.empty()
    }

    PEAR(
        ch_trimmed_reads
    ).assembled.set { ch_merged_reads }
    ch_versions = ch_versions.mix(PEAR.out.versions)

    FASTQC_pear(ch_merged_reads)
    MULTIQC_pear(
        FASTQC_pear.out.zip.collect{ it[1] },
        "decon",
        [], [], [], [], []
    )
    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_pear.out.report)

    MERGE_MULTIQC_STATS(
        MULTIQC_reads.out.multiqc_stats,
        ch_multiqc_stats_trim,
        MULTIQC_pear.out.multiqc_stats
    )

    SEQKIT_SAMPLE(
        ch_merged_reads
    ).reads.set { ch_subsampled_reads }
    
    emit:
    untrimmed           = reads
    trimmed             = ch_trimmed_reads
    merged              = ch_merged_reads
    subsampled_reads    = ch_subsampled_reads
    multiqc_report      = ch_multiqc_files
    read_stats          = MERGE_MULTIQC_STATS.out.read_stats
    versions            = ch_versions
}
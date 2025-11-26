/*

    READ TRIM + QC + MERGE

*/
include { CUTADAPT } from '../../modules/nf-core/cutadapt.nf'
include { SEQKIT_SAMPLE } from '../../modules/nf-core/seqkit/sample.nf'
include { PEAR } from '../../modules/nf-core/pear.nf'
include { FLASH } from '../../modules/nf-core/flash.nf'

include { FASTQC as FASTQC_reads } from '../../modules/nf-core/fastqc.nf'
include { FASTQC as FASTQC_trim } from '../../modules/nf-core/fastqc.nf'
include { FASTQC as FASTQC_merged } from '../../modules/nf-core/fastqc.nf'

workflow PREPROCESSING {
    take:
    reads

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions = Channel.empty()    

    FASTQC_reads(reads, "raw")

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_reads.out.zip.collect{ it[1] })
    ch_versions = ch_versions.mix(FASTQC_reads.out.versions)

    if (!params.bypass_trim) {
        CUTADAPT(
            reads
        ).reads.set { ch_trimmed_reads }
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{ it[1] })

        FASTQC_trim(ch_trimmed_reads, "trim")
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_trim.out.zip.collect{ it[1] })

    } else {
        ch_trimmed_reads = reads
    }

    if (!params.singleEnd) {
        if (params.read_merger_tool == "pear") {
            PEAR(
                ch_trimmed_reads
            ).assembled.set { ch_merged_reads }
            ch_versions = ch_versions.mix(PEAR.out.versions)
        }
        
        if (params.read_merger_tool == "flash") {
            FLASH(
                ch_trimmed_reads
            ).merged.set { ch_merged_reads }
            ch_versions = ch_versions.mix(FLASH.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(FLASH.out.histogram.collect{ it[1] })
        }

        FASTQC_merged(ch_merged_reads, params.read_merger_tool)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_merged.out.zip.collect{ it[1] })

    } else {
        ch_merged_reads = ch_trimmed_reads
    }
    
    SEQKIT_SAMPLE(
        ch_merged_reads
    ).reads.set { ch_subsampled_reads }
    
    emit:
    untrimmed           = reads
    trimmed             = ch_trimmed_reads
    merged              = ch_merged_reads
    subsampled          = ch_subsampled_reads
    multiqc_files       = ch_multiqc_files
    versions            = ch_versions
}
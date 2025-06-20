/*

    DENOISE VIA DADA2 (BATCH-MODE OPTIONAL)

*/
include { DADA2 } from '../../modules/local/dada2.nf'
include { CREATE_MAPPING } from '../../modules/local/create_mapping.nf'

workflow DENOISE {
    take:
    merged_reads

    main:
    ch_versions = Channel.empty()    

    CREATE_MAPPING(
        merged_reads
    ).dada2_mapping_file.set { ch_mapping }   

    DADA2(
        ch_mapping
    )
    ch_versions = ch_versions.mix(DADA2.out.versions)

    emit:
    dada2_mapping_file      = ch_mapping
    counts_table            = DADA2.out.seq_table
    denoise_stats           = DADA2.out.denoising_stats
    sequences               = DADA2.out.rep_seqs_fasta
    dada2_errors            = DADA2.out.dada2_errors
    versions                = ch_versions
}
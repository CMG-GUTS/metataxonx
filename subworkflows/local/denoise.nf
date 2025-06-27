/*

    DENOISE VIA DADA2 (BATCH-MODE OPTIONAL)

*/
include { DADA2 } from '../../modules/local/dada2.nf'
include { CREATE_MAPPING } from '../../modules/local/create_mapping.nf'

workflow DENOISE {
    take:
    meta_ch

    main:
    ch_versions = Channel.empty()   

    // Create Qiime2 mapping file
    mapping_rows_ch = meta_ch.map { meta, reads -> [ meta.id, reads ] }
    CREATE_MAPPING(
        mapping_rows_ch.collect(flat:false)
    ).mapping.set{ mapping_ch }

    DADA2(mapping_ch)
    ch_versions = ch_versions.mix(DADA2.out.versions)

    emit:
    dada2_mapping_file      = mapping_ch
    counts_table            = DADA2.out.seq_table
    denoise_stats           = DADA2.out.denoising_stats
    sequences               = DADA2.out.rep_seqs_fasta
    dada2_errors            = DADA2.out.dada2_errors
    versions                = ch_versions
}
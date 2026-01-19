/*

    ADDING TAXONOMY TO BIOM

*/
include { QIIME_IMPORT_EXPORT }     from    '../../modules/local/qiime_import_export.nf'
include { ASSIGN_TAXONOMY }         from    '../../modules/local/assign_taxonomy.nf'
include { COMBINE_TAXONOMY }        from    '../../modules/local/combine_taxonomy.nf'

workflow TAXONOMY {
    take:
    dada2_counts_table
    dada2_sequences
    classifier

    main:
    ch_versions = Channel.empty()    

    QIIME_IMPORT_EXPORT(
        dada2_counts_table,
        dada2_sequences
    )
    ch_versions = ch_versions.mix(QIIME_IMPORT_EXPORT.out.versions)

    ASSIGN_TAXONOMY(
        classifier,
        QIIME_IMPORT_EXPORT.out.repseqs_qza
    )
    ch_versions = ch_versions.mix(ASSIGN_TAXONOMY.out.versions)

    COMBINE_TAXONOMY(
        QIIME_IMPORT_EXPORT.out.biom_no_taxonomy,
        ASSIGN_TAXONOMY.out.taxonomy_tsv
    )
    ch_versions = ch_versions.mix(COMBINE_TAXONOMY.out.versions)

    emit:
    biom_without_taxonomy           = QIIME_IMPORT_EXPORT.out.biom_no_taxonomy
    biom_with_taxonomy              = COMBINE_TAXONOMY.out.biom_taxonomy
    taxonomy_file_tsv               = ASSIGN_TAXONOMY.out.taxonomy_tsv
    taxonomy_file_qza               = ASSIGN_TAXONOMY.out.taxonomy_qza
    qiime_asv_table                 = QIIME_IMPORT_EXPORT.out.asv_qza
    qiime_sequences_qza             = QIIME_IMPORT_EXPORT.out.repseqs_qza
    versions                        = ch_versions
}
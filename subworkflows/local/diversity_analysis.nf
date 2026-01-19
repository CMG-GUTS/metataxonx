/*

    COMPUTATION OF BETA & ALPHA DIVERSITY

*/
include { FASTTREE }                    from    '../../modules/local/fasttree.nf'
include { MINMAX }                      from    '../../modules/local/minmax.nf'
include { ALPHA_RAREFACTION }           from    '../../modules/local/alpha_rarefaction.nf'
include { CORE_DIVERSITY }              from    '../../modules/local/core_diversity.nf'
include { MERGE_ALPHA_DIVERSITY }       from    '../../modules/local/merge_alpha_diversity.nf'
include { BETA_RAREFACTION }            from    '../../modules/local/beta_rarefaction.nf'


workflow DIVERSITY_ANALYSIS {
    take:
    metadata
    qiime_sequences
    biom_without_taxonomy
    qiime_asv_table
    sample_size

    main:
    ch_versions = Channel.empty()   

    FASTTREE(qiime_sequences)
    ch_versions = ch_versions.mix(FASTTREE.out.versions)

    MINMAX( biom_without_taxonomy, sample_size )
    ch_versions = ch_versions.mix(MINMAX.out.versions)

    ALPHA_RAREFACTION(
        qiime_asv_table,
        FASTTREE.out.rooted_tree_qza,
        MINMAX.out.maxcount,
        sample_size
    )
    ch_versions = ch_versions.mix(ALPHA_RAREFACTION.out.versions)

    CORE_DIVERSITY(
        metadata,
        FASTTREE.out.rooted_tree_qza,
        qiime_asv_table,
        MINMAX.out.mincount,
        sample_size
    )
    ch_versions = ch_versions.mix(CORE_DIVERSITY.out.versions)

    MERGE_ALPHA_DIVERSITY( CORE_DIVERSITY.out.alpha, sample_size )
    ch_versions = ch_versions.mix(MERGE_ALPHA_DIVERSITY.out.versions)

    BETA_RAREFACTION(
        qiime_asv_table,
        metadata,
        MINMAX.out.mincount,
        FASTTREE.out.rooted_tree_qza,
        sample_size
    )
    ch_versions = ch_versions.mix(BETA_RAREFACTION.out.versions)

    emit:
    rooted_tree_newick              = FASTTREE.out.rooted_tree_newick.ifEmpty(null)
    rooted_tree_qza                 = FASTTREE.out.rooted_tree_qza.ifEmpty(null)
    alpha_div_metrics               = ALPHA_RAREFACTION.out.alpha_div_metrics.ifEmpty(null)
    merged_alpha_div                = MERGE_ALPHA_DIVERSITY.out.alpha_div.ifEmpty(null)
    alpha_div_shannon               = ALPHA_RAREFACTION.out.shannon_file.ifEmpty(null)
    qiime_alpha_file                = ALPHA_RAREFACTION.out.qiime_alpha_file.ifEmpty(null)
    beta_div_metrics                = CORE_DIVERSITY.out.beta_div_metrics.ifEmpty(null)
    beta_rarefaction_qiime          = BETA_RAREFACTION.out.beta_rarefaction_qiime.ifEmpty(null)
    wunifrac_matrix                 = CORE_DIVERSITY.out.wunifrac_matrix.ifEmpty(null)
    weighted_unifrac_tree           = BETA_RAREFACTION.out.weighted_unifrac_tree.ifEmpty(null)
    versions                        = ch_versions
}
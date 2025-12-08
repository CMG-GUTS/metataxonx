/*

    COMPUTATION OF BETA & ALPHA DIVERSITY

*/
include { FASTTREE } from '../../modules/local/fasttree.nf'
include { MINMAX } from '../../modules/local/minmax.nf'
include { ALPHA_RAREFACTION } from '../../modules/local/alpha_rarefaction.nf'
include { CORE_DIVERSITY } from '../../modules/local/core_diversity.nf'
include { MERGE_ALPHA_DIVERSITY } from '../../modules/local/merge_alpha_diversity.nf'
include { BETA_RAREFACTION } from '../../modules/local/beta_rarefaction.nf'


workflow DIVERSITY_ANALYSIS {
    take:
    metadata
    qiime_sequences
    biom_without_taxonomy
    qiime_asv_table

    main:
    ch_versions = Channel.empty()   

    FASTTREE( qiime_sequences )
    ch_versions = ch_versions.mix(FASTTREE.out.versions)

    MINMAX( biom_without_taxonomy )
    ch_versions = ch_versions.mix(MINMAX.out.versions)

    MINMAX.out.maxcount.view()
    MINMAX.out.mincount.view()

    ALPHA_RAREFACTION(
        qiime_asv_table,
        FASTTREE.out.rooted_tree_qza,
        MINMAX.out.maxcount
    )
    ch_versions = ch_versions.mix(ALPHA_RAREFACTION.out.versions)

    CORE_DIVERSITY(
        metadata,
        FASTTREE.out.rooted_tree_qza,
        qiime_asv_table,
        MINMAX.out.mincount
    )
    ch_versions = ch_versions.mix(CORE_DIVERSITY.out.versions)

    MERGE_ALPHA_DIVERSITY( CORE_DIVERSITY.out.alpha )
    ch_versions = ch_versions.mix(MERGE_ALPHA_DIVERSITY.out.versions)

    BETA_RAREFACTION(
        qiime_asv_table,
        metadata,
        MINMAX.out.mincount,
        FASTTREE.out.rooted_tree_qza
    )
    ch_versions = ch_versions.mix(BETA_RAREFACTION.out.versions)

    emit:
    rooted_tree_newick              = FASTTREE.out.rooted_tree_newick
    rooted_tree_qza                 = FASTTREE.out.rooted_tree_qza
    alpha_div_metrics               = ALPHA_RAREFACTION.out.alpha_div_metrics
    merged_alpha_div                = MERGE_ALPHA_DIVERSITY.out.alpha_div
    alpha_div_shannon               = ALPHA_RAREFACTION.out.shannon_file
    qiime_alpha_file                = ALPHA_RAREFACTION.out.qiime_alpha_file
    beta_div_metrics                = CORE_DIVERSITY.out.beta_div_metrics
    beta_rarefaction_qiime          = BETA_RAREFACTION.out.beta_rarefaction_qiime
    wunifrac_matrix                 = CORE_DIVERSITY.out.wunifrac_matrix
    weighted_unifrac_tree           = BETA_RAREFACTION.out.weighted_unifrac_tree
    versions                        = ch_versions
}
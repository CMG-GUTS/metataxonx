/*

    

*/
include { CREATE_ANALYSIS_MAPPING } from '../../modules/local/create_analysis_mapping.nf'
include { BIOTAVIZ } from '../../modules/local/biotaviz.nf'
include { SANKEYPLOTS } from '../../modules/local/sankeyplots.nf'
include { OMICFLOW } from '../../modules/local/omicsflow.nf'


workflow REPORT {
    take:
    biom
    metadata
    rooted_tree_newick
    weighted_unifrac_file
    shannon_file
    reads_stats_file
    dada2_errors

    main:
    ch_versions = Channel.empty() 

    if (metadata) {
        CREATE_ANALYSIS_MAPPING(
            metadata
        ).mapping.set{ metadata_ch }
    }

    BIOTAVIZ(biom)
    ch_versions = ch_versions.mix(BIOTAVIZ.out.versions)        

    SANKEYPLOTS(
        metadata_ch,
        BIOTAVIZ.out.biotaviz_relative
    )
    ch_versions = ch_versions.mix(SANKEYPLOTS.out.versions)

    OMICFLOW(
        metadata_ch,
        biom,
        rooted_tree_newick,
        weighted_unifrac_file,
        shannon_file,
        Channel.empty(),
        SANKEYPLOTS.out.sankey_image,
        dada2_errors
    )
    ch_versions = ch_versions.mix(OMICFLOW.out.versions)

    emit:
    biotaviz_clean                  = BIOTAVIZ.out.biotaviz_clean
    biotaviz_relAbun                = BIOTAVIZ.out.biotaviz_relative
    asv_table_with_taxonomy         = BIOTAVIZ.out.asv_table_with_taxonomy_test
    sankey_html                     = SANKEYPLOTS.out.sankey_html
    sankey_png                      = SANKEYPLOTS.out.sankey_image
    report                          = OMICFLOW.out.report
    versions                        = ch_versions
}
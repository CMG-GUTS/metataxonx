/*

    REPORT CREATION

*/

include { CREATE_ANALYSIS_MAPPING }     from    '../../modules/local/create_analysis_mapping.nf'
include { BIOTAVIZ }                    from    '../../modules/local/biotaviz.nf'
include { SANKEYPLOTS }                 from    '../../modules/local/sankeyplots.nf'
include { MULTIQC }                     from    '../../modules/local/multiqc.nf'
include { OMICFLOW }                    from    '../../modules/local/omicflow.nf'
include { softwareVersionsToYAML }      from    '../../subworkflows/nf-core/nf_pipeline_utils.nf'

workflow REPORT {
    take:
    biom
    metadata
    sample_size
    rooted_tree_newick
    ch_multiqc_files
    ch_versions

    main:
    omicflow_report = Channel.empty()
    sankeyplots_ch = Channel.empty()  

    BIOTAVIZ(biom)
    ch_versions = ch_versions.mix(BIOTAVIZ.out.versions)

    if (metadata) {
        CREATE_ANALYSIS_MAPPING(
            metadata
        ).mapping.set{ metadata_ch }

        sankeyplots_ch = SANKEYPLOTS(
            metadata,
            BIOTAVIZ.out.biotaviz_relative
        )
        ch_versions = ch_versions.mix(SANKEYPLOTS.out.versions)

        OMICFLOW(
            metadata_ch,
            biom,
            rooted_tree_newick,
            sample_size
        ).report.set{ omicflow_report }
        ch_versions = ch_versions.mix(OMICFLOW.out.versions)
    }

    // Create software versions yaml file
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'metataxonx_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC(
        ch_multiqc_files.collect(),
        params.multiqc_config,
        [], [], [], []
    )

    emit:
    biotaviz_clean                  = BIOTAVIZ.out.biotaviz_clean
    biotaviz_relAbun                = BIOTAVIZ.out.biotaviz_relative
    asv_table_with_taxonomy         = BIOTAVIZ.out.asv_table_with_taxonomy_test
    sankey_html                     = sankeyplots_ch.sankey_html.ifEmpty(null)
    sankey_png                      = sankeyplots_ch.sankey_image.ifEmpty(null)
    technical_report                = MULTIQC.out.report
    analysis_report                 = omicflow_report
    versions                        = ch_versions
}
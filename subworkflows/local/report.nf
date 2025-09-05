/*

    REPORT CREATION

*/

include { CREATE_ANALYSIS_MAPPING } from '../../modules/local/create_analysis_mapping.nf'
include { BIOTAVIZ } from '../../modules/local/biotaviz.nf'
include { SANKEYPLOTS } from '../../modules/local/sankeyplots.nf'
include { MULTIQC } from '../../modules/local/multiqc.nf'
include { OMICFLOW } from '../../modules/local/omicflow.nf'

import org.yaml.snakeyaml.Yaml

workflow REPORT {
    take:
    biom
    metadata
    rooted_tree_newick
    ch_multiqc_files
    ch_versions

    main:

    BIOTAVIZ(biom)
    ch_versions = ch_versions.mix(BIOTAVIZ.out.versions)        

    if (metadata) {
        CREATE_ANALYSIS_MAPPING(
            metadata
        ).mapping.set{ metadata_ch }

        sankeplots_ch = SANKEYPLOTS(
            metadata_ch,
            BIOTAVIZ.out.biotaviz_relative
        )
        ch_versions = ch_versions.mix(SANKEYPLOTS.out.versions)

        OMICFLOW(
            metadata_ch,
            biom,
            rooted_tree_newick
        ).report.set{ omicflow_report }
        ch_versions = ch_versions.mix(OMICFLOW.out.versions)

    } else {
        omicflow_report = Channel.empty()
        sankeplots_ch = Channel.empty()
    }

    // Combine all versions
    ch_versions_parsed = ch_versions.collect().map { fileList ->
        def all_versions = []
        fileList.each { filePath ->
            File f = filePath.toFile()
            if (f.exists()) {
                def versions = new Yaml().load(f.text)
                all_versions += versions
            }
        }
        return all_versions.unique()
    }

    MULTIQC(
        ch_multiqc_files,
        params.multiqc_config,
        ch_versions_parsed,
        [], [], [], []
    )

    emit:
    biotaviz_clean                  = BIOTAVIZ.out.biotaviz_clean
    biotaviz_relAbun                = BIOTAVIZ.out.biotaviz_relative
    asv_table_with_taxonomy         = BIOTAVIZ.out.asv_table_with_taxonomy_test
    sankey_html                     = sankeplots_ch.sankey_html.ifEmpty(null)
    sankey_png                      = sankeplots_ch.sankey_image.ifEmpty(null)
    technical_report                = MULTIQC.out.report
    analysis_report                 = omicflow_report
    versions                        = ch_versions
}
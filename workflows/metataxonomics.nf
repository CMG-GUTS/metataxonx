/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { save_output } from        '../lib/utils.groovy'

include { CHECK_INPUT } from        '../subworkflows/local/check_input.nf'
include { PREPROCESSING } from      '../subworkflows/local/preprocessing.nf'
include { DENOISE } from            '../subworkflows/local/denoise.nf'
include { TAXONOMY } from           '../subworkflows/local/taxonomy.nf'
include { DIVERSITY_ANALYSIS } from '../subworkflows/local/diversity_analysis.nf'
include { REPORT } from             '../subworkflows/local/report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METATAXONOMICS {

    CHECK_INPUT ()

    if (params.input || params.reads) {

        PREPROCESSING(
            CHECK_INPUT.out.meta
        )
        
        if (params.save_trim_reads & !params.bypass_trim) save_output(PREPROCESSING.out.trimmed, "preprocessing/trimmed")
        // Save output that cannot be bypassed
        if (params.save_merged_reads) save_output(PREPROCESSING.out.merged, "preprocessing/merged")
        if (params.save_sampled_reads) save_output(PREPROCESSING.out.subsampled_reads, "preprocessing/subsampled")

        DENOISE( PREPROCESSING.out.subsampled_reads )

        if (params.save_denoise_stats) {
            save_output(DENOISE.out.denoise_stats, "denoise/stats")
            save_output(DENOISE.out.dada2_errors, "denoise/errors")
        }
        if (params.save_denoise_sequences) save_output(DENOISE.out.sequences, "denoise/sequences")
        if (params.save_denoise_table) save_output(DENOISE.out.counts_table, "denoise/counts_without_taxonomy")

        TAXONOMY(
            DENOISE.out.counts_table,
            DENOISE.out.sequences,
            CHECK_INPUT.out.classifier
        )

        if (params.save_biom_files) {
            save_output(TAXONOMY.out.biom_without_taxonomy, "taxonomy/biom")
            save_output(TAXONOMY.out.biom_with_taxonomy, "taxonomy/biom")
        }

        if (!params.bypass_post_analysis) {
            DIVERSITY_ANALYSIS(
                DENOISE.out.dada2_mapping_file,
                TAXONOMY.out.qiime_sequences_qza,
                TAXONOMY.out.biom_without_taxonomy,
                TAXONOMY.out.qiime_asv_table
            )

            if (params.save_alpha_div_files) {
                save_output(DIVERSITY_ANALYSIS.out.alpha_div_metrics, "alpha_diversity/metrics")
                save_output(DIVERSITY_ANALYSIS.out.merged_alpha_div, "alpha_diversity/metrics")
            }
            if (params.save_beta_div_files) {
                save_output(DIVERSITY_ANALYSIS.out.beta_div_metrics, "beta_diversity/metrics")
            }
            if (params.save_tree_files) {
                save_output(DIVERSITY_ANALYSIS.out.rooted_tree_newick, "diversity_analysis/tree")
                save_output(DIVERSITY_ANALYSIS.out.weighted_unifrac_tree, "diversity_analysis/tree")
            }

            // Save qiime artifacts
            if (params.save_qiime_artifacts) {
                save_output(TAXONOMY.out.taxonomy_file_tsv, "qiime_artifact/taxonomy")
                save_output(TAXONOMY.out.taxonomy_file_qza, "qiime_artifact/taxonomy")
                save_output(TAXONOMY.out.qiime_asv_table, "qiime_artifact/taxonomy")
                save_output(TAXONOMY.out.qiime_sequences_qza, "qiime_artifact/taxonomy")
                save_output(DIVERSITY_ANALYSIS.out.rooted_tree_qza, "qiime_artifact/tree")
                save_output(DIVERSITY_ANALYSIS.out.qiime_alpha_file, "qiime_artifact/alpha_div")
                save_output(DIVERSITY_ANALYSIS.out.beta_rarefaction_qiime, "qiime_artifact/beta_div")
            }
        }

        ch_versions = PREPROCESSING.out.versions
            .mix(DENOISE.out.versions)
            .mix(TAXONOMY.out.versions)
            .mix(DIVERSITY_ANALYSIS.out.versions)
            .collect()

        ch_multiqc_files = PREPROCESSING.out.multiqc_files
            .mix(DENOISE.out.dada2_errors)
            .mix(DENOISE.out.denoise_stats)
            .collect()

        if (!params.bypass_report) {
            REPORT(
                TAXONOMY.out.biom_with_taxonomy,
                CHECK_INPUT.out.metadata,
                DIVERSITY_ANALYSIS.out.rooted_tree_newick,
                ch_multiqc_files,
                ch_versions
            )

            if (params.save_biotaviz_files) {
                save_output(REPORT.out.biotaviz_clean, "biotiaviz")
                save_output(REPORT.out.biotaviz_relAbun, "biotiaviz")
                save_output(REPORT.out.asv_table_with_taxonomy, "biotiaviz")
            }

            if (params.save_sankey_plot) {
                save_output(REPORT.out.sankey_html, "sankey")
                save_output(REPORT.out.sankey_png, "sankey")
            }

            if (params.save_final_reports) {
                save_output(REPORT.out.analysis_report, "report")
                save_output(REPORT.out.technical_report, "report")
            }
        }        
    }
    
}
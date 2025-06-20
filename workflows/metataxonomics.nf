/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_INPUT } from '../subworkflows/local/check_input.nf'
include { PREPROCESSING } from '../subworkflows/local/preprocessing.nf'
include { DENOISE } from '../subworkflows/local/denoise.nf'
include { TAXONOMY } from '../subworkflows/local/taxonomy.nf'
include { DIVERSITY_ANALYSIS } from '../subworkflows/local/diversity_analysis.nf'
include { REPORT } from '../subworkflows/local/report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAPIPE {
    
    CHECK_INPUT ()

    PREPROCESSING(
        CHECK_INPUT.out.meta
    )
    
    if (params.save_trim_reads & !params.bypass_trim) save_output(PREPROCESSING.out.trimmed, "preprocessing/trimmed")
    // Save output that cannot be bypassed
    if (params.save_merged_reads) save_output(PREPROCESSING.out.merged, "preprocessing/merged")
    if (params.save_sampled_reads) save_output(PREPROCESSING.out.subsampled_reads, "preprocessing/subsampled")
    if (params.save_multiqc_reports) {
        save_output(PREPROCESSING.out.multiqc_report, "preprocessing/multiqc")
        save_output(PREPROCESSING.out.read_stats, "preprocessing/multiqc")
    }

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
        save_output(DIVERSITY_ANALYSIS.out.beta_rarefaction_files, "beta_diversity/metrics")
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

    REPORT(

    )
}

def save_output(input_ch, sub_dir_name) {
    input_ch.map { item ->
        def (meta, files) = (item instanceof List && item.size() == 2) ? [item[0], item[1]] : [null, item]
        def outDir = file("${params.outdir}/${sub_dir_name}")
        outDir.mkdir()
        if (files.size() == 2) {
            files.each { inputFile ->
               file(inputFile).copyTo(file("${outDir}/${file(inputFile).getName()}"))
            }
        } else {
            file(files).copyTo(file("${outDir}/${file(files).getName()}"))
        }
    }
}
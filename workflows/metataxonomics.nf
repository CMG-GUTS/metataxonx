/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_INPUT } from '../subworkflows/local/check_input.nf'
include { PREPROCESSING } from '../subworkflows/local/preprocessing.nf'
include { DENOISE } from '../subworkflows/local/denoise.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAPIPE {
    
    CHECK_INPUT ()

    PREPROCESSING(
        CHECK_INPUT.out.meta,
        CHECK_INPUT.out.classifier
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
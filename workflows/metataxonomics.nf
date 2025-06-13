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

    DENOISE(

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
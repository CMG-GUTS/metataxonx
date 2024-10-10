process multiqc_untrimmed {
// run multiqc
    container "$projectDir/containers/singularity/pyrrr.sif"
    publishDir "${params.outdir}/untrimmed", mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc --force --interactive .
    """
}
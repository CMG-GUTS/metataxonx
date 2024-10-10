process multiqc_pear {
// run multiqc
    container "$projectDir/containers/singularity/pyrrr.sif"
    publishDir "${params.outdir}/multiqc_pear", mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc ./
    """
}
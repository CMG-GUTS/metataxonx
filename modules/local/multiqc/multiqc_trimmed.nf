process multiqc_trimmed {
// run multiqc
    container "$projectDir/containers/singularity/pyrrr.sif"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc --force --interactive .
    """
}
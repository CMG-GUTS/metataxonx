process multiqc_trimmed {
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
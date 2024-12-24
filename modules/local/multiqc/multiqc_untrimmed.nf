process multiqc_untrimmed {
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
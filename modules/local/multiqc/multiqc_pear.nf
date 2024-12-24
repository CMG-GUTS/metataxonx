process multiqc_pear {
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
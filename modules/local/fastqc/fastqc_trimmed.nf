process fastqc_trimmed {
    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    ls * | parallel -j ${params.cpus} --verbose "fastqc {}"
    """
}
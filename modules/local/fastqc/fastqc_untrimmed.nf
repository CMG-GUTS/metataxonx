process fastqc_untrimmed {
// run fastqc
    container "$projectDir/containers/singularity/fastqc.sif"

    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    ls * | parallel -j ${params.cpus} --verbose "fastqc {}"
    """
}
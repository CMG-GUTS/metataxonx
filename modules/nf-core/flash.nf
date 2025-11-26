process FLASH {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.extendedFrags.fastq.gz"), emit: merged
    tuple val(meta), path("*.notCombined_*.fastq.gz"), emit: notcombined
    tuple val(meta), path("*.hist")                  , emit: histogram
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$reads" == "${prefix}.extendedFrags.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$reads" == "${prefix}.notCombined_1.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$reads" == "${prefix}.notCombined_2.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    flash \\
        $args \\
        -o ${prefix} \\
        -z \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flash: \$(echo \$(flash --version 2>&1) | sed 's/^.*FLASH v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.extendedFrags.fastq.gz
    touch ${prefix}.hist
    touch ${prefix}.notCombined_1.fastq.gz
    touch ${prefix}.notCombined_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flash: stub-version
    END_VERSIONS
    """
}
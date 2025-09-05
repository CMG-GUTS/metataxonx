process MERGE_ALPHA_DIVERSITY {
    label 'process_single'

    input:
    path(alpha)

    output:
    path "alpha_diversity.txt"          , emit: alpha_div
    path "versions.yml"                 , emit: versions

    script:
    """
    qiime metadata tabulate \\
        --m-input-file *vector.qza \\
        --o-visualization combined-alpha-metadata.qzv

    qiime tools export \\
        --input-path combined-alpha-metadata.qzv \\
        --output-path data

    cp data/metadata.tsv alpha_diversity.txt

    qiime_version=\$(qiime --version | head -1 | sed -E 's/.*version ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$qiime_version
    END_VERSIONS
    """

    stub:
    """
    touch alpha_diversity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version)
    END_VERSIONS
    """
}
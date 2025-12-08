process COMBINE_TAXONOMY {
    label 'process_single'

    input:
    path(biom)
    path(taxonomy)

    output:
    path "biom_with_taxonomy.biom"          , emit: biom_taxonomy
    path "versions.yml"                     , emit: versions

    script:
    """
    biom add-metadata \\
        -i $biom \\
        -o biom_with_taxonomy.biom \\
        --observation-metadata-fp $taxonomy \\
        --sc-separated taxonomy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biom: \$(biom --version | sed -e "s/biom, version //g")
    END_VERSIONS
    """

    stub:
    """
    touch biom_with_taxonomy.biom

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biom: stub-version
    END_VERSIONS
    """
}
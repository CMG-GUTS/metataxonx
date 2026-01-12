process SANKEYPLOTS {
    label 'process_single'

    input:
    path(mapping)
    path(biotaviz)

    output:
    path "*.html"           , emit: sankey_html
    path "*.png"            , emit: sankey_image
    path "versions.yml"     , emit: versions

    script:
    """
    sankey-file-prep \\
        --taxa-filter 0.01 \\
        --sample-repeat false \\
        --combine-rankstat false \\
        -m ${mapping} \\
        -i ${biotaviz}

    sankey-diagram-html-generator biotaviz_sankey_prepfile-AverageAllSamples.csv
    sankey-diagram-png-generator biotaviz_sankey_prepfile-AverageAllSamples.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | cut -d ' ' -f3)
        python: \$(python3.11 --version 2>&1 | sed -e "s/Python //g")
    END_VERSIONS
    """

    stub:
    """
    touch image.png
    touch image.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub-version
        python: stub-version
    END_VERSIONS
    """
}
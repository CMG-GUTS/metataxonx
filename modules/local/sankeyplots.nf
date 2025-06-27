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
    python3.11 $projectDir/bin/python/sankey-file-prep.py \\
        --taxa-filter 0.01 \\
        --sample-repeat false \\
        --combine-rankstat false \\
        -m ${mapping} \\
        -i ${biotaviz}

    Rscript $projectDir/bin/R/sankey-diagram-html-generator.R biotaviz_sankey_prepfile-AverageAllSamples.csv
    Rscript $projectDir/bin/R/sankey-diagram-png-generator.R biotaviz_sankey_prepfile-AverageAllSamples.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """

    stub:
    """
    touch image.png
    touch image.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """
}
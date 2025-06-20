process OMICFLOW {
    label 'process_low'

    input:
    path(metadata_clean)
    path(biom_taxonomy)
    path(rooted_tree_newick)
    path(wunifrac)
    path(shannon)
    path(read_stats)
    path(sankey_image)
    path(errProfile)

    output:
    path "report.html"              , emit: report
    path "versions.yml"             , emit: versions

    script:
    """
    Rscript $projectDir/bin/R/OmicFlow/autoFlow.R \\
        --metadata ${metadata_clean} \\
        --biom ${biom_taxonomy} \\
        --tree ${rooted_tree_newick} \\
        --cpus ${task.cpus} \\
        --i-beta-div ${wunifrac} \\
        --i-alpha-div ${shannon}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """

    stub:
    """
    touch report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """
}
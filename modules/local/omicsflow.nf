process omics_analysis {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(metadata_clean)
    file(biom_taxonomy)
    file(rooted_tree_newick)
    file(wunifrac)
    file(shannon)
    file(read_stats)
    file(sankey_image)
    file(errProfile)

    output:
    file("report.html")

    script:
    """
    Rscript $projectDir/bin/R/OmicFlow/autoFlow.R \
        --metadata ${metadata_clean} \
        --biom ${biom_taxonomy} \
        --tree ${rooted_tree_newick} \
        --cpus ${params.cpus} \
        --i-beta-div ${wunifrac} \
        --i-alpha-div ${shannon}
    """
}
process omics_analysis {
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(metadata_clean)
    file(biom_taxonomy)
    file(rooted_tree_newick)
    file(sequences_fasta)
    file(read_stats)
    file(sankey_image)
    file(shannon_file)

    output:
    file("report.html")

    script:
    """
    Rscript $projectDir/bin/R/OmicFlow/00_main.R \
        --metadata ${metadata_clean} \
        --biom ${biom_taxonomy} \
        --tree ${rooted_tree_newick}
    """
}
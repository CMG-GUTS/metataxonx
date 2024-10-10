process combine_taxonomy_biom {
// add taxonomy info to biom file
    container "$projectDir/containers/singularity/qiime2.sif"

    publishDir params.outdir, mode: 'copy'

    input:
    path(biomfile)
    path(taxonomyfile)

    output:
    path("biom_with_taxonomy.biom"), emit: biom_taxonomy

    script:
    """
    biom add-metadata \
        -i ${biomfile} \
        -o biom_with_taxonomy.biom \
        --observation-metadata-fp ${taxonomyfile} \
        --sc-separated taxonomy
    """
}
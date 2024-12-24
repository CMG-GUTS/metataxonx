process combine_taxonomy_biom {
    input:
    path(biomfile)
    path(taxonomyfile)

    output:
    path("biom_with_taxonomy.biom"), emit: biom_taxonomy

    publishDir params.outdir, mode: 'copy'

    script:
    """
    biom add-metadata \
        -i ${biomfile} \
        -o biom_with_taxonomy.biom \
        --observation-metadata-fp ${taxonomyfile} \
        --sc-separated taxonomy
    """
}
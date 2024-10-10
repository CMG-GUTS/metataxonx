process beta_rarefaction {
    // beta rarefaction & tree building
    container "$projectDir/containers/singularity/qiime2.sif"
    cpus = params.cpus

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.newick", mode: 'copy'
    publishDir params.outdir, pattern: "beta_rarefaction/*", mode: 'copy'

    input:
    file(table)
    file(metadata)
    file(mincount)
    file(tree)

    output:
    file("beta_rarefaction.qzv")
    file("weighted_unifrac_tree.newick")
    path("beta_rarefaction/*")

    script:
    """

    nice -${params.niceness} qiime diversity beta-rarefaction \
        --i-table ${table} \
        --i-phylogeny ${tree} \
        --p-metric weighted_unifrac \
        --p-clustering-method upgma \
        --m-metadata-file ${metadata} \
        --p-sampling-depth "\$(<${mincount})"  \
        --o-visualization beta_rarefaction.qzv

    qiime tools export \
        --input-path beta_rarefaction.qzv \
		--output-path beta_rarefaction

    cp beta_rarefaction/sample-clustering-upgma.tre ./weighted_unifrac_tree.newick
    """
}
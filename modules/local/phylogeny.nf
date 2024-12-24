process phylogeny {
    input:
    file(sequences)

    output:
    path("rooted_tree.newick"), emit: rooted_tree_newick
    path("rooted-tree.qza"), emit: rooted_tree_qza

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.newick", mode: 'copy'

    script:
    """
    nice -${params.niceness} qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${sequences} \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --p-n-threads ${params.cpus} \
        --o-rooted-tree rooted-tree.qza

    qiime tools export \
        --input-path rooted-tree.qza  \
		--output-path ./

	mv tree.nwk rooted_tree.newick
    """
}
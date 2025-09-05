process FASTTREE {
    label 'process_high'

    input:
    path(sequences)

    output:
    path "rooted_tree.newick"       , emit: rooted_tree_newick
    path "rooted-tree.qza"          , emit: rooted_tree_qza
    path "versions.yml"             , emit: versions

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \\
        --i-sequences ${sequences} \\
        --o-alignment aligned-rep-seqs.qza \\
        --o-masked-alignment masked-aligned-rep-seqs.qza \\
        --o-tree unrooted-tree.qza \\
        --p-n-threads ${task.cpus} \\
        --o-rooted-tree rooted-tree.qza

    qiime tools export \\
        --input-path rooted-tree.qza  \\
        --output-path ./

    mv tree.nwk rooted_tree.newick

    qiime_version=\$(qiime --version | head -1 | sed -E 's/.*version ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$qiime_version
    END_VERSIONS
    """

    stub:
    """
    touch rooted_tree.newick
    touch rooted-tree.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version)
    END_VERSIONS
    """
}
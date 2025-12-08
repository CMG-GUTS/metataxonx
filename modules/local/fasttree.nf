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
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime phylogeny align-to-tree-mafft-fasttree \\
        --i-sequences ${sequences} \\
        --o-alignment aligned-rep-seqs.qza \\
        --o-masked-alignment masked-aligned-rep-seqs.qza \\
        --o-tree unrooted-tree.qza \\
        --p-n-threads ${task.cpus} \\
        --o-rooted-tree rooted-tree.qza

    qiime tools export \\
        --input-path rooted-tree.qza  \\
        --output-path phylogenetic_tree

    mv phylogenetic_tree/tree.nwk rooted_tree.newick

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1 | sed -e "s/q2cli version //g")
    END_VERSIONS
    """

    stub:
    """
    touch rooted_tree.newick
    touch rooted-tree.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: stub-version
    END_VERSIONS
    """
}
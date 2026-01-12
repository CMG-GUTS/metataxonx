process BETA_RAREFACTION {
    label 'process_medium'

    input:
    path(table)
    path(metadata)
    path(mincount)
    path(tree)

    output:
    path "beta_rarefaction.qzv"                 , emit: beta_rarefaction_qiime
    path "weighted_unifrac_tree.newick"         , emit: weighted_unifrac_tree
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime diversity beta-rarefaction \\
        --i-table ${table} \\
        --i-phylogeny ${tree} \\
        --p-metric weighted_unifrac \\
        --p-clustering-method upgma \\
        --m-metadata-file ${metadata} \\
        --p-sampling-depth "\$(<${mincount})"  \\
        --o-visualization beta_rarefaction.qzv

    qiime tools export \\
        --input-path beta_rarefaction.qzv \\
        --output-path beta_rarefaction

    cp beta_rarefaction/sample-clustering-upgma.tre ./weighted_unifrac_tree.newick

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1 | sed -e "s/q2cli version //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch beta_rarefaction.qzv
    touch weighted_unifrac_tree.newick
    mkdir beta_rarefaction
    touch beta_rarefaction/example.txt   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: stub-version
    END_VERSIONS
    """
}
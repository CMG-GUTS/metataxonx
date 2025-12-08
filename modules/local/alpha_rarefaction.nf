process ALPHA_RAREFACTION {
    label 'process_medium'

    input:
    path(table)
    path(phylogeny)
    path(maxcount)

    output:
    path "alpha_rarefaction.qzv"            , emit: qiime_alpha_file
    path "alpha_rarefaction/*.csv"          , emit: alpha_div_metrics
    path "alpha_rarefaction/shannon.csv"    , emit: shannon_file
    path "versions.yml"                     , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime diversity alpha-rarefaction \\
        --i-table ${table} \\
        --i-phylogeny ${phylogeny} \\
        --p-max-depth "\$(<${maxcount})" \\
        --o-visualization alpha_rarefaction.qzv

    qiime tools export \\
        --input-path alpha_rarefaction.qzv \\
        --output-path alpha_rarefaction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1 | sed -e "s/q2cli version //g")
    END_VERSIONS
    """

    stub:
    """
    touch alpha_rarefaction.qzv
    mkdir alpha_rarefaction
    touch alpha_rarefaction/shannon.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: stub-version
    END_VERSIONS
    """
}
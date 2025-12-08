process MERGE_ALPHA_DIVERSITY {
    label 'process_single'

    input:
    path(alpha)

    output:
    path "alpha_diversity.txt"          , emit: alpha_div
    path "versions.yml"                 , emit: versions

    script:
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime metadata tabulate \\
        --m-input-file *vector.qza \\
        --o-visualization combined-alpha-metadata.qzv

    qiime tools export \\
        --input-path combined-alpha-metadata.qzv \\
        --output-path data

    cp data/metadata.tsv alpha_diversity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1 | sed -e "s/q2cli version //g")
    END_VERSIONS
    """

    stub:
    """
    touch alpha_diversity.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: stub-version
    END_VERSIONS
    """
}
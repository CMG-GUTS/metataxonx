process CORE_DIVERSITY {   
    label 'process_medium'

    input:
    path(metadata)
    path(tree)
    path(table)
    path(mincount)

    output:
    path "diversity_core/*vector.qza"               , emit: alpha
    path "*distance*.txt"                           , emit: beta_div_metrics
    path "weighted_unifrac_distance_matrix.qza.txt" , emit: wunifrac_matrix
    path "versions.yml"                             , emit: versions
    
    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"
    
    qiime diversity core-metrics-phylogenetic \\
        --m-metadata-file ${metadata} \\
        --i-phylogeny ${tree} \\
        --i-table ${table} \\
        --p-sampling-depth "\$(<${mincount})" \\
        --output-dir diversity_core \\
        --p-n-jobs-or-threads ${task.cpus} \\
        --quiet

    for i in diversity_core/*distance*; do \\
            qiime tools export \\
                --input-path \$i \\
                --output-path \$i.data; \\
            mv \$i.data/distance-matrix.tsv ./\$i.txt; done

    mv diversity_core/*.txt .

    qiime_version=\$(qiime --version | head -1 | sed -E 's/.*version ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$qiime_version
    END_VERSIONS
    """

    stub:
    """
    mkdir diversity_core
    touch diversity_core/alpha_vector.qza
    touch weighted_unifrac_distance_matrix.qza.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: stub-version
    END_VERSIONS
    """
}
process core_diversity {    
    input:
    file(metadata)
    file(tree)
    file(table)
    file(mincount)

    output:
    path("diversity_core/*vector.qza"), emit: alpha
    path("*.txt")
    path("weighted_unifrac_distance_matrix.qza.txt"), emit: wunifrac_matrix

    publishDir "${params.outdir}/distance_matrices/", pattern: "*.txt", mode: "copy"
    
    script:
    """
    nice -${params.niceness} qiime diversity core-metrics-phylogenetic \
		--m-metadata-file ${metadata} \
		--i-phylogeny ${tree} \
		--i-table ${table} \
        --p-sampling-depth "\$(<${mincount})" \
		--output-dir diversity_core \
		--p-n-jobs-or-threads ${params.cpus} \
		--quiet

	for i in diversity_core/*distance*; do \
            qiime tools export \
                --input-path \$i \
		        --output-path \$i.data; \
		    mv \$i.data/distance-matrix.tsv ./\$i.txt; done

    mv diversity_core/*.txt .
    """
}
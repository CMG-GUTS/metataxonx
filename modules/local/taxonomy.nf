process assign_taxonomy {
    // assign taxonomy using sklearn qiime2
    container "$projectDir/containers/singularity/qiime2.sif"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.tsv", mode: 'copy'

    input:
    file(classifier)
    file(sequences)

    output:
    path("taxonomy.tsv"), emit: taxonomy_tsv
    path("taxonomy_sklearn.qza"), emit: taxonomy_qza

    script:
    """
    nice -${params.niceness} qiime feature-classifier classify-sklearn \
         --i-classifier ${classifier} \
         --i-reads ${sequences} \
         --o-classification taxonomy_sklearn.qza \
         --p-n-jobs ${params.cpus}

    qiime metadata tabulate \
         --m-input-file taxonomy_sklearn.qza \
         --o-visualization taxonomy_sklearn.qzv

    qiime tools export --input-path taxonomy_sklearn.qza  \
         --output-path taxonomy

    qiime tools export --input-path taxonomy_sklearn.qzv  \
         --output-path taxonomy

    cat taxonomy/taxonomy.tsv | \
    sed 's/Feature ID/#OTUID/' | \
    sed 's/Taxon/taxonomy/' | \
    sed 's/Consensus/consensus/' > taxonomy/taxonomy_relabeled.tsv
    mv taxonomy/taxonomy_relabeled.tsv ./taxonomy.tsv
    """
}
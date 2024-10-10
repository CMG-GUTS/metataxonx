process merge_alpha_diversity {
    // merge alpha diversity metrics
    container "$projectDir/containers/singularity/qiime2.sif"
    cpus = params.cpus

    publishDir params.outdir, mode: 'copy'

    input:
    file(alpha)

    output:
    path("alpha_diversity.txt")

    script:
    """
    qiime metadata tabulate \
        --m-input-file *vector.qza \
        --o-visualization combined-alpha-metadata.qzv

    qiime tools export \
        --input-path combined-alpha-metadata.qzv \
		--output-path data

	cp data/metadata.tsv alpha_diversity.txt
    """
}
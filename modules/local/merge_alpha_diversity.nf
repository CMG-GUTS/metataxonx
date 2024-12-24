process merge_alpha_diversity {
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
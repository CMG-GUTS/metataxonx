process alpha_rarefaction {
    output:
    path("alpha_rarefaction.qzv")
    path("alpha_rarefaction/*")
    path("alpha_rarefaction/shannon.csv"), emit: shannon_file

    input:
    file(table)
    file(phylogeny)
    file(maxcount)

    output:
    file("alpha_rarefaction.qzv")
    file("alpha_rarefaction/*.csv")

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.csv", mode: 'copy'
    publishDir params.outdir, pattern: "alpha_rarefaction/*", mode: 'copy'

    script:
    """
    nice -${params.niceness} qiime diversity alpha-rarefaction \
        --i-table ${table} \
        --i-phylogeny ${phylogeny} \
        --p-max-depth "\$(<${maxcount})" \
        --o-visualization alpha_rarefaction.qzv

    qiime tools export \
        --input-path alpha_rarefaction.qzv \
	--output-path alpha_rarefaction
    """
}
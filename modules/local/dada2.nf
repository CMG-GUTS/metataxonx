process dada2_denoise {
    input:
    file(mapping)

    output:
    path("seq-tab.tsv"), emit: seq_table
    path("denoising-stats.tsv"), emit: denoising_stats
    path("rep-seqs.fna"), emit: rep_seqs_fasta
    path("dada_report.txt")
    path("*.png")
    path("errProfile_1.png"), emit: dada2_errors

    publishDir "${params.outdir}/qiime_artefacts/", mode: "copy"

    script:
    def novaseq = params.novaseq == true ? "--novaseq" : ""
    """
    Rscript $projectDir/bin/R/run-dada2-batch/parallel_dada2.R \
        --metadata ${mapping} \
        --batch_n ${params.batch_size} \
        --cpus ${params.cpus} \
        --p-trunc-q 2 \
        --p-max-ee 6 \
        --p-min-fold-parent-over-abundance 2 \
        --p-chimera-method consensus \
        ${novaseq} \
        > dada_report.txt
    """
}
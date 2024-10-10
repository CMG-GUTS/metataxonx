process dada2_denoise {
    // denoise with dada2
    container "$projectDir/containers/singularity/parallel_dada2.sif"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{rds}", mode: "copy"

    input:
    file(mapping)

    output:
    path("seq-tab.tsv"), emit: seq_table
    path("denoising-stats.tsv"), emit: denoising_stats
    path("rep-seqs.fna"), emit: rep_seqs_fasta
    path("dada_report.txt")

    script:
    """
    Rscript $projectDir/bin/R/run-dada2-batch/parallel_dada2.R \
        --metadata ${mapping} \
        --batch_n ${params.batch_size} \
        --cpus ${params.cpus} \
        --p-trunc-q 2 \
        --p-max-ee 6 \
        --p-min-fold-parent-over-abundance 2 \
        --p-chimera-method consensus \
        > dada_report.txt
    """
}
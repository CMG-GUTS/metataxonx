process qiime_import_export {
    // import sequences from fastq and metadata file
    container "$projectDir/containers/singularity/qiime2.sif"

    publishDir "${params.outdir}/qiime_artefacts/", mode: "copy"

    input:
    file(seq_table)
    file(rep_seqs_fasta)

    output:
    // path('representative_sequences.fasta'), emit: sequences_fasta
    path("table.qza"), emit: asv_qza
    path('asv_table_no_taxonomy.biom'), emit: biom_no_taxonomy
    path('rep-seqs.qza'), emit: repseqs_qza

    script:
    """
    # Converts phyloseq biom to qiime compatible biom file
    biom convert -i ${seq_table} -o table.biom --to-hdf5

    # Convert RDS files from parallel_dada2 into qiime2 formats
    qiime tools import \
    --input-path table.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path table.qza    

    # convert absolute (raw) count table to biom format
    qiime tools export \
        --input-path table.qza \
        --output-path ./
    mv feature-table.biom asv_table_no_taxonomy.biom

    # representative sequences
    qiime tools import \
    --input-path ${rep_seqs_fasta} \
    --type 'FeatureData[Sequence]' \
    --output-path rep-seqs.qza
    """
}
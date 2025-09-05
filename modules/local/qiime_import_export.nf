process QIIME_IMPORT_EXPORT {
    tag "BATCH"
    label 'process_low' 

    input:
    path(seq_table)
    path(rep_seqs_fasta)

    output:
    // path('representative_sequences.fasta'), emit: sequences_fasta
    path "table.qza"                        , emit: asv_qza
    path "asv_table_no_taxonomy.biom"       , emit: biom_no_taxonomy
    path "rep-seqs.qza"                     , emit: repseqs_qza
    path "versions.yml"                     , emit: versions

    script:
    """
    # Converts phyloseq biom to qiime compatible biom file
    biom convert -i ${seq_table} -o table.biom --to-hdf5

    # Convert RDS files from parallel_dada2 into qiime2 formats
    qiime tools import \\
    --input-path table.biom \\
    --type 'FeatureTable[Frequency]' \\
    --input-format BIOMV210Format \\
    --output-path table.qza    

    # convert absolute (raw) count table to biom format
    qiime tools export \\
        --input-path table.qza \\
        --output-path ./
    mv feature-table.biom asv_table_no_taxonomy.biom

    # representative sequences
    qiime tools import \\
    --input-path ${rep_seqs_fasta} \\
    --type 'FeatureData[Sequence]' \\
    --output-path rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1)
        biom: \$(biom --version | head -1)
    END_VERSIONS

    sed -i.bak -E '
    /^ *biom:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    /^ *qiime:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    ' versions.yml

    """

    stub:
    """
    touch table.qza
    touch asv_table_no_taxonomy.biom
    touch rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version)
        biom: \$(biom --version)
    END_VERSIONS
    """
}
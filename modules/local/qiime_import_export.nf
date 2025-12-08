process QIIME_IMPORT_EXPORT {
    tag "BATCH"
    label 'process_low' 

    input:
    path(seq_table)
    path(rep_seqs_fasta)

    output:
    path "asvs_without_taxonomy.qza"        , emit: asv_qza
    path "asvs_without_taxonomy.biom"       , emit: biom_no_taxonomy
    path "rep-seqs.qza"                     , emit: repseqs_qza
    path "versions.yml"                     , emit: versions

    script:
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # Converts phyloseq biom to qiime compatible biom file
    biom convert -i ${seq_table} -o asvs_without_taxonomy.biom --to-hdf5

    # Convert RDS files from parallel_dada2 into qiime2 formats
    qiime tools import \\
    --input-path asvs_without_taxonomy.biom \\
    --type 'FeatureTable[Frequency]' \\
    --input-format BIOMV210Format \\
    --output-path asvs_without_taxonomy.qza    

    # representative sequences
    qiime tools import \\
    --input-path ${rep_seqs_fasta} \\
    --type 'FeatureData[Sequence]' \\
    --output-path rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1 | sed -e "s/q2cli version //g")
        biom: \$(biom --version | sed -e "s/biom, version //g")
    END_VERSIONS

    """

    stub:
    """
    touch asvs_without_taxonomy.qza
    touch asvs_without_taxonomy.biom
    touch rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: stub-version
        biom: stub-version
    END_VERSIONS
    """
}
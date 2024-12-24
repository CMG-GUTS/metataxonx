process pear {
    input:
    file(fastqfile)
    file(metadata_clean)

    output:
    path("assembled/*.fastq.gz"), emit: assembled
    path("pear_counts.txt"), emit: counts
    path("mapping.txt"), emit: qiime_metadata

    script:
    """
    nice -${params.niceness} python3.11 $projectDir/bin/python/run_pear.py \
        -m ${metadata_clean} \
        --cpu ${params.cpus} \
        --sub-sampling ${params.read_sample_size}
    """
}
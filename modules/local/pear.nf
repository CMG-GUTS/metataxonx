process pear {
    // run Paired-End reAd mergeR
    container "$projectDir/containers/singularity/pyrrr.sif"

    input:
    file(fastqfile)
    file(metadata_clean)

    output:
    path("assembled/*.fastq.gz"), emit: assembled
    path("pear_counts.txt"), emit: counts
    path("mapping.txt"), emit: qiime_metadata

    script:
    """
    nice -${params.niceness} python3.11 $projectDir/bin/python/run_pear.py -m ${metadata_clean} --cpu ${params.cpus}
    """
}
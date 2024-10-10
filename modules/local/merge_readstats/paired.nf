process merge_readstats_paired  {
    // insert cutadapt and pear stats in dada2 stats
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(dada2_stats)
    file(cutadapt_stats)
    file(pear_stats)

    output:
    path("read_stats.tsv"), emit: read_stats
    
    script:
    """
    python3.11 $projectDir/bin/python/merge_readstats.py -d ${dada2_stats} -c ${cutadapt_stats} -p ${pear_stats} -s ${params.read_end}
    """
}
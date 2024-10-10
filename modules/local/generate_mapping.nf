process generate_mapping {
    // generate mapping file, including the paths required for Qiime import
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir params.outdir, mode: 'copy'
    
    input:
    file(reads)

    output:
    path("metadata.tsv"), emit: metadata_clean
    
    script: 
    """
    python3.11 $projectDir/bin/python/dummy_metadata.py -s ${params.read_end} -o metadata.tsv
    """
}
process generate_mapping {    
    input:
    file(reads)

    output:
    path("metadata.tsv"), emit: metadata_clean

    publishDir params.outdir, mode: 'copy'
    
    script: 
    """
    python3.11 $projectDir/bin/python/dummy_metadata.py -s ${params.read_end} -o metadata.tsv
    """
}
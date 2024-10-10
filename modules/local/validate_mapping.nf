process validate_mapping {
// validate existing mapping file
    container "$projectDir/containers/singularity/pyrrr.sif"

    input:
    file(metadata)
    file(filename)

    output:
    path("metadata_clean.tsv"), emit: metadata_clean

    script:
    """
    python3.11 $projectDir/bin/python/validate_mapping.py -i ${metadata} -s ${params.read_end}
    """
}
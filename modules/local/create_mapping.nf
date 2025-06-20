process CREATE_MAPPING {
    tag "$meta.id"
    label 'process_single'

    input:
        tuple val(meta), path(reads)
    
    output:
        tuple val(meta), path("mapping.csv"),   emit: dada2_mapping_file
    
    script:
        def args = task.ext.args ?: ''
        meta.id = 'batch'
        """
        echo "sample-id,absolute-filepath" > mapping.csv
        echo "batch,\"$(realpath "${reads[0]}")\"" >> mapping.csv
        """

    stub:
        def args = task.ext.args ?: ''
        meta.id = 'batch'
        """
        touch mapping.csv
        """
}


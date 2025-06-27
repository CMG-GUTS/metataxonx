process CREATE_MAPPING {
    label 'process_single'

    input:
    val(mapping_rows) // A list of lists!

    output:
    path("mapping.csv")             , emit: mapping

    script:
    def rows_str = mapping_rows.collect { row -> "${row[0]},${row[1]}" }.join('\n')
    """
    echo "sample_id,absolute-filepath" > mapping.csv
    echo "${rows_str}" >> mapping.csv
    """

    stub:
    """
    touch mapping.csv
    """
}

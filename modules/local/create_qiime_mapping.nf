process CREATE_QIIME_MAPPING {
    label 'process_single'

    input:
    val(mapping_rows) // A list of lists!

    output:
    path("mapping.tsv")             , emit: mapping

    script:
    def rows_str = mapping_rows.collect { row -> "${row[0]}\t${row[1]}" }.join('\n')
    """
    echo "sample-id\tabsolute-filepath" > mapping.tsv
    echo "${rows_str}" >> mapping.tsv
    """

    stub:
    """
    touch mapping.tsv
    """
}

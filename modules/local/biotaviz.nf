process BIOTAVIZ {
    label 'process_single'

    input:
    path(biomfile)

    output:
    path "biotaviz_clean_absolute.txt"          , emit: biotaviz_clean
    path "biotaviz_clean_relative.txt"          , emit: biotaviz_relative
    path "asv_table_with_taxonomy.txt"          , emit: asv_table_with_taxonomy_test
    path "versions.yml"                         , emit: versions

    script:
    """
    python3.11 $projectDir/bin/python/biom2biotaviz.py \\
        -i ${biomfile} \\
        -o foo.txt
    python3.11 $projectDir/bin/python/clean_biom_txt.py \\
        -i  ${biomfile}.txt \\
        -o biom_clean_absolute.txt
    python3.11 $projectDir/bin/python/biom2biotaviz.py \\
        -i biom_clean_absolute.txt \\
        -o biotaviz_clean_absolute.txt \\
        -t
    python3.11 $projectDir/bin/python/Biotaviz_counts_to_abundance.py \\
        -i biotaviz_clean_absolute.txt \\
        -o biotaviz_clean_relative.txt \\
        -r "${params.root_taxon}"
    cp biom_clean_absolute.txt asv_table_with_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3.11 --version)
    END_VERSIONS

    sed -i.bak -E '
    /^ *python:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    ' versions.yml
    """

    stub:
    """
    touch biotaviz_clean_absolute.txt
    touch biotaviz_clean_relative.txt
    touch asv_table_with_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3.11 --version)
    END_VERSIONS
    """
}
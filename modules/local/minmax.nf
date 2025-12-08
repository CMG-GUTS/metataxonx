process MINMAX {
    label 'process_single'

    input:
    path biom_file

    output:
    path "maxcount.txt"     , emit: maxcount
    path "mincount.txt"     , emit: mincount
    path "versions.yml"     , emit: versions

    script:
    """
    python3.11 << 'EOF'

    from biom import load_table

    biom_table = load_table("${biom_file}")
    sample_sums = biom_table.matrix_data.sum(axis=0)

    with open('maxcount.txt', 'w') as max_file:
        print(int(sample_sums.max()), file = max_file)

    with open('mincount.txt', 'w') as min_file:
        print(int(sample_sums.min()), file = min_file)
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3.11 --version 2>&1 | sed -e "s/Python //g")
    END_VERSIONS
    """

    stub:
    """
    touch maxcount.txt
    touch mincount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub-version
    END_VERSIONS
    """
}
process MINMAX {
    label 'process_single'

    input:
    path(featuretable)

    output:
    path "maxcount.txt"     , emit: maxcount
    path "mincount.txt"     , emit: mincount
    path "versions.yml"     , emit: versions

    script:
    """
    biom convert -i ${featuretable} -o ${featuretable}.txt --to-tsv
    python3.11 $projectDir/bin/python/table_minmax.py ${featuretable}.txt maximum > maxcount.txt
    python3.11 $projectDir/bin/python/table_minmax.py ${featuretable}.txt minimum > mincount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version)
        biom: \$(biom --version)
    END_VERSIONS
    """

    stub:
    """
    touch maxcount.txt
    touch mincount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version)
        biom: \$(biom --version)
    END_VERSIONS
    """
}
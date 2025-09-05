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
        python: \$(python3.11 --version)
        biom: \$(biom --version)
    END_VERSIONS

    sed -i.bak -E '
    /^ *python:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    /^ *biom:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    ' versions.yml
    """

    stub:
    """
    touch maxcount.txt
    touch mincount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3.11 --version)
        biom: \$(biom --version)
    END_VERSIONS
    """
}
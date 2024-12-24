process minmax {
    input:
    file(featuretable)

    output:
    path("maxcount.txt"), emit: maxcount
    path("mincount.txt"), emit: mincount

    script:
    """
    biom convert -i ${featuretable} -o ${featuretable}.txt --to-tsv
    python3.11 $projectDir/bin/python/table_minmax.py ${featuretable}.txt maximum > maxcount.txt
    python3.11 $projectDir/bin/python/table_minmax.py ${featuretable}.txt minimum > mincount.txt
    """
}
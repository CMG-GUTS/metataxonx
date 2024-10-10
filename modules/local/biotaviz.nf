process biom_to_biotaviz {
// generate biotaviz from biom-with-taxonomy
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir params.outdir, mode: 'copy'

    input:
    file(biomfile)

    output:
    path("biotaviz_clean_absolute.txt")
    path("biotaviz_clean_relative.txt"), emit: biotaviz
    path("asv_table_with_taxonomy.txt")

    script:
    """
    python3.11 $projectDir/bin/python/biom2biotaviz.py \
        -i ${biomfile} \
        -o foo.txt
    python3.11 $projectDir/bin/python/clean_biom_txt.py \
        -i  ${biomfile}.txt \
        -o biom_clean_absolute.txt
    python3.11 $projectDir/bin/python/biom2biotaviz.py \
        -i biom_clean_absolute.txt \
        -o biotaviz_clean_absolute.txt \
        -t
    python3.11 $projectDir/bin/python/Biotaviz_counts_to_abundance.py \
        -i biotaviz_clean_absolute.txt \
        -o biotaviz_clean_relative.txt \
        -r "${params.root_taxon}"
    cp biom_clean_absolute.txt asv_table_with_taxonomy.txt
    """
}
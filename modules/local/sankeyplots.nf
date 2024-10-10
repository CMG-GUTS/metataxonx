process sankeyplots  {
    // merge alpha diversity metrics
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir "${params.outdir}/sankeyplots", mode: 'copy'

    input:
    file(mapping)
    file(biotaviz)

    output:
    path("*.html")
    path("*.png"), emit: sankey_image

    script:
    """
    python3.11 $projectDir/bin/python/sankey-file-prep.py --taxa-filter 0.01 --sample-repeat false --combine-rankstat false -m ${mapping} -i ${biotaviz}
    Rscript $projectDir/bin/R/sankey-diagram-html-generator.R biotaviz_sankey_prepfile-AverageAllSamples.csv
    Rscript $projectDir/bin/R/sankey-diagram-png-generator.R biotaviz_sankey_prepfile-AverageAllSamples.html
    """
}
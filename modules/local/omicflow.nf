process OMICFLOW {
    label 'process_low'

    input:
    path(metadata_clean)
    path(biom_taxonomy)
    path(rooted_tree_newick)

    output:
    path "report.html"              , emit: report
    path "versions.yml"             , emit: versions

    script:
    def tree = rooted_tree_newick ? ",treeData = '${rooted_tree_newick}'" : ''
    """
    cat << 'EOF' > omicflow.R
    library('OmicFlow')
    library('ggplot2')

    set.seed(999)
    data.table::setDTthreads(${task.cpus})
    taxa <- metagenomics\$new(
        metaData = "${metadata_clean}",
        biomData = "${biom_taxonomy}"
        ${tree}
    )

    taxa\$feature_subset(Kingdom == "Bacteria")
    taxa\$normalize()
    taxa\$autoFlow(
        normalize = FALSE,
        threads = ${task.cpus}
    )
    EOF

    Rscript omicflow.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)
        OmicFlow: \$(Rscript -e 'cat(as.character(packageVersion("OmicFlow")))')
    END_VERSIONS

    # Rewrite the R version
    sed -i.bak -E '
    /^ *R:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    ' versions.yml

    """

    stub:
    """
    touch report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub-version
        OmicFlow: stub-version
    END_VERSIONS
    """
}
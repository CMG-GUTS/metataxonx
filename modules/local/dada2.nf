process DADA2 {
    label 'process_high'

    input:
    path(mapping)

    output:
    path("seq-tab.tsv")             , emit: seq_table
    path("denoising-stats.tsv")     , emit: denoising_stats
    path("rep-seqs.fna")            , emit: rep_seqs_fasta
    path("errProfile*.png")         , optional:true, emit: dada2_errors
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def novaseq = params.novaseq ? '--novaseq' : ''
    def denoise = params.bypass_denoise ? '--skip-denoise' : ''
    """
    parallel_dada2 \\
        --metadata ${mapping} \\
        --cpus ${task.cpus} \\
        ${novaseq} \\
        ${denoise} \\
        ${args} \\
        > dada_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)
        dada2: \$(Rscript -e 'cat(as.character(packageVersion("dada2")))')
    END_VERSIONS

    # Rewrite the R version
    sed -i.bak -E '
    /^ *R:/ s/(: *).*\\b([0-9]+\\.[0-9]+\\.[0-9]+)\\b.*/\\1 \\2/
    ' versions.yml

    # Rename and resize image for multiqc
    shopt -s nullglob
    files=(errProfile_*)
    if (( \${#files[@]} )); then
        rename s'/.png/_mqc.png/' *errProfile_*
        mogrify -resize 800x800 *.png
    fi
    shopt -u nullglob

    """
    
    stub:
    def args = task.ext.args ?: ''
    def novaseq = params.novaseq == true ? "--novaseq" : ""
    """
    touch seq-tab.tsv
    touch denoising-stats.tsv
    touch rep-seqs.fna
    touch errProfile.png
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub-version
        dada2: stub-version
    END_VERSIONS
    """
}
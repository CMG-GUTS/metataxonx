process cutadapt {
    input:
    file(filename)
    file(metadata_clean)
    
    output:
    path("trimmed/*"), emit: trimmed
    path("cutadapt_counts.txt"), emit: counts
    path("mapping.txt"), emit: qiime_metadata
    
    script:
    def fwd_primer = params.fwd_primer == false ? "$projectDir/assets/primers.fasta" : params.fwd_primer
    def rev_primer = params.rev_primer == false ? "$projectDir/assets/primers.fasta" : params.rev_primer
    def fwd_adapter = params.fwd_adapter == false ? "$projectDir/assets/adapters.fasta" : params.fwd_adapter
    def rev_adapter = params.rev_adapter == false ? "$projectDir/assets/adapters.fasta" : params.rev_adapter

    """
    python3.11 $projectDir/bin/python/run_cutadapt.py \
    -fp ${fwd_primer} -rp ${rev_primer} \
    -fa ${fwd_adapter} -ra ${rev_adapter} \
    -c ${params.cpus} -m ${metadata_clean} -s ${params.read_end}
    """
}
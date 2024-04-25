log.info "Microbiomics - N F"
log.info "====================================="
log.info "reads                  : ${params.reads}"
log.info "cpus                   : ${params.cpus}"
log.info "output                 : ${params.outdir}"
log.info "metadata file          : ${params.metadata}"
log.info "taxonomy classifier    : ${params.classifier}"
log.info "root taxon             : ${params.root_taxon}"
log.info "run_pear               : ${params.run_pear}"
log.info "nanopore               : ${params.nanopore}"
log.info "niceness               : ${params.niceness}"
log.info "forward primer         : ${params.fwd_primer}"
log.info "reverse primer         : ${params.rev_primer}"
log.info "====================================="
log.info "\n"

/*
Qiim2 workflow
Loosly based on CMBI - RadboudUMC DSL1 code
*/

// add pear stats to dada2 stats file
// add cutadapt stats to stats file

process generate_mapping {
    // generate mapping file, including the paths required for Qiime import
    container "script_dependencies:v1.0"

    publishDir params.outdir, mode: 'copy'
    
    input:
    file(reads)

    output:
    path("mapping.txt")
    
    script:
    if (params.nanopore == "yes" || params.run_pear == "yes") {
	singleton = "yes"
    }
    else {
	singleton = "no"
    }    
    """
    python3 $projectDir/bin/python/dummy_metadata.py ${params.reads} -o mapping.txt -n -p "\\\$PWD" -a ${singleton}
    """
}

process cutadapt {
    // run paired cutadapt
    container "script_dependencies:v1.0"
    
    input:
    file(filename)
    
    output:
    path("trimmed/*"), emit: trimmed
    path("cutadapt_counts.txt"), emit: counts
    
    script:
    """
    python3 $projectDir/bin/python/run_cutadapt.py -f ${params.fwd_primer} -r ${params.rev_primer} -c ${params.cpus} *.fastq.gz
    """
}


process validate_mapping {
// validate existing mapping file
    container "script_dependencies:v1.0"

    input:
    file(metadata)
    file(filename)

    output:
    path("mapping_for_nextflow.txt")

    script:
    """
    python3 $projectDir/bin/python/validate_mapping.py -i ${metadata} -d
    """
}

process fastqc {
// run fastqc
    container "staphb/fastqc:0.11.9"

    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    fastqc ${fastqfile}
    """
}

process multiqc {
// run multiqc
    container "script_dependencies:v1.0"
    publishDir params.outdir, mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc ./
    """
}

process fastqc_pear {
// run fastqc
    container "staphb/fastqc:0.11.9"

    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    fastqc ${fastqfile}
    """
}

process multiqc_pear {
// run multiqc
    container "script_dependencies:v1.0"
    publishDir "${params.outdir}/multiqc_pear", mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc ./
    """
}

process pear {
    // run Paired-End reAd mergeR
    container "qiimehelpers:v1.2"

    input:
    file(fastqfile)

    output:
    path("assembled/*.fastq.gz"), emit: assembled
    path("pear_counts.txt"), emit: counts

    script:
    """
    nice -${params.niceness} JOS_run_pear.py
    """
}

process qiime_import {
    // import sequences from fastq en metadata file
    container "quay.io/qiime2/core:2020.8"

    publishDir "${params.outdir}/qiime_artefacts/", mode: "copy"

    input:
    file(metadata)
    file(fastqfile)
    val(sampletype)
    val(inputformat)

    output:
    file("demux.qza")

    script:
    """
    qiime tools import \
         --type '${sampletype}' \
         --input-path ${metadata} \
         --output-path demux.qza \
         --input-format ${inputformat}
    """
}

process dada2_denoise {
    // denoise with dada2
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza}", mode: "copy"
    publishDir params.outdir, pattern: "*.{txt, fasta, biom}", mode: 'copy'

    input:
    file(artifact)

    output:
    path('dada_report.txt')
    path('representative_sequences.fasta')
    path('denoising_stats.txt'), emit: dada2stats
    path("table.qza"), emit: asv_qza
    path("rep-seqs.qza"), emit: repseqs_qza
    path("denoising-stats.qza"), emit: dada2stats_qza
    path('asv_table_no_taxonomy.biom'), emit: biom_no_taxonomy

    script:
    if (params.run_pear == "yes" || params.nanopore == "yes") {
        dadacall = "nice -${params.niceness} qiime dada2 denoise-single \
            --i-demultiplexed-seqs ${artifact} \
            --p-trim-left 0 \
            --p-n-threads 0 \
            --p-chimera-method \'consensus\' \
            --p-trunc-len 0 \
            --p-trunc-q 2 \
            --p-max-ee 3 \
            --p-min-fold-parent-over-abundance 2 \
            --o-table table.qza \
            --o-representative-sequences rep-seqs.qza \
            --o-denoising-stats denoising-stats.qza \
            --verbose \
            > dada_report.txt"
    }
    else {
        dadacall ="nice -${params.niceness} qiime dada2 denoise-paired \
            --i-demultiplexed-seqs ${artifact} \
            --p-trim-left-f 0 \
            --p-trim-left-r 0 \
            --p-n-threads 0 \
            --p-chimera-method \'consensus\' \
            --p-trunc-len-f 0 \
            --p-trunc-len-r 0 \
            --p-trunc-q 2 \
            --p-max-ee-f 2 \
            --p-max-ee-r 3 \
            --p-min-fold-parent-over-abundance 2 \
            --o-table table.qza \
            --o-representative-sequences rep-seqs.qza \
            --o-denoising-stats denoising-stats.qza \
            --verbose \
            > dada_report.txt"
    }

    """
    # run dada2_denoise
    ${dadacall}

    # convert absolute (raw) count table to biom format
    qiime tools export \
        --input-path table.qza \
        --output-path ./
    mv feature-table.biom asv_table_no_taxonomy.biom

    #denoising stats
    qiime tools export \
        --input-path denoising-stats.qza \
        --output-path ./
    mv stats.tsv denoising_stats.txt

    #create representative sequences fasta file in qiime artifact
    qiime feature-table tabulate-seqs \
        --i-data rep-seqs.qza \
        --o-visualization rep-seqs.qzv

    #convert artifact to .tsv file
    qiime tools export \
        --input-path rep-seqs.qzv \
        --output-path ./
    mv sequences.fasta representative_sequences.fasta
    """
}

process assign_taxonomy {
    // assign taxonomy using sklearn qiime2
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.tsv", mode: 'copy'

    input:
    file(classifier)
    file(sequences)

    output:
    path("taxonomy.tsv"), emit: taxonomy_tsv
    path("taxonomy_sklearn.qza"), emit: taxonomy_qza

    script:
    """
    nice -${params.niceness} qiime feature-classifier classify-sklearn \
         --i-classifier ${classifier} \
         --i-reads ${sequences} \
         --o-classification taxonomy_sklearn.qza \
         --p-n-jobs ${params.cpus}

    qiime metadata tabulate \
         --m-input-file taxonomy_sklearn.qza \
         --o-visualization taxonomy_sklearn.qzv

    qiime tools export --input-path taxonomy_sklearn.qza  \
         --output-path taxonomy

    qiime tools export --input-path taxonomy_sklearn.qzv  \
         --output-path taxonomy

    cat taxonomy/taxonomy.tsv | \
    sed 's/Feature ID/#OTUID/' | \
    sed 's/Taxon/taxonomy/' | \
    sed 's/Consensus/consensus/' > taxonomy/taxonomy_relabeled.tsv
    mv taxonomy/taxonomy_relabeled.tsv ./taxonomy.tsv
    """
}

process combine_taxonomy_biom {
// add taxonomy info to biom file
    container "quay.io/qiime2/core:2020.8"

    publishDir params.outdir, mode: 'copy'

    input:
    path(biomfile)
    path(taxonomyfile)

    output:
    path("biom_with_taxonomy.biom")

    script:
    """
    biom add-metadata \
        -i ${biomfile} \
        -o biom_with_taxonomy.biom \
        --observation-metadata-fp ${taxonomyfile} \
        --sc-separated taxonomy
    """
}

process biom_to_biotaviz {
// generate biotaviz from biom-with-taxonomy
    container "script_dependencies:v1.0"

    publishDir params.outdir, mode: 'copy'

    input:
    file(biomfile)

    output:
    path("biotaviz_clean_absolute.txt")
    path("biotaviz_clean_relative.txt"), emit: biotaviz
    path("asv_table_with_taxonomy.txt")

    script:
    """
    python3 $projectDir/bin/python/biom2biotaviz.py \
        -i ${biomfile} \
        -o foo.txt
    python3 $projectDir/bin/python/clean_biom_txt.py \
        -i  ${biomfile}.txt \
        -o biom_clean_absolute.txt
    python3 $projectDir/bin/python/biom2biotaviz.py \
        -i biom_clean_absolute.txt \
        -o biotaviz_clean_absolute.txt \
        -t
    python3 $projectDir/bin/python/Biotaviz_counts_to_abundance.py \
        -i biotaviz_clean_absolute.txt \
        -o biotaviz_clean_relative.txt \
        -r "${params.root_taxon}"
    cp biom_clean_absolute.txt asv_table_with_taxonomy.txt
    """
}

process phylogeny {
    // align ASVs
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.newick", mode: 'copy'

    input:
    file(sequences)

    output:
    path("rooted_tree.newick")
    path("rooted-tree.qza"), emit: rooted_tree_qza

    script:
    """
    nice -${params.niceness} qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${sequences} \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --p-n-threads ${params.cpus} \
        --o-rooted-tree rooted-tree.qza

    qiime tools export \
        --input-path rooted-tree.qza  \
		--output-path ./

	mv tree.nwk rooted_tree.newick
    """
}

process alpha_rarefaction {
// alpha rarefaction
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.csv", mode: 'copy'
    publishDir params.outdir, pattern: "alpha_rarefaction/*", mode: 'copy'

    output:
    path("alpha_rarefaction.qzv")
    path("alpha_rarefaction/*")

    input:
    file(table)
    file(phylogeny)
    file(maxcount)

    output:
    file("alpha_rarefaction.qzv")
    file("alpha_rarefaction/*.csv")

    script:
    """
    nice -${params.niceness} qiime diversity alpha-rarefaction \
        --i-table ${table} \
        --i-phylogeny ${phylogeny} \
        --p-max-depth "\$(<${maxcount})" \
        --o-visualization alpha_rarefaction.qzv

    qiime tools export \
        --input-path alpha_rarefaction.qzv \
	--output-path alpha_rarefaction
    """
}

process minmax {
    // get the minimum and maximum readcounts
    container "script_dependencies:v1.0"

    input:
    file(featuretable)

    output:
    path("maxcount.txt"), emit: maxcount
    path("mincount.txt"), emit: mincount

    script:
    """
    biom convert -i ${featuretable} -o ${featuretable}.txt --to-tsv
    python3 $projectDir/bin/python/table_minmax.py ${featuretable}.txt maximum > maxcount.txt
    python3 $projectDir/bin/python/table_minmax.py ${featuretable}.txt minimum > mincount.txt
    """
}

process core_diversity {
    // div diversity metrics (alpha + beta)
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir "${params.outdir}/distance_matrices/", pattern: "*.txt", mode: "copy"
    
    input:
    file(metadata)
    file(tree)
    file(table)
    file(mincount)

    output:
    path("diversity_core/*vector.qza"), emit: alpha
    path("*.txt")

    script:
    """
    nice -${params.niceness} qiime diversity core-metrics-phylogenetic \
		--m-metadata-file ${metadata} \
		--i-phylogeny ${tree} \
		--i-table ${table} \
        --p-sampling-depth "\$(<${mincount})" \
		--output-dir diversity_core \
		--p-n-jobs-or-threads ${params.cpus} \
		--quiet

	for i in diversity_core/*distance*; do \
            qiime tools export \
                --input-path \$i \
		        --output-path \$i.data; \
		    mv \$i.data/distance-matrix.tsv ./\$i.txt; done

    mv diversity_core/*.txt .
    """
}

process merge_alpha_diversity {
    // merge alpha diversity metrics
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir params.outdir, mode: 'copy'

    input:
    file(alpha)

    output:
    path("alpha_diversity.txt")

    script:
    """
    qiime metadata tabulate \
        --m-input-file *vector.qza \
        --o-visualization combined-alpha-metadata.qzv

    qiime tools export \
        --input-path combined-alpha-metadata.qzv \
		--output-path data

	cp data/metadata.tsv alpha_diversity.txt
    """
}

process beta_rarefaction {
    // beta rarefaction & tree building
    container "quay.io/qiime2/core:2020.8"
    containerOptions "--cpus ${params.cpus}"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.newick", mode: 'copy'
    publishDir params.outdir, pattern: "beta_rarefaction/*", mode: 'copy'

    input:
    file(table)
    file(metadata)
    file(mincount)
    file(tree)

    output:
    file("beta_rarefaction.qzv")
    file("weighted_unifrac_tree.newick")
    path("beta_rarefaction/*")

    script:
    """

    nice -${params.niceness} qiime diversity beta-rarefaction \
        --i-table ${table} \
        --i-phylogeny ${tree} \
        --p-metric weighted_unifrac \
        --p-clustering-method upgma \
        --m-metadata-file ${metadata} \
        --p-sampling-depth "\$(<${mincount})"  \
        --o-visualization beta_rarefaction.qzv

    qiime tools export \
        --input-path beta_rarefaction.qzv \
		--output-path beta_rarefaction

    cp beta_rarefaction/sample-clustering-upgma.tre ./weighted_unifrac_tree.newick
    """
}


process sankeyplots  {
    // merge alpha diversity metrics
    container "script_dependencies:v1.0"

    publishDir "${params.outdir}/sankeyplots", mode: 'copy'

    input:
    file(mapping)
    file(biotaviz)

    output:
    path("*.html")
    path("*.png")

    script:
    """
    python3 $projectDir/bin/python/Biotaviz2sankey.py -m ${mapping} -i ${biotaviz}
    """
}

process merge_readstats_both  {
    // insert cutadapt and pear stats in dada2 stats
    container "script_dependencies:v1.0"

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(dada2_stats)
    file(cutadapt_stats)
    file(pear_stats)

    output:
    path("read_stats.txt")
    
    script:
    """
    python3 $projectDir/bin/python/merge_readstats.py -d ${dada2_stats} -c ${cutadapt_stats} -p ${pear_stats} -o read_stats.txt
    """
}


process merge_readstats_pear  {
    // insert pear stats in dada2 stats
    container "script_dependencies:v1.0"

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(dada2_stats)
    file(pear_stats)

    output:
    path("read_stats.txt")
    
    script:
    """
    python3 $projectDir/bin/python/merge_readstats.py -d ${dada2_stats} -p ${pear_stats} -o read_stats.txt
    """
}


process merge_readstats_cutadapt  {
    // insert cutadapt stats in dada2 stats
    container "script_dependencies:v1.0"

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(dada2_stats)
    file(cutadapt_stats)

    output:
    path("read_stats.txt")
    
    script:
    """
    python3 $projectDir/bin/python/merge_readstats.py -d ${dada2_stats} -c ${cutadapt_stats} -o read_stats.txt
    """
}



workflow {

    sample_ch = Channel.fromPath("${params.reads}",  checkIfExists:true)
    classifier_ch = Channel.fromPath(params.classifier)

    // mapping file stuff
    if (params.metadata) { // was a mapping file provided?
	    mapping_ch = Channel.fromPath("${params.metadata}")
	    validate_mapping(mapping_ch, sample_ch.collect())
	    mapping_ch = validate_mapping.out
    }
    else { // if not, generate one
	    generate_mapping(sample_ch.collect())
	    mapping_ch = generate_mapping.out
    }

    // QC stuff
    fastqc(sample_ch.collect())
    multiqc(fastqc.out.collect())

    // cutadapt
    if (params.fwd_primer != "" && params.rev_primer != "") {
	cutadapt(sample_ch.collect())
	sample_ch = cutadapt.out.trimmed
    }
    
    // Qiime stuff
    // assembly (if any), import, and DADA2 are pairedness-dependent
    if (params.run_pear == "yes") { 
        sampletype = "'SampleData[SequencesWithQuality]'"
        inputformat = 'SingleEndFastqManifestPhred33V2'
        pear(sample_ch.collect())
        sample_ch = pear.out.assembled
	fastqc_pear(sample_ch)
	multiqc_pear(fastqc_pear.out.collect())
    }
    else if (params.nanopore == "yes"){
        sampletype = "'SampleData[SequencesWithQuality]'"
        inputformat = 'SingleEndFastqManifestPhred33V2'	
    }
    else {
        sampletype = "'SampleData[PairedEndSequencesWithQuality]'"
        inputformat = 'PairedEndFastqManifestPhred33V2'
    }
    
    qiime_import(mapping_ch, sample_ch.collect(), sampletype, inputformat)
    dada2_denoise(qiime_import.out)

    // below bit: conditional readstats merger. A nextflow conditional channel input option would be nice.
    if (params.fwd_primer != "" && params.run_pear == "yes"){
	merge_readstats_both(dada2_denoise.out.dada2stats, cutadapt.out.counts, pear.out.counts)
    }
    if (params.fwd_primer != "" && params.run_pear == "no"){
	merge_readstats_cutadapt(dada2_denoise.out.dada2stats, cutadapt.out.counts)
    }
    if (params.fwd_primer == "" && params.run_pear == "yes"){
	merge_readstats_pear(dada2_denoise.out.dada2stats, pear.out.counts)
    }
    
    assign_taxonomy(classifier_ch, dada2_denoise.out.repseqs_qza)
    combine_taxonomy_biom(dada2_denoise.out.biom_no_taxonomy, assign_taxonomy.out.taxonomy_tsv)
    phylogeny(dada2_denoise.out.repseqs_qza)
    biom_to_biotaviz(combine_taxonomy_biom.out)
    minmax(dada2_denoise.out.biom_no_taxonomy)
    alpha_rarefaction(dada2_denoise.out.asv_qza, phylogeny.out.rooted_tree_qza, minmax.out.maxcount)
    core_diversity(mapping_ch, phylogeny.out.rooted_tree_qza, dada2_denoise.out.asv_qza, minmax.out.mincount)
    merge_alpha_diversity(core_diversity.out.alpha.collect())
    beta_rarefaction(dada2_denoise.out.asv_qza, mapping_ch, minmax.out.mincount, phylogeny.out.rooted_tree_qza)

    // post-qiime visualization and statistics
    sankeyplots(mapping_ch, biom_to_biotaviz.out.biotaviz)
}

log.info "Microbiomics - N F"
log.info "====================================="
log.info "reads                  : ${params.reads}"
log.info "cpus                   : ${params.cpus}"
log.info "output                 : ${params.outdir}"
log.info "metadata file          : ${params.metadata}"
log.info "taxonomy classifier    : ${params.classifier}"
log.info "root taxon             : ${params.root_taxon}"
log.info "batch_size             : ${params.batch_size}"
log.info "nanopore               : ${params.nanopore}"
log.info "niceness               : ${params.niceness}"
log.info "Sequence reading       : ${params.seq_read}"
log.info "forward primer         : ${params.fwd_primer}"
log.info "reverse primer         : ${params.rev_primer}"
log.info "forward adapter        : ${params.fwd_adapter}"
log.info "reverse adapter        : ${params.rev_adapter}"
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
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir params.outdir, mode: 'copy'
    
    input:
    file(reads)

    output:
    path("metadata.tsv"), emit: metadata_clean
    
    script: 
    """
    python3.11 $projectDir/bin/python/dummy_metadata.py -s ${params.seq_read} -o metadata.tsv
    """
}

process cutadapt {
    // run paired cutadapt
    container "$projectDir/containers/singularity/pyrrr.sif"
    
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
    -c ${params.cpus} -m ${metadata_clean} -s ${params.seq_read}
    """
}


process validate_mapping {
// validate existing mapping file
    container "$projectDir/containers/singularity/pyrrr.sif"

    input:
    file(metadata)
    file(filename)

    output:
    path("metadata_clean.tsv"), emit: metadata_clean

    script:
    """
    python3.11 $projectDir/bin/python/validate_mapping.py -i ${metadata} -s ${params.seq_read}
    """
}

process fastqc_trimmed {
// run fastqc
    container "$projectDir/containers/singularity/fastqc.sif"

    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    ls * | parallel -j ${params.cpus} --verbose "fastqc {}"
    """
}

process multiqc_trimmed {
// run multiqc
    container "$projectDir/containers/singularity/pyrrr.sif"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc --force --interactive .
    """
}

process fastqc_untrimmed {
// run fastqc
    container "$projectDir/containers/singularity/fastqc.sif"

    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    ls * | parallel -j ${params.cpus} --verbose "fastqc {}"
    """
}

process multiqc_untrimmed {
// run multiqc
    container "$projectDir/containers/singularity/pyrrr.sif"
    publishDir "${params.outdir}/untrimmed", mode: 'copy'
    
    input:
    file(qcresults)

    output:
    path("multiqc*")

    script:
    """
    multiqc --force --interactive .
    """
}

process fastqc_pear {
// run fastqc
    container "$projectDir/containers/singularity/fastqc.sif"

    input:
    file(fastqfile)

    output:
    path("*.{zip,html}")

    script:
    """
    ls * | parallel -j ${params.cpus} --verbose "fastqc {}"
    """
}

process multiqc_pear {
// run multiqc
    container "$projectDir/containers/singularity/pyrrr.sif"
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
    container "$projectDir/containers/singularity/pyrrr.sif"

    input:
    file(fastqfile)
    file(metadata_clean)

    output:
    path("assembled/*.fastq.gz"), emit: assembled
    path("pear_counts.txt"), emit: counts
    path("mapping.txt"), emit: qiime_metadata

    script:
    """
    nice -${params.niceness} python3.11 $projectDir/bin/python/run_pear.py -m ${metadata_clean} --cpu ${params.cpus}
    """
}

process dada2_denoise {
    // denoise with dada2
    container "$projectDir/containers/singularity/parallel_dada2.sif"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{rds}", mode: "copy"

    input:
    file(mapping)

    output:
    path("seq-tab.tsv"), emit: seq_table
    path("denoising-stats.tsv"), emit: denoising_stats
    path("rep-seqs.fna"), emit: rep_seqs_fasta
    path("dada_report.txt")

    script:
    """
    Rscript $projectDir/bin/R/parallel_dada2.R --metadata ${mapping} --batch_n ${params.batch_size} --cpus ${params.cpus} > dada_report.txt
    """
}

process qiime_import_export {
    // import sequences from fastq and metadata file
    container "$projectDir/containers/singularity/qiime2.sif"

    publishDir "${params.outdir}/qiime_artefacts/", mode: "copy"

    input:
    file(seq_table)
    file(rep_seqs_fasta)

    output:
    // path('representative_sequences.fasta'), emit: sequences_fasta
    path("table.qza"), emit: asv_qza
    path('asv_table_no_taxonomy.biom'), emit: biom_no_taxonomy
    path('rep-seqs.qza'), emit: repseqs_qza

    script:
    """
    # Converts phyloseq biom to qiime compatible biom file
    biom convert -i ${seq_table} -o table.biom --to-hdf5

    # Convert RDS files from parallel_dada2 into qiime2 formats
    qiime tools import \
    --input-path table.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path table.qza    

    # convert absolute (raw) count table to biom format
    qiime tools export \
        --input-path table.qza \
        --output-path ./
    mv feature-table.biom asv_table_no_taxonomy.biom

    # representative sequences
    qiime tools import \
    --input-path ${rep_seqs_fasta} \
    --type 'FeatureData[Sequence]' \
    --output-path rep-seqs.qza
    """
}

process assign_taxonomy {
    // assign taxonomy using sklearn qiime2
    container "$projectDir/containers/singularity/qiime2.sif"

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
    container "$projectDir/containers/singularity/qiime2.sif"

    publishDir params.outdir, mode: 'copy'

    input:
    path(biomfile)
    path(taxonomyfile)

    output:
    path("biom_with_taxonomy.biom"), emit: biom_taxonomy

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

process phylogeny {
    // align ASVs
    container "$projectDir/containers/singularity/qiime2.sif"

    publishDir "${params.outdir}/qiime_artefacts/", pattern: "*.{qza, qzv}", mode: "copy"
    publishDir params.outdir, pattern: "*.newick", mode: 'copy'

    input:
    file(sequences)

    output:
    path("rooted_tree.newick"), emit: rooted_tree_newick
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
    container "$projectDir/containers/singularity/qiime2.sif"
    cpus = params.cpus

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
    container "$projectDir/containers/singularity/pyrrr.sif"

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

process core_diversity {
    // div diversity metrics (alpha + beta)
    container "$projectDir/containers/singularity/qiime2.sif"

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
    container "$projectDir/containers/singularity/qiime2.sif"
    cpus = params.cpus

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
    container "$projectDir/containers/singularity/qiime2.sif"
    cpus = params.cpus

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


process merge_readstats_paired  {
    // insert cutadapt and pear stats in dada2 stats
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(dada2_stats)
    file(cutadapt_stats)
    file(pear_stats)

    output:
    path("read_stats.tsv"), emit: read_stats
    
    script:
    """
    python3.11 $projectDir/bin/python/merge_readstats.py -d ${dada2_stats} -c ${cutadapt_stats} -p ${pear_stats} -s ${params.seq_read}
    """
}

process merge_readstats_single  {
    // insert cutadapt stats in dada2 stats
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file(dada2_stats)
    file(cutadapt_stats)

    output:
    path("read_stats.tsv"), emit: read_stats
    
    script:
    """
    python3.11 $projectDir/bin/python/merge_readstats.py -d ${dada2_stats} -c ${cutadapt_stats} -s ${params.seq_read}
    """
}

process omics_analysis {
    container "$projectDir/containers/singularity/pyrrr.sif"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(metadata_clean)
    file(biom_taxonomy)
    file(rooted_tree_newick)
    file(sequences_fasta)
    file(read_stats)
    file(sankey_image)

    output:
    file("report.html")

    script:
    """
    Rscript $projectDir/bin/R/OmicFlow/00_main.R \
        --metadata ${metadata_clean} \
        --biom ${biom_taxonomy} \
        --tree ${rooted_tree_newick} \
        --refseq ${sequences_fasta} \
        --outdir ${params.outdir}
    """
}

// importing java API utils
import java.nio.file.Files
import java.nio.file.Paths

workflow {
    try {
        // Implements automatic directory file fetching
        def path = Paths.get(params.reads).toAbsolutePath().toString()
        def basename = new File(path).getName()

        if (Files.isDirectory(Paths.get(path))) {
            sample_ch = Channel.fromPath("${params.reads}/*.fastq.gz",  checkIfExists:true)
        } else if (!basename.isEmpty()){
            sample_ch = Channel.fromPath("${params.reads}",  checkIfExists:true)
        } else {
            throw new IllegalArgumentException("No matching path found.")
        }
    } catch (Exception e) {
        println("Error: ${e.message}")
    }
    trimmed_ch = sample_ch
    classifier_ch = Channel.fromPath(params.classifier)

    // mapping file stuff
    if (params.metadata != "") { // was a mapping file provided?
        mapping_ch = Channel.fromPath("${params.metadata}")
	    validate_mapping(mapping_ch, sample_ch.collect())
	    mapping_ch = validate_mapping.out
    }
    else { // if not, generate one
	    generate_mapping(sample_ch.collect())
	    mapping_ch = generate_mapping.out
    }

    // cutadapt
    if (params.seq_read != "") {
        cutadapt(sample_ch.collect(), mapping_ch.metadata_clean)
        trimmed_ch = cutadapt.out.trimmed
        qiime_mapping = cutadapt.out.qiime_metadata        
    }
        
    // Quality control on both untrimmed and trimmed reads
    // combined_ch = sample_ch.toList().merge(trimmed_ch)

    // QC
    if (params.fwd_primer == "skip" && params.fwd_adapter == "skip") {
        fastqc_untrimmed(sample_ch.collect())
        multiqc_untrimmed(fastqc_untrimmed.out.collect())
    } else {
        fastqc_untrimmed(sample_ch.collect())
        multiqc_untrimmed(fastqc_untrimmed.out.collect())
        fastqc_trimmed(trimmed_ch.collect())
        multiqc_trimmed(fastqc_trimmed.out.collect())
    }    
    // Qiime stuff
    // assembly (if any), import, and DADA2 are pairedness-dependent
    if (params.seq_read == "paired") {
        pear(trimmed_ch.collect(), mapping_ch.metadata_clean)
        fastqc_pear(trimmed_ch.collect())
        multiqc_pear(fastqc_pear.out.collect())
        trimmed_ch = pear.out.assembled
        qiime_mapping = pear.out.qiime_metadata
    }
    
    dada2_denoise(qiime_mapping)
    qiime_import_export(dada2_denoise.out.seq_table, dada2_denoise.out.rep_seqs_fasta)

    // readstat merger
    if (params.seq_read == "paired"){
	    merge_readstats_paired(dada2_denoise.out.denoising_stats, cutadapt.out.counts, pear.out.counts)
        reads_stats_file = merge_readstats_paired.out.read_stats
    } else if (params.seq_read == "single") {
        merge_readstats_single(dada2_denoise.out.denoising_stats, cutadapt.out.counts)
        reads_stats_file = merge_readstats_single.out.read_stats
    }
    
    assign_taxonomy(classifier_ch, qiime_import_export.out.repseqs_qza)
    combine_taxonomy_biom(qiime_import_export.out.biom_no_taxonomy, assign_taxonomy.out.taxonomy_tsv)
    phylogeny(qiime_import_export.out.repseqs_qza)
    biom_to_biotaviz(combine_taxonomy_biom.out)
    minmax(qiime_import_export.out.biom_no_taxonomy)
    alpha_rarefaction(qiime_import_export.out.asv_qza, phylogeny.out.rooted_tree_qza, minmax.out.maxcount)
    core_diversity(mapping_ch.metadata_clean, phylogeny.out.rooted_tree_qza, qiime_import_export.out.asv_qza, minmax.out.mincount)
    merge_alpha_diversity(core_diversity.out.alpha.collect())
    beta_rarefaction(qiime_import_export.out.asv_qza, mapping_ch.metadata_clean, minmax.out.mincount, phylogeny.out.rooted_tree_qza)

    // post-qiime visualization and statistics
    sankeyplots(mapping_ch.metadata_clean, biom_to_biotaviz.out.biotaviz)
    omics_analysis(mapping_ch.metadata_clean, combine_taxonomy_biom.out.biom_taxonomy, phylogeny.out.rooted_tree_newick, dada2_denoise.out.rep_seqs_fasta, reads_stats_file, sankeyplots.out.sankey_image)
}

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
log.info "Sequence reading       : ${params.read_end}"
log.info "forward primer         : ${params.fwd_primer}"
log.info "reverse primer         : ${params.rev_primer}"
log.info "forward adapter        : ${params.fwd_adapter}"
log.info "reverse adapter        : ${params.rev_adapter}"
log.info "====================================="
log.info "\n"


/*
=========================================================================
IMPORT API utils
=========================================================================
*/

import java.nio.file.Files
import java.nio.file.Paths


/*
=========================================================================
IMPORT MODULES
=========================================================================
*/

// metadata validation
include { generate_mapping } from './modules/local/generate_mapping'
include { validate_mapping } from './modules/local/validate_mapping'

// QC modules, to be merged into a single process as python script
include { fastqc_pear } from './modules/local/fastqc/fastqc_pear'
include { fastqc_trimmed } from './modules/local/fastqc/fastqc_trimmed'
include { fastqc_untrimmed } from './modules/local/fastqc/fastqc_untrimmed'
include { multiqc_pear } from './modules/local/multiqc/multiqc_pear'
include { multiqc_trimmed } from './modules/local/multiqc/multiqc_trimmed'
include { multiqc_untrimmed } from './modules/local/multiqc/multiqc_untrimmed'

// Merge and trim/denoise modules
include { dada2_denoise } from './modules/local/dada2'
include { cutadapt } from './modules/local/cutadapt'
include { pear } from './modules/local/pear'

// readstats
include { merge_readstats_single } from './modules/local/merge_readstats/single'
include { merge_readstats_paired } from './modules/local/merge_readstats/paired'

// qiime2 modules
include { qiime_import_export } from './modules/local/import_export'
include { assign_taxonomy } from './modules/local/taxonomy'
include { combine_taxonomy_biom } from './modules/local/combine_taxonomy'
include { phylogeny } from './modules/local/phylogeny'

// alpha and beta diversity
include { minmax } from './modules/local/minmax'
include { alpha_rarefaction } from './modules/local/alpha_rarefaction'
include { merge_alpha_diversity } from './modules/local/merge_alpha_diversity'
include { beta_rarefaction } from './modules/local/beta_rarefaction'
include { core_diversity } from './modules/local/core_diversity'

// visualization modules
include { biom_to_biotaviz } from './modules/local/biotaviz'
include { sankeyplots } from './modules/local/sankeyplots'
include { omics_analysis } from './modules/local/omicsflow'
    
/*
=========================================================================
RUN MAIN WORKFLOW
=========================================================================
*/

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
    if (params.read_end != "") {
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
    if (params.read_end == "paired") {
        pear(trimmed_ch.collect(), mapping_ch.metadata_clean)
        fastqc_pear(trimmed_ch.collect())
        multiqc_pear(fastqc_pear.out.collect())
        trimmed_ch = pear.out.assembled
        qiime_mapping = pear.out.qiime_metadata
    }
    
    dada2_denoise(qiime_mapping)
    qiime_import_export(dada2_denoise.out.seq_table, dada2_denoise.out.rep_seqs_fasta)

    // readstat merger
    if (params.read_end == "paired"){
	    merge_readstats_paired(dada2_denoise.out.denoising_stats, cutadapt.out.counts, pear.out.counts)
        reads_stats_file = merge_readstats_paired.out.read_stats
    } else if (params.read_end == "single") {
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
    omics_analysis(mapping_ch.metadata_clean, combine_taxonomy_biom.out.biom_taxonomy, phylogeny.out.rooted_tree_newick, dada2_denoise.out.rep_seqs_fasta, reads_stats_file, sankeyplots.out.sankey_image, alpha_rarefaction.out.shannon_file)
}

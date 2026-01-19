/*

    PIPELINE INITIALISATION
    
*/

include { paramsSummaryMap          }   from    'plugin/nf-validation'
include { fromSamplesheet           }   from    'plugin/nf-validation'

include { metapipeLogo              }   from    '../nf-core/nf_pipeline_utils.nf'
include { completionEmail           }   from    '../nf-core/nf_pipeline_utils.nf'
include { completionSummary         }   from    '../nf-core/nf_pipeline_utils.nf'
include { dashedLine                }   from    '../nf-core/nf_pipeline_utils.nf'
include { imNotification            }   from    '../nf-core/nf_pipeline_utils.nf'
include { UTILS_NEXTFLOW_PIPELINE }     from    '../nf-core/utils_nextflow_pipeline.nf'
include { UTILS_NFVALIDATION_PLUGIN }   from    '../nf-core/utils_nfvalidation_plugin.nf'


/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        false
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = metapipeLogo(monochrome_logs)
    post_help_text = """
    =============================================
    ${workflow.manifest.name.toUpperCase()} v${workflow.manifest.version}
    =============================================

    Usage:
        nextflow run ${workflow.manifest.mainScript} \\
            -profile <docker/singularity> \\
            -work-dir <workdir> \\
            --input <samplesheet.csv> \\
            --outdir <outdir>

    Description:
        ${workflow.manifest.description ?: "A 16S metataxonomics pipeline"}

Data Input Options:
        --input                     Path to samplesheet CSV file
        --reads                     Path to input fastq files
        --outdir                    Path to output directory
        --singleEnd                 Using flag sets singleEnd to true

    Profiles:
        -profile                    Configuration profile to use. 
                                    Available: docker, singularity

    Key Pipeline Options:
        --read_merger_tool          Paired-end Merger tool to use (default: flash)
        --classifier_name           Name of classifier to be downloaded ('standard', 'diverse', 'human-stool', 'greengenes', default: 'standard')
        --classifier_path           Path to save and download the classifier (default: $projectDir/classifiers)
        --classifier_custom         Path to a custom classifier compatible to sklearn version 1.4.0 (default: null)
        --classifier_test_set       Path to a dummy file of sequences to check the custom classifier compatibility (default: $projectDir/assets/classifier_test_query.fna).
        --standard_adapters         Path to a collection of adapter sequences (default: $projectDir/assets/adapters.fasta)
        --standard_primers          Path to a collection of primer sequences (default: $projectDir/assets/primers.fasta)
        --custom_primer_1           Custom primer sequence forward (default: none)
        --custom_primer_2           Custom primer sequence reverse (default: none)
        --custom_adapter_1          Custom adapter sequence forward (default: none)
        --custom_adapter_2          Custom adapter sequence reverse (default: none)

    Software non-default settings
        --multiqc_config            Path to custom multiqc yaml file (default: '$projectDir/assets/multiqc_config.yaml')
        --seqkit_opt                Seqkit sub-sampling options (default: '-s100 -n 100000')
        --pear_opt                  Paired-end Read Merger (PEAR) options (default: '-y 25G -j 32 -q 30 -v 35 -p 0.0001')
        --flash_opt                 FLASH paired-end read merger options (default: '--min-overlap 10 --max-mismatch-densit 0.2')
        --cutadapt_opt              Cutadapt Read Trimming options (default: '--discard-untrimmed --minimum-length 20')
        --dada2_opt                 DADA2 Denoising options (default: '--p-trunc-q 2 --p-max-ee 6 --p-min-fold-parent-over-abundance 2 --p-chimera-method consensus --batch 500')
        --novaseq                   DADA2 NovaSeq denoising (default: false)
        --root_taxon                Biotaviz feature subset (default: 'domain - Bacteria')

    Process Bypass Options:
        --bypass_trim               Skip read trimming (default: false)
        --bypass_denoise            Skip DADA2 denoising (default: false)
        --bypass_post_analysis      Skip Alpha and Beta diversity analysis (default: false)
        --bypass_report             Skip MultiQC and OmicFlow reports (default: false)

    File Saving Options:
        --save_trim_reads           Save trimmed reads (default: false)
        --save_sampled_reads        Save interleaved reads (default: false)
        --save_merged_reads         Save paired-end merged reads (default: true)

        --save_denoise_stats        Save DADA2 denoising stats (default: true)
        --save_denoise_sequences    Save representative sequences (default: true)
        --save_denoise_table        Save abundance table (default: true)

        --save_biom_files           Save BIOM format file (default: true)
        --save_tree_files           Save NEWICK format tree from fasttree (default: true)
        --save_alpha_div_files      Save alpha diversity metrics (default: true)
        --save_beta_div_files       Save beta diversity metrics (default: true)
        --save_qiime_artifacts      Save qiime2 artifacts (default: true)

        --save_biotaviz_files       Save BiotaViz format (default: false)
        --save_sankey_plot          Save sankey plot of taxa (default: true)
        --save_final_report         Save report file (default: true)

    Resources Options:
        --process_low_cpu           Number of cores to allocate for process with low cpu requirement (default: 4)
        --process_med_cpu           Number of cores to allocate for process with medium cpu requirement (default: 8)
        --process_high_cpu          Number of cores to allocate for process with high cpu requirement (default: 16)
        --cpus                      Maximum number of cores to allocate for the global pipeline scope (default: 32)

    For advanced resource customization, see the configuration file:
        ${projectDir}/conf/base.config
    or supply your own config file with the -c option:
        -c path/to/your.config

    Other Options:
        --help          Show this help message and exit
        --version       Show the pipeline version and exit

    For additional customization, see the configuration file: 
        ${projectDir}/nextflow.config

    Example Command:
        nextflow run ${workflow.manifest.mainScript} -profile docker --reads '*_{1,2}.fastq.gz' --outdir results
        nextflow run ${workflow.manifest.mainScript} -profile singularity --input samplesheet.csv --outdir results

    For more information, visit: ${workflow.manifest.homePage}
    """.stripIndent()
    
    workflow_command = "nextflow run ${workflow.manifest.mainScript} -profile <docker/singularity/test> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Custom validation for pipeline parameters
    //
    // validateParameters()

    emit:
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        if (email_on_fail) {
            // placeholder
        }
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
// def validateParameters() {
//     if ( !params.reads && !params.input) {
//         error("Missing reads and input declaration, one is required.")
//     }

//     if ( !params.outdir ) {
//         error("Missing output directory declaration: --outdir` is required.")
//     }
// }


//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
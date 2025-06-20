/*

    COMPUTATION OF BETA & ALPHA DIVERSITY

*/
include { BIOTAVIZ } from '../../modules/local/biotaviz.nf'
include { SANKEYPLOTS } from '../../modules/local/sankeyplots.nf'
include { OMICFLOW } from '../../modules/local/omicsflow.nf'


workflow REPORT {
    take:

    main:
    ch_versions = Channel.empty() 

    BIOTAVIZ()  
    ch_versions = ch_versions.mix(BIOTAVIZ.out.versions)

    SANKEYPLOTS()
    ch_versions = ch_versions.mix(SANKEYPLOTS.out.versions)

    OMICFLOW()
    ch_versions = ch_versions.mix(OMICFLOW.out.versions)

    emit:
    report                          = OMICFLOW.out.report
    versions                        = ch_versions
}
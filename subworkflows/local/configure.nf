/*

    Checks if classifier is downloaded
    
*/
include { CLASSIFIER_DOWNLOAD }     from    '../../modules/local/classifier/download.nf'
include { CLASSIFIER_CHECK }         from    '../../modules/local/classifier/check.nf'        
include { ensureDir }               from    '../../lib/utils.groovy'

workflow CONFIGURE {
    main:
    classifier_ch = Channel.empty()
    ch_versions = Channel.empty()

    if (!params.classifier_custom) {
        def classifier_path = ensureDir(params.classifier_path)

        CLASSIFIER_DOWNLOAD(
            classifier_path,
            params.classifier_name
        ).classifier.set{ classifier_ch }
    } else {
        classifier_ch = Channel
            .fromPath(params.classifier_custom)
            .ifEmpty { exit 1, 'Cannot find file: ${params.classifier_custom}\n'}
        
        CLASSIFIER_CHECK(
            classifier_ch,
            params.classifier_test_set
        )
        ch_versions = ch_versions.mix(CLASSIFIER_CHECK.out.versions)
    }

    emit:
    classifier          =   classifier_ch
    versions            =   ch_versions
}
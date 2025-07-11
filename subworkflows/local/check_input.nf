/*

    Fetches sample reads and puts them into the right structure
    
*/

include { samplesheetToList } from 'plugin/nf-schema'

workflow CHECK_INPUT {

    main:
    if (params.input) {
        def input = file("${params.input}", checkIfExists: true)
        def schema = file("${projectDir}/assets/schema_input.json", checkIfExists: true)
        def sample_ch = Channel.fromList(samplesheetToList(input, schema))
        metadata_ch = Channel.fromList(samplesheetToMetadata(input))

        meta_ch = sample_ch.map { arrayList ->
            def sample = arrayList[0]
            def files = params.singleEnd ? arrayList[1] : arrayList[1..2]
            def meta = [:]
            meta.id = sample.id
            meta.single_end = params.singleEnd
            return tuple(meta, files)
        }

    log.info "meta channel from samplesheet"

    } else if (params.reads) {
        sample_ch = Channel
            .fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2, checkIfExists: true)
            .ifEmpty { exit 1, 'Cannot find any reads matching: ${reads}\n'}
        meta_ch = sample_ch.map { arrayList ->
            def sample_id = arrayList[0]
            def files = arrayList[1]
            def meta = [:]
            meta.id = sample_id
            meta.single_end = params.singleEnd
            return tuple(meta, files)
        }

        metadata_ch = Channel.empty()
        log.info "meta channel from directory"
    }

    classifier_ch = Channel
        .fromPath(params.classifier)
        .ifEmpty { exit 1, 'Cannot find file: ${params.classifier}\n'}


    emit:
    meta        = meta_ch
    metadata    = metadata_ch
    classifier  = classifier_ch
}

def samplesheetToMetadata(input) {
    def rows = []
    input.withReader { reader ->
        def headers = reader.readLine().split(',').collect { it.trim() }
        reader.eachLine { line ->
            def values = line.split(',').collect { it.trim() }
            def row = [:]
            headers.eachWithIndex { h, i -> row[h] = values[i] }
            rows << row
        }
    }

    def isNumeric = { str ->
        str ==~ /^-?\d+(\.\d+)?$/
    }

    // Check types for each column
    def columnTypes = [:]
    if (rows) {
        def headers = rows[0].keySet()
        headers.each { col ->
            def values = rows.collect { it[col] }
            def allNumeric = values.every { v -> isNumeric(v) }
            def allString = values.every { v -> v instanceof String && !isNumeric(v) }
            columnTypes[col] = allNumeric ? 'numeric' : (allString ? 'string' : 'mixed')
        }
    }

    return [rows]
}
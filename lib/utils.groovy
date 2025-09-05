def save_output(input_ch, sub_dir_name) {
    input_ch.map { item ->
        def (meta, files) = (item instanceof List && item.size() == 2) ? [item[0], item[1]] : [null, item]
        def outDir = file("${params.outdir}/${sub_dir_name}")
        outDir.mkdir()

        // files is a list
        if (files instanceof List) {
            files.each { inputFile ->
               file(inputFile).copyTo(file("${outDir}/${file(inputFile).getName()}"))
            }
        } else {
            file(files).copyTo(file("${outDir}/${file(files).getName()}"))
        }
    }
}
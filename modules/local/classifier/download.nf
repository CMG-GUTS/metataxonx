process CLASSIFIER_DOWNLOAD {
    tag "$db_name"
    label 'process_single'

    input:
    path db_dir
    val db_name

    output:
    path("${db_dir}/${db_name}_classifier.qza"),   emit: classifier

    when:
    task.ext.when == null || task.ext.when

    script:
    def FILE = "${db_dir}/${db_name}_classifier.qza"
    """
    set -e

    if [ -f "${FILE}" ]; then
        echo "File ${FILE} already exists, skipping download."
        exit 0
    fi

    # Define URLs and expected SHA256 sums for known db_name values
    case "${db_name}" in
        "standard")
            URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
            SHA256_EXPECTED="c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616"
            ;;
        "diverse")
            URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-diverse-weighted-classifier.qza"
            SHA256_EXPECTED="decfae408061fab8ff2fec7dac1fe2a2e0041581589715062cc789bd4f9933db"
            ;;
        "human-stool")
            URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-human-stool-weighted-classifier.qza"
            SHA256_EXPECTED="db9e3c0105b1b9173deaa8bd828113b422c467443587cc8be3aed2e6f7cc995f"
            ;;
        "greengenes")
            URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza"
            SHA256_EXPECTED="5c5900c10ad601b12cf0d2280d52cc0f7852627e89070595cc3899a6d690ee27"
            ;;
        *)
            echo "Unknown db_name: ${db_name}"
            exit 1
            ;;
    esac

    # Download the classifier file into db_dir
    mkdir -p ${db_dir}
    wget \$URL -O ${FILE}

    # Calculate SHA256 and compare to expected
    SHA256_ACTUAL=\$(sha256sum ${FILE} | cut -d ' ' -f1)
    if [ "\$SHA256_ACTUAL" != "\$SHA256_EXPECTED" ]; then
        echo "SHA256 checksum mismatch for ${FILE}"
        echo "Expected: \$SHA256_EXPECTED"
        echo "Actual  : \$SHA256_ACTUAL"
        exit 1
    fi

    echo "Download and checksum verification succeeded for ${FILE}"
    """

    stub:
    def FILE = "${db_dir}/${db_name}_classifier.qza"
    """
    mkdir -p ${db_dir}
    touch ${FILE}
    """
}

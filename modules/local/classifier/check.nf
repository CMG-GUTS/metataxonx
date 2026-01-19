process CLASSIFIER_CHECK {
    label 'process_single'

    input:
    path classifier_path
    path test_data

    output:
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime tools import \\
        --type 'FeatureData[Sequence]' \\
        --input-path ${test_data} \\
        --output-path test_query.qza

    qiime feature-classifier classify-sklearn \\
        --i-classifier ${classifier_path} \\
        --i-reads test_query.qza \\
        --o-classification test_output.qza \\
        --p-n-jobs 1 2>&1 | tee qiime.log

    if grep -q "scikit-learn.*does not match" qiime.log; then
        echo "SKLEARN INCOMPATIBLE VERSION: ${classifier_path}, please use classifiers compatible to scikit-learn 1.4.2 version."
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime: \$(qiime --version | head -1 | sed -e "s/q2cli version //g")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
     "${task.process}":
          qiime: stub-version
     END_VERSIONS
    """
}

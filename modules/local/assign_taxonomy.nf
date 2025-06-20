process ASSIGN_TAXONOMY {
     tag "BATCH"
     label 'process_high'

     input:
     path(classifier)
     path(sequences)

     output:
     path "taxonomy.tsv"           , emit: taxonomy_tsv
     path "taxonomy_sklearn.qza"   , emit: taxonomy_qza
     path "versions.yml"           , emit: versions

     script:
     """
     qiime feature-classifier classify-sklearn \\
          --i-classifier $classifier \\
          --i-reads $sequences \\
          --o-classification taxonomy_sklearn.qza \\
          --p-n-jobs ${task.cpus}

     qiime metadata tabulate \\
          --m-input-file taxonomy_sklearn.qza \\
          --o-visualization taxonomy_sklearn.qzv

     qiime tools export --input-path taxonomy_sklearn.qza  \\
          --output-path taxonomy

     qiime tools export --input-path taxonomy_sklearn.qzv  \\
          --output-path taxonomy

     cat taxonomy/taxonomy.tsv | \\
     sed 's/Feature ID/#OTUID/' | \\
     sed 's/Taxon/taxonomy/' | \\
     sed 's/Consensus/consensus/' > taxonomy/taxonomy_relabeled.tsv
     mv taxonomy/taxonomy_relabeled.tsv ./taxonomy.tsv
     
     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
          qiime: \$(qiime --version)
     END_VERSIONS
     """

     stub:
     """
     touch taxonomy.tsv
     touch taxonomy_sklearn.qza

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
          qiime: \$(qiime --version)
     END_VERSIONS
     """
}
process SLEUTH_QC {
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::r-base=4.0 bioconda::bioconductor-deseq2=1.28.0 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://plateau-gao/rnaseq/rnaseq_sleuth:latest' :
        'littleplateau/rnaseq_sleuth' }"    
    
    input:
    path ("kallisto/*")
    path (tx2gene)
    path (pca_header_multiqc)
    path (clustering_header_multiqc)

    output:
    path "*.pdf"                , optional: true, emit: pdf
    path "*.RData"              , optional: true, emit: rdata
    path "*.pca.vals.txt"       , optional: true, emit: pca_txt
    path "*.pca.vals.mqc.tsv"   , optional: true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "norm_factors"         , optional:true, emit: norm_factors
    path "versions.yml"         , emit: versions
   
    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    def label_upper = args2.toUpperCase()
    """
    Rscript ${projectDir}/bin/sleuth_qc.r \\
        --h5_dir kallisto \\
        --tx2gene $tx2gene \\
        --outdir ./ \\
        --cores $task.cpus \\
        $args

    if [ -f "R_sessionInfo.log" ]; then
        sed "s/sleuth_pca/${label_lower}_sleuth_pca/g" <$pca_header_multiqc >tmp.txt
        sed -i -e "s/Sleuth PCA/${label_upper} Sleuth PCA/g" tmp.txt
        cat tmp.txt *.pca.vals.txt > ${label_lower}.pca.vals_mqc.tsv

        sed "s/sleuth_clustering/${label_lower}_sleuth_clustering/g" <$clustering_header_multiqc >tmp.txt
        sed -i -e "s/Sleuth sample/${label_upper} Sleuth sample/g" tmp.txt
        cat tmp.txt *.sample.dists.txt > ${label_lower}.sample.dists_mqc.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-sleuth: \$(Rscript -e "library(sleuth); cat(as.character(packageVersion('sleuth')))")
    END_VERSIONS    
    """
}
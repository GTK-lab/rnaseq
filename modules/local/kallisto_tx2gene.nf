process KALLISTO_TX2GENE {
    tag "$gtf"
    label "process_low"

    conda (params.enable_conda ? "bioconda::bioconductor-genomicfeatures=1.46.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-genomicfeatures:1.46.1--r41hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-genomicfeatures:1.46.1--r41hdfd78af_0' }"

    input:
    path gtf

    output:
    path "*.csv"       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    """
    kallisto_tx2gene.r ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-GenomicFeatures: \$(Rscript -e "library(GenomicFueatures); cat(as.character(packageVersion('GenomicFeatures')))")
    END_VERSIONS
    """

}
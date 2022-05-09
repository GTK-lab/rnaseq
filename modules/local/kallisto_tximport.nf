process KALLISTO_TXIMPORT {
    label "process_medium"

    conda (params.enable_conda ? "bioconda::bioconductor-tximeta=1.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0' :
        'quay.io/biocontainers/bioconductor-tximeta:1.8.0--r40_0' }"

    input:
    path ("salmon/*")
    path (tx2gene)

    output:
    path "*gene_tpm.tsv"                   , emit: tpm_gene
    path "*gene_counts.tsv"                , emit: counts_gene
    path "*gene_counts_length_scaled.tsv"  , emit: counts_gene_length_scaled
    path "*gene_counts_scaled.tsv"         , emit: counts_gene_scaled
    path "*transcript_tpm.tsv"             , emit: tpm_transcript
    path "*transcript_counts.tsv"          , emit: counts_transcript
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when  == null || task.ext.when

    script:
    """
    Rscript ../../bin/kallisto_tximport.r \\
        NULL \\
        kallisto \\
        kallisto.merged
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS    
    """
}
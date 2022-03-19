process KALLISTO_ABUNDANCE {
    tag "${task.process} Attempt_${task.attempt}_cpus_${task.cpus}_mem_${task.memory}"
    container 'quay.io/biocontainers/bioconductor-tximeta:1.8.0--r40_0' 
    
    publishDir "${params.outdir}/${task.process}",
            mode: 'copy', pattern: "*.{csv,txt,yml}"

    cpus { 12 * task.attempt }
    memory { 24.GB * task.attempt }

    input:
    path  tx2gene
    path ("*")
    
    output:
    path("*.{csv,txt}")
    path "versions.yml"                  , emit: versions
       
    script: 
    """
    mkdir kallisto
    mv *.kallisto_abundance.txt kallisto/
    FILES=`find -L ./ -name "*_abundance.txt"`
    echo \$FILES > check_files.txt

    kallisto_tximport.R \\
        $tx2gene    \\
        ./kallisto \\
        1>stdout.txt \\
        2>stderr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        container: "${task.container}"
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}

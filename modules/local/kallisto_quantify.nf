process KALLISTO_QUANTIFY {
    tag "$meta.id"
    label "process_high"

    conda(params.enable_conda? "bioconda::kallisto=0.48.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h0d531b0_1' :
        'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_0' }"
    
    input:
    tuple val(meta), path(reads)
    path (kallisto_index)

    output:
    tuple val(meta), path('*.kallisto_abundance.*')          , emit: kallisto_abundance
    tuple val(meta), path('*.kallisto_stderr.txt')           , emit: kallisto_stderr
    tuple val(meta), path('*.kallisto_stdout.txt')           , emit: kallisto_stdout
    tuple val(meta), path('*.kallisto_json_log.txt')         , emit: kallisto_json_log
    path "versions.yml"                                      , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        kallisto quant \\
            --single \\
            $agrs \\
            -i ${kallisto_index} \\
            -t ${task.cpus} \\
            -o kallisto_results \\
            ${prefix}.fastq.gz \\
            1> ${prefix}_kallisto_stdout.txt \\
            2> ${prefix}_kallisto_stderr.txt

        if $bootstrap>0
        then   
            mv kallisto_results/*.h5 ${prefix}.kallisto_abundance.h5
        else
            mv kallisto_results/*.tsv ${prefix}.kallisto_abundance.tsv
        fi
        cp kallisto_results/run_info.json ${prefix}.kallisto_json_log.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

        kallisto quant \\
            $args \\
            -i ${kallisto_index} \\
            -t ${task.cpus} \\
            -o kallisto_results \\
            ${reads[0]} \\
            ${reads[1]} \\
            1> ${prefix}.kallisto_stdout.txt \\
            2> ${prefix}.kallisto_stderr.txt

        if $bootstrap>0
        then   
            mv kallisto_results/*.h5 ${prefix}.kallisto_abundance.h5
        else
            mv kallisto_results/*.tsv ${prefix}.kallisto_abundance.tsv
        fi
        cp kallisto_results/run_info.json ${prefix}.kallisto_json_log.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
        END_VERSIONS
        """
    }
}
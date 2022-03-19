process KALLISTO_QUANT {
    tag "$sample_id Attempt_${task.attempt}_cpus_${task.cpus}_mem_${task.memory}"
    container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'

    publishDir "${params.outdir}/${task.process}",
        mode: 'copy', pattern: "*.{txt,yml}"

    cpus { 12 * task.attempt }
    memory { 24.GB * task.attempt }
    
    input:
    path(kallisto_index)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.kallisto_abundance.txt')        , emit: kallisto_abundance
    tuple val(sample_id), path('*.kallisto_stderr.txt')  , emit: kallisto_stderr_mqc
    tuple val(sample_id), path('*.kallisto_stdout.txt')  , emit: kallisto_stdout
    tuple val(sample_id), path('*.kallisto_json_log.txt')     , emit: kallisto_json_log
    path "versions.yml"                                  , emit: versions

    script:
    if (params.single_end) {
        """
        [ ! -f  ${sample_id}_1.fastq.gz ] && ln -s $reads ${sample_id}_1.fastq.gz

        kallisto quant \\
            -i ${kallisto_index} \\
            -o kallisto_results \\
            --seed=1 \\
            --plaintext \\
            -t ${task.cpus} \\
            ${sample_id}_1.fastq.gz \\
            1> ${sample_id}.kallisto_stdout.txt \\
            2> ${sample_id}.kallisto_stderr.txt

        mv kallisto_results/*.tsv ${sample_id}.kallisto_abundance.txt
        cp kallisto_results/run_info.json ${sample_id}.kallisto_json_log.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
            container: "${task.container}"
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${sample_id}_1.fastq.gz ] && ln -s ${reads[0]} ${sample_id}_1.fastq.gz
        [ ! -f  ${sample_id}_2.fastq.gz ] && ln -s ${reads[1]} ${sample_id}_2.fastq.gz

        kallisto quant \\
            -i ${kallisto_index} \\
            -o kallisto_results \\
            --seed=1 \\
            --plaintext \\
            -t ${task.cpus} \\
            ${sample_id}_1.fastq.gz \\
            ${sample_id}_2.fastq.gz \\
            1> ${sample_id}.kallisto_stdout.txt \\
            2> ${sample_id}.kallisto_stderr.txt

        mv kallisto_results/*.tsv ${sample_id}.kallisto_abundance.txt
        cp kallisto_results/run_info.json ${sample_id}.kallisto_json_log.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
            container: "${task.container}"
        END_VERSIONS
        """
    }
}

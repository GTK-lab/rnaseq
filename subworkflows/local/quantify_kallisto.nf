/*
===============================================================================
    KALLISTO QUANTIFICATION
===============================================================================
*/
include { QUANTIFY        } from '../modules/kallisto_quantify'

workflow KALLISTO {
    take:
    reads             // [val(meta), path[reads]]
    kallisto_index    // path(index)
    
    main:
    ch_versions    = Channel.empty()
    index          = Channel.empty()


    kallisto_h5    = Channel.empty()
    kallisto_stderr = Channel.empty()
    kallisto_stdout = Channel.empty()
    kallisto_log    = Channel.empty()

    QUANTIFY(reads, index)
    kallisto_h5    = QUANTIFY.out.kallisto_abundance_h5
    kallisto_stderr = QUANTIFY.out.kallisto_stderr_mqc
    kallisto_stdout = QUANTIFY.out.kallisto_stdout
    kallisto_log    = QUANTIFY.out.kallisto_json_log
    ch_versions     = ch_versions.mix(QUANTIFY.out.versions)

    emit:
    index

    kallisto_h5
    kallisto_stderr
    kallisto_stdout
    kallisto_log

    versions = ch_versions     // channel: [ versions.yml ]
}


/*
===============================================================================
    KALLISTO QUANTIFICATION
===============================================================================
*/
include { KALLISTO_QUANTIFY       } from '../../modules/local/kallisto_quantify'
include { KALLISTO_TX2GENE        } from '../../modules/local/kallisto_tx2gene'
include { KALLISTO_TXIMPORT       } from '../../modules/local/kallisto_tximport'

include { KALLISTO_SUMMARIZEDEXPERIMENT as KALLISTO_SE_GENE               } from '../../modules/local/kallisto_summarizedexperiment'
include { KALLISTO_SUMMARIZEDEXPERIMENT as KALLISTO_SE_GENE_LENGTH_SCALED } from '../../modules/local/kallisto_summarizedexperiment'
include { KALLISTO_SUMMARIZEDEXPERIMENT as KALLISTO_SE_GENE_SCALED        } from '../../modules/local/kallisto_summarizedexperiment'
include { KALLISTO_SUMMARIZEDEXPERIMENT as KALLISTO_SE_GENE_TRANSCRIPT    } from '../../modules/local/kallisto_summarizedexperiment'


workflow QUANTIFY_KALLISTO {
    take:
    reads             // [val(meta), path[reads]]
    index             // kallisto_index
    gtf               // genome.gtf

    main:
    ch_versions    = Channel.empty()

    KALLISTO_QUANTIFY (reads, index)
    ch_versions = ch_versions.mix(KALLISTO_QUANTIFY.out.versions)

    emit:
    index

    kallisto_h5
    kallisto_stderr
    kallisto_stdout
    kallisto_log

    versions = ch_versions     // channel: [ versions.yml ]
}


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

    KALLISTO_TX2GENE ( gtf )
    ch_versions = ch_versions.mix(KALLISTO_TX2GENE.out.versions)

    KALLISTO_TXIMPORT ( KALLISTO_QUANTIFY.out.results.collect{it[1]}, KALLISTO_TX2GENE.out.tsv.collect() )
    ch_versions = ch_versions.mix(KALLISTO_TXIMPORT.out.versions)

    KALLISTO_SE_GENE (
        SALMON_TXIMPORT.out.counts_gene,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )
    ch_versions = ch_versions.mix(KALLISTO_SE_GENE.out.versions)

    KALLISTO_SE_GENE_LENGTH_SCALED (
        SALMON_TXIMPORT.out.counts_gene_length_scaled,
        SALMON_TXIMPORT.out.tpm_gene_length_scaled,
        SALMON_TX2GENE.out.tsv.collect()
    )

    KALLISTO_SE_GENE_SCALED (
        SALMON_TXIMPORT.out.counts_gene_scaled,
        SALMON_TXIMPORT.out.tpm_gene_scaled,
        SALMON_TX2GENE.out.tsv.collect()
    )

    KALLISTO_SE_TRANSCRIPT (
        SALMON_TXIMPORT.out.counts_transcript,
        SALMON_TXIMPORT.out.tpm_transcript,
        SALMON_TX2GENE.out.tsv.collect()

    emit:
    results                       = KALLISTO_QUANTIFY.out.results                   // channel: [ val(meta), results_dir ]

    tpm_gene                      = KALLISTO_TXIMPORT.out.tpm_gene                  // channel: [ val(meta), counts ]
    counts_gene                   = KALLISTO_TXIMPORT.out.counts_gene               // channel: [ val(meta), counts ]
    counts_gene_length_scaled     = KALLISTO_TXIMPORT.out.counts_gene_length_scaled // channel: [ val(meta), counts ]
    counts_gene_scaled            = KALLISTO_TXIMPORT.out.counts_gene_scaled        // channel: [ val(meta), counts ]
    tpm_transcript                = KALLISTO_TXIMPORT.out.tpm_transcript            // channel: [ val(meta), counts ]
    counts_transcript             = KALLISTO_TXIMPORT.out.counts_transcript         // channel: [ val(meta), counts ]

    merged_gene_rds               = KALLISTO_SE_GENE.out.rds                        //    path: *.rds
    merged_gene_rds_length_scaled = KALLISTO_SE_GENE_LENGTH_SCALED.out.rds          //    path: *.rds
    merged_gene_rds_scaled        = KALLISTO_SE_GENE_SCALED.out.rds                 //    path: *.rds

    merged_counts_transcript      = KALLISTO_TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript         = KALLISTO_TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv
    merged_transcript_rds         = KALLISTO_SE_TRANSCRIPT.out.rds                  //    path: *.rds

    versions = ch_versions     // channel: [ versions.yml ]
}


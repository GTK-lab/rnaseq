nextflow.enable.dsl=2


ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)


include { PREPARE_GENOME     } from './subworkflows/local/prepare_genome'
include { QUANTIFY_KALLISTO  } from './subworkflows/local/quantify_kallisto'

include { FASTQC_UMITOOLS_TRIMGALORE } from './subworkflows/nf-core/fastqc_umitools_trimgalore'
include { MARK_DUPLICATES_PICARD     } from './subworkflows/nf-core/mark_duplicates_picard'

include { INPUT_CHECK    } from './subworkflows/local/input_check'
include { CAT_FASTQ                   } from './modules/nf-core/modules/cat/fastq/main'

include { DESEQ2_QC as DESEQ2_QC_KALLISTO }  from './modules/local/deseq2_qc'
include { SLEUTH_QC as SLEUTH_QC_KALLISTO }  from './modules/local/sleuth_qc'

def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

def prepareToolIndices  = []
if (!params.skip_bbsplit)   { prepareToolIndices << 'bbsplit'             }
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }
if (params.pseudo_aligner)  { prepareToolIndices << params.pseudo_aligner }

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

workflow TEST_KALLISTO{
    ch_versions = Channel.empty()

    
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME(
        prepareToolIndices,
        biotype
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))


    FASTQC_UMITOOLS_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_trimming,
        params.umi_discard_read
    )
   ch_versions = ch_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.versions)


    ch_filtered_reads = FASTQC_UMITOOLS_TRIMGALORE.out.reads 

    ch_kallisto_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc          = Channel.empty()
    ch_pseudoaligner_clustering_multiqc   = Channel.empty()
    if(params.pseudo_aligner == 'kallisto') {
        QUANTIFY_KALLISTO(
            ch_filtered_reads,
            PREPARE_GENOME.out.kallisto_index,
            PREPARE_GENOME.out.gtf
        )
        ch_kallisto_multiqc = QUANTIFY_KALLISTO.out.results
        ch_versions         = ch_versions.mix(QUANTIFY_KALLISTO.out.versions)
    
        if ( !params.skip_qc ) {
            if (params.kallisto_bootstrap > 0) {
                SLEUTH_QC_KALLISTO (
                    QUANTIFY_KALLISTO.out.results,
                    QUANTIFY_KALLISTO.out.tx2gene
                    // ch_pca_header_multiqc, 
                    // ch_clustering_header_multiqc               
                )
                // ch_pseudoaligner_pca_multiqc        = SLEUTH_QC_KALLISTO.out.pca_multiqc
                // ch_pseudoaligner_clustering_multiqc = SLEUTH_QC_KALLISTO.out.dists_multiqc
                ch_versions = ch_versions.mix(SLEUTH_QC_KALLISTO.out.versions)
            } else if (!params.skip_deseq2_qc) {
                DESEQ2_QC_KALLISTO (
                    QUANTIFY_KALLISTO.out.counts_gene_length_scaled,
                    ch_pca_header_multiqc,
                    ch_clustering_header_multiqc
                )
                ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_KALLISTO.out.pca_multiqc
                ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_KALLISTO.out.dists_multiqc
                ch_versions = ch_versions.mix(DESEQ2_QC_KALLISTO.out.versions)
            }
        }
    }
}

workflow {
    TEST_KALLISTO()
}

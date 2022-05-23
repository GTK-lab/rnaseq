nextflow.enable.dsl=2

include { PREPARE_GENOME     } from './subworkflows/local/prepare_genome'
include { QUANTIFY_KALLISTO  } from './subworkflows/local/quantify_kallisto'

include { FASTQC_UMITOOLS_TRIMGALORE } from './subworkflows/nf-core/fastqc_umitools_trimgalore'
include { MARK_DUPLICATES_PICARD     } from './subworkflows/nf-core/mark_duplicates_picard'

include { INPUT_CHECK    } from './subworkflows/local/input_check'
include { CAT_FASTQ                   } from './modules/nf-core/modules/cat/fastq/main'

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
   // ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
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

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    FASTQC_UMITOOLS_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_trimming,
        params.umi_discard_read
    )
   // ch_versions = ch_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.versions)


    ch_filtered_reads = FASTQC_UMITOOLS_TRIMGALORE.out.reads 

    if(params.pseudo_aligner == 'kallisto') {
        QUANTIFY_KALLISTO(
            ch_filtered_reads,
            PREPARE_GENOME.out.kallisto_index,
            PREPARE_GENOME.out.gtf
        )
    }   
}

workflow {
    TEST_KALLISTO()
}
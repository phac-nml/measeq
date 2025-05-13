/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_measeq_pipeline'
include { SETUP_REFERENCE_DATA    } from '../subworkflows/local/setup_reference_data'
include { GENERATE_CONSENSUS      } from '../subworkflows/local/generate_consensus'
include { ARTIC_GET_MODELS        } from '../modules/local/artic/get_model/main'
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { NANOQ                   } from '../modules/nf-core/nanoq/main'
include { MINIMAP2_ALIGN          } from '../modules/local/minimap2/main'
include { NEXTCLADE_DATASETGET    } from '../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_N450   } from '../modules/nf-core/nextclade/run/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_CUSTOM } from '../modules/nf-core/nextclade/run/main'
include { ADJUST_FASTA_HEADER     } from '../modules/local/artic/subcommands/main'
include { SAMTOOLS_DEPTH          } from '../modules/nf-core/samtools/depth/main'
include { MAKE_SAMPLE_QC_CSV      } from '../modules/local/qc/sample/main'
include { MAKE_FINAL_QC_CSV       } from '../modules/local/qc/summary/main'
include { GENERATE_REPORT         } from '../subworkflows/local/generate_report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MEASEQ {

    take:
    ch_input_fastqs // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_metadata = params.metadata ? file(params.metadata, type: 'file', checkIfExists: true) : []
    Channel
        .value(file(params.custom_nextclade_dataset, type: 'dir', checkIfExists: true))
        .set { ch_custom_nextclade_dataset }

    //
    // MODULE: Setup nextclade dataset
    //
    NEXTCLADE_DATASETGET(
        params.nextclade_dataset_name,
        params.nextclade_dataset_tag
    )
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

    //
    // WORKFLOW: Reference Setup
    //
    SETUP_REFERENCE_DATA(
        params.reference,
        params.primer_bed,
        NEXTCLADE_DATASETGET.out.dataset
    )
    ch_reference                = SETUP_REFERENCE_DATA.out.reference
    ch_reference_fai            = SETUP_REFERENCE_DATA.out.fai
    ch_primer_bed               = SETUP_REFERENCE_DATA.out.primer_bed
    ch_split_amp_pools_bed      = SETUP_REFERENCE_DATA.out.split_amp_pools_bed
    ch_ref_n450                 = SETUP_REFERENCE_DATA.out.ref_n450
    ch_strain                   = SETUP_REFERENCE_DATA.out.strain
    ch_versions                 = ch_versions.mix(SETUP_REFERENCE_DATA.out.versions)

    //
    // MODULE: Model download if not a local one
    //
    if ( params.local_model ) {
        ch_model = Channel.value(file(params.local_model, type: 'dir', checkIfExists: true))
    } else {
        ARTIC_GET_MODELS(
            params.model
        )
        ch_model = ARTIC_GET_MODELS.out.model
        ch_versions = ch_versions.mix(ARTIC_GET_MODELS.out.versions)
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC(
        ch_input_fastqs
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Nanoq
    //
    NANOQ(
        ch_input_fastqs,
        'fastq'
    )
    ch_versions = ch_versions.mix(NANOQ.out.versions.first())

    //
    // MODULE: Run Minimap2
    //
    MINIMAP2_ALIGN(
        NANOQ.out.reads,
        ch_reference
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //
    // WORKFLOW: Generate Consensus
    //
    GENERATE_CONSENSUS(
        ch_reference,
        ch_reference_fai,
        MINIMAP2_ALIGN.out.bam_bai,
        ch_model,
        ch_primer_bed,
        ch_split_amp_pools_bed
    )
    ch_consensus = GENERATE_CONSENSUS.out.consensus
    ch_bam_bai   = GENERATE_CONSENSUS.out.bam_bai
    ch_vcf       = GENERATE_CONSENSUS.out.vcf
    ch_versions  = ch_versions.mix(GENERATE_CONSENSUS.out.versions)

    //
    // MODULE: Nextclade Run
    //
    NEXTCLADE_RUN_N450(
        ch_consensus,
        NEXTCLADE_DATASETGET.out.dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN_N450.out.versions.first())

    // Using the N450 nextclade align, create renamed, easy to find N450 output
    ADJUST_FASTA_HEADER(
        NEXTCLADE_RUN_N450.out.fasta_aligned,
        ch_reference,
        '.N450',
        '-N450'
    )
    ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions.first())

    //
    // QC
    //
    SAMTOOLS_DEPTH(
        ch_bam_bai.map{ it -> [it[0], it[1]] },
        [[:], []] // Empty as want to run whole depth
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    NEXTCLADE_RUN_CUSTOM(
        ch_consensus,
        ch_custom_nextclade_dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN_CUSTOM.out.versions.first())

    MAKE_SAMPLE_QC_CSV(
        ch_bam_bai
            .join(ch_consensus, by: [0])
            .join(SAMTOOLS_DEPTH.out.tsv, by: [0])
            .join(NANOQ.out.stats, by: [0])
            .join(NEXTCLADE_RUN_N450.out.csv, by: [0])
            .join(NEXTCLADE_RUN_CUSTOM.out.csv, by: [0])
            .join(ch_vcf, by: [0]),
        ch_strain,
        ch_primer_bed.ifEmpty([])
    )
    ch_versions = ch_versions.mix(MAKE_SAMPLE_QC_CSV.out.versions.first())

    MAKE_FINAL_QC_CSV(
        MAKE_SAMPLE_QC_CSV.out.csv
            .map{ it -> it[1] }
            .collectFile(keepHeader: true, skip: 1, name: 'concat.qc.csv'),
        ch_metadata,
        params.neg_control_threshold,
        params.neg_ctrl_substrings
    )
    ch_versions = ch_versions.mix(MAKE_FINAL_QC_CSV.out.versions)

    //
    // Report
    //
    GENERATE_REPORT(
        ch_reference,
        ch_reference_fai,
        ch_bam_bai,
        ch_strain,
        SAMTOOLS_DEPTH.out.tsv,
        MAKE_FINAL_QC_CSV.out.csv
    )
    ch_versions = ch_versions.mix(GENERATE_REPORT.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'measeq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

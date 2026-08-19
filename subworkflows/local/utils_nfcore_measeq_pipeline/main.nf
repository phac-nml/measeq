//
// Subworkflow with functionality specific to the phac-nml/measeq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:
    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Summarize and Validate Params
    //
    if (validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    //
    // Create channel from input file provided through params.input
    //
    def processedIDs = [] as Set
    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                // Meta ID assignment
                if (!meta.id) {
                    meta.id = meta.irida_id
                } else {
                    meta.id = meta.id.replaceAll(/[^A-Za-z0-9_.\-]/, '_')
                }

                // Ensure ID is unique by appending meta.irida_id if needed
                meta.id = makeUniqueId(meta.id, meta.irida_id, processedIDs)

                // Add the ID to the set of processed IDs
                processedIDs << meta.id

                // File assignment
                if (!fastq_2) {
                    return tuple(meta.id, meta + [ single_end:true ], [ fastq_1 ])
                } else {
                    return tuple(meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ])
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return tuple(meta, fastqs.flatten())
        }

    emit:
    samplesheet = ch_samplesheet    // channel: samplesheet read in from --input // [ meta(id, single-end), fastqs(f1,f2) ]
    versions    = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir          // path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Make the meta.id unique by appending meta.irida_id, if needed
//
def makeUniqueId(id, irida_id, processedIDs) {
    if (processedIDs.contains(id)) {
        return makeUniqueId("${id}_${irida_id}", irida_id, processedIDs)
    }
    return id
}

//
// Subworkflow with functionality specific to the phac-nml/measeq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { PREDICT_GENOTYPE          } from '../../../modules/local/predict_genotype/main.nf'
include { STAGE_FILE_IRIDANEXT      } from '../../../modules/local/custom/stage_file_iridanext/main.nf'

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
    Channel
        .value(file("$projectDir/assets/reference/n450/measles_N450_genotypes.mmi", type: 'file', checkIfExists: true))
        .set { ch_measles_n450_mmi }
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
    // Help
    //
    if (params.help) {
        log.info paramsHelp("nextflow run phac-nml/measeq -profile <profile> --input samplesheet.csv --platform <illumina|nanopore> --outdir <outdir>")
        exit 0
    }

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
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
                // Meta ID assignment
                if (!meta.id) {
                    meta.id = meta.irida_id
                } else {
                    meta.id = meta.id.replaceAll(/[^A-Za-z0-9_.\-]/, '_')
                }

                // Ensure ID is unique by appending meta.irida_id if needed
                //  Note that nextflow does not like while loops
                while (processedIDs.contains(meta.id)) {
                    meta.id = "${meta.id}_${meta.irida_id}"
                }
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
        .set { ch_samplesheet }

    // Prepare Function for reference detection
    // Function: Get FASTA header and use that as ref_id
    def fastaHeaderId = { Path fasta ->
        def header = fasta.withReader { r ->
            String line
            while ((line = r.readLine()) != null) {
                line = line.trim()
                if (line) return line
            }
            return null
        }
        assert header?.startsWith('>') : "Reference FASTA (${fasta}) must start with a header line"
        return header.substring(1).tokenize()[0]
    }

    //
    // Modify samplesheet according to params.reference / params.primer_bed
    // or predict genotype when no reference is provided.
    //
    if ( params.reference ) {
        // When params.reference is NOT null
        ch_reference = Channel
            .value( file(params.reference, type: 'file', checkIfExists: true) )
            .map { fasta ->
                ref_id = fastaHeaderId(fasta)
                tuple([ id: ref_id, irida_id: ref_id ], fasta)
            }

        // Create a channel with primers if the argument is passed or ouput an empty channel
        if ( params.amplicon && params.primer_bed || params.primer_bed ) {
            ch_primer_bed = ch_reference
                .map { meta_ref, fasta ->
                    tuple(meta_ref, file(params.primer_bed, type: 'file', checkIfExists: true))
                }
        } else if( params.amplicon && !params.primer_bed ) {
            error "Please provide a file with primers using '--primer_bed' when running with '--amplicon'"
        } else {
            ch_primer_bed = Channel.empty()
        }

        // Attach ref_id to every sample
        ch_samples = ch_samplesheet
            .combine(ch_reference.map { meta_ref, fasta -> meta_ref.id } )
            .map { meta, fastqs, meta_ref ->
                tuple(meta + [ ref_id: meta_ref ], fastqs)
            }

    } else {
        // When params.reference is null
        //
        // MODULE: predict genotype and assign reference
        //
        PREDICT_GENOTYPE(
            ch_samplesheet,
            ch_measles_n450_mmi
        )

        // Create predictions.csv file and save in a way so IRIDA Next can get it
        STAGE_FILE_IRIDANEXT(
            PREDICT_GENOTYPE.out.csv
                .collectFile(name: 'predictions.csv', keepHeader: true, skip: 1, cache: false, sort: true)
        )

        // Modify module input to create samples, reference, and primers channels
        ch_pred = PREDICT_GENOTYPE.out.samplesheet
            .map { meta, fastqs, genotype ->
                def ref_path = params."${genotype}_ref" ?: params.default_ref
                def fasta = file(ref_path, type: 'file', checkIfExists: true)
                def ref_id = fastaHeaderId(fasta)
                def primer_bed
                if ( !params.amplicon ) {
                    primer_bed = ""
                } else if ( params."${genotype}_ref" && params.amplicon ){
                    if ( params."${genotype}_bed" ) {
                        primer_bed = file(params."${genotype}_bed", type: 'file', checkIfExists: true)
                    } else {
                        error "You have specified --amplicon but the ${genotype} predicted genotype doesn't have a valid primer bed file. Provide --${genotype}_bed or remove --amplicon."
                    }
                } else if ( ref_path == params.default_ref && params.amplicon ) {
                    primer_bed = file(params.default_bed, type: 'file', checkIfExists: true)
                }
                tuple(meta + [ ref_id: ref_id ], fastqs, fasta, primer_bed)
            }

        // Create samples channel
        ch_samples = ch_pred.map { meta, fastqs, _fasta, _primer_bed ->
            tuple(meta, fastqs)
        }

        // Create reference channel
        ch_reference = ch_pred
            .map { _meta, _fastqs, fasta, _primer_bed ->
                tuple([ id: fastaHeaderId(fasta), irida_id: fastaHeaderId(fasta) ], fasta)
            }
            .unique { it[0].id }

        // Create primer_bed channel if needed
        if ( params.amplicon ) {
            ch_primer_bed = ch_pred
                .map { _meta, _fastqs, fasta, primer_bed ->
                    tuple([ id: fastaHeaderId(fasta), irida_id: fastaHeaderId(fasta) ], primer_bed)
                }
                .unique { it[0].id }
        } else {
            ch_primer_bed = Channel.empty()
        }
    }

    emit:
    samplesheet = ch_samples        // channel: samplesheet read in from --input // tuple meta[id, single-end, ref_id], fastqs[f1,f2]
    reference   = ch_reference      // channel: reference created from --reference or predicted // tuple meta[ref_id], fasta
    primer_bed  = ch_primer_bed     // channel: primer bed file created from --primer_bed or predicted // tuple meta[ref_id], primer_bed : null
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
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

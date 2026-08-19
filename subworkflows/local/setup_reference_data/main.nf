/*
Subworkflow with functionality to setup reference data
    - Create required files using reference sequence
    - Create amplicon files if using amplicon data
    - Genotype and get the N450 of the reference
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GENERATE_REF_INTERMEDIATES } from '../../../modules/local/input_utils/main'
include { SPLIT_AMPLICON_REGION      } from '../../../modules/local/input_utils/main'
include { GENERATE_AMPLICON_BED      } from '../../../modules/local/input_utils/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_REFERENCE } from '../../../modules/nf-core/nextclade/run/main'
include { ADJUST_FASTA_HEADER        } from '../../../modules/local/adjust_fasta_header/main'
include { EXTRACT_GENOTYPE           } from '../../../modules/local/custom/extract_genotype/main'
include { PREDICT_GENOTYPE           } from '../../../modules/local/predict_genotype/main.nf'
include { STAGE_FILE_IRIDANEXT       } from '../../../modules/local/custom/stage_file_iridanext/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SETUP_REFERENCE_DATA {

    take:
    ch_samplesheet      // channel: samplesheet read in from --input // [ meta(id, single-end), fastqs(f1,f2) ]
    nextclade_dataset   // channel: [ dataset ]

    main:
    ch_versions            = Channel.empty()
    ch_amplicon_bed        = Channel.empty()
    ch_split_amp_pools_bed = Channel.empty()
    ch_measles_n450_mmi    = Channel
        .value(file("$projectDir/assets/reference/n450/measles_N450_genotypes.mmi", type: 'file', checkIfExists: true))

    //
    // Modify samplesheet according to params.reference / params.primer_bed
    // or predict genotype when no reference is provided.
    //
    if ( params.reference ) {
        // When params.reference is NOT null
        ch_reference = Channel
            .value( file(params.reference, type: 'file', checkIfExists: true) )
            .map { fasta ->
                def ref_id = fastaHeaderId(fasta)
                tuple([ id: ref_id, irida_id: ref_id ], fasta)
            }

        // Create a channel with primers if the argument is passed or ouput an empty channel
        if ( params.amplicon && params.primer_bed || params.primer_bed ) {
            ch_primer_bed = ch_reference
                .map { meta_ref, _fasta ->
                    tuple(meta_ref, file(params.primer_bed, type: 'file', checkIfExists: true))
                }
        } else if( params.amplicon && !params.primer_bed ) {
            error "Please provide a file with primers using '--primer_bed' when running with '--amplicon'"
        } else {
            ch_primer_bed = Channel.empty()
        }

        // Attach ref_id to every sample
        ch_samples = ch_samplesheet
            .combine(ch_reference.map { meta_ref, _fasta -> meta_ref.id } )
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
                def ref_path = params["${genotype}_ref"] ?: params.default_ref
                def fasta = file(ref_path, type: 'file', checkIfExists: true)
                def ref_id = fastaHeaderId(fasta)
                def primer_bed
                if ( !params.amplicon ) {
                    primer_bed = ""
                } else if ( params["${genotype}_ref"] && params.amplicon ){
                    if ( params["${genotype}_bed"] ) {
                        primer_bed = file(params["${genotype}_bed"], type: 'file', checkIfExists: true)
                    } else {
                        error "You have specified --amplicon but the ${genotype} predicted genotype doesn't have a valid primer bed file. Provide --${genotype}_bed or remove --amplicon."
                    }
                } else if ( ref_path == params.default_ref && params.amplicon ) {
                    primer_bed = file(params.default_bed, type: 'file', checkIfExists: true)
                }
                tuple(meta + [ ref_id: ref_id ], fastqs, fasta, primer_bed)
            }

        // Create samples channel
        ch_samples = ch_pred
            .map { meta, fastqs, _fasta, _primer_bed ->
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

    //
    // MODULE: Generate intermediate reference files
    //
    GENERATE_REF_INTERMEDIATES(
        ch_reference
    )
    ch_versions = ch_versions.mix(GENERATE_REF_INTERMEDIATES.out.versions)

    //
    // Generate needed files if we have a primer scheme
    //
    if ( params.amplicon || params.primer_bed ) {
        //
        // MODULE: Generate amplicon bed for reporting
        //
        GENERATE_AMPLICON_BED(
            ch_primer_bed
        )
        ch_versions = ch_versions.mix(GENERATE_AMPLICON_BED.out.versions)
        ch_amplicon_bed = GENERATE_AMPLICON_BED.out.bed

        //
        // MODULE: Generate Bed files for clair3 and reporting
        //
        SPLIT_AMPLICON_REGION(
            ch_amplicon_bed
        )
        // Bed files are POOL.REF_ID.bed, we need to get the POOL
        //  Thus using a replace for everything after the first .
        ch_split_amp_pools_bed = SPLIT_AMPLICON_REGION.out.bed
            .flatMap { meta, bed_files ->
                bed_files.collect { bedF ->
                    def pool_name = bedF.name.split((/\./))[0]
                    tuple(meta + [ pool: pool_name ], file(bedF))
                }
            }

    }
    //
    // MODULE: Run nextclade on the reference to get the genotype and N450
    //
    NEXTCLADE_RUN_REFERENCE(
        ch_reference,
        nextclade_dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN_REFERENCE.out.versions)

    //
    // MODULE: Adjust fasta header to make the final output easier to work with downstream
    //
    ADJUST_FASTA_HEADER(
        NEXTCLADE_RUN_REFERENCE.out.fasta_aligned,
        [['id': ''], []],
        '.N450',
        '-N450'
    )

    //
    // MODULE: Extract genotype
    //
    EXTRACT_GENOTYPE(
        NEXTCLADE_RUN_REFERENCE.out.csv
    )
    ch_genotype = EXTRACT_GENOTYPE.out.genotype
        .map { meta, genotype -> tuple(meta.id, genotype.trim()) }

    // Add genotype to the reference channel
    ch_reference = ch_reference
        .map { meta, fasta -> tuple(meta.id, meta, fasta) }
        .join(ch_genotype, by: [0])
        .map { _meta_id, meta, fasta, genotype -> tuple(meta + [genotype: genotype], fasta) }

    emit:
    samples             = ch_samples                                 // channel: [ meta(id, single-end, ref_id), fastqs(f1,f2) ]
    primer_bed          = ch_primer_bed                              // channel: [ meta(ref_id), primer_bed ] ?: null
    reference           = ch_reference                               // channel: [ meta(ref_id, genotype), fasta ]
    fai                 = GENERATE_REF_INTERMEDIATES.out.fai         // channel: [ meta(ref_id), *.fai ]
    refstats            = GENERATE_REF_INTERMEDIATES.out.refstats    // channel: [ meta(ref_id), refstats ]
    genome_bed          = GENERATE_REF_INTERMEDIATES.out.genome_bed  // channel: [ meta(ref_id), genome.bed ]
    amplicon_bed        = ch_amplicon_bed                            // channel: [ meta(ref_id), amplicon.bed ]
    split_amp_pools_bed = ch_split_amp_pools_bed                     // channel: [ meta(ref_id, pool), *.bed ]
    ref_n450            = ADJUST_FASTA_HEADER.out.consensus          // channel: [ meta(ref_id), consensus ]
    versions            = ch_versions                                // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get FASTA header and use that as ref_id
//
def fastaHeaderId(fasta) {
    def header = fasta.withReader { r ->
        r.find { it.trim() }?.trim()
    }
    if (!header?.startsWith('>')) {
        error "Reference FASTA (${fasta}) must start with a header line"
    }
    return header.drop(1).tokenize()[0]
}

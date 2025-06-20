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
include { ADJUST_FASTA_HEADER        } from '../../../modules/local/artic/subcommands/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SETUP_REFERENCE_DATA {

    take:
    reference           // value: path to reference file or null
    bed                 // value: path to primer bed file or null
    nextclade_dataset   // channel: [ dataset ]

    main:
    ch_versions            = Channel.empty()
    ch_primer_bed          = bed ? Channel.value(file(bed, type: 'file', checkIfExists: true)) : Channel.empty()
    ch_split_amp_pools_bed = Channel.empty()
    ch_amplicon_bed        = Channel.empty()

    //
    // Logic - Need to have a reference file input at minimum
    //
    if ( !reference ) {
        error "Please provide a measles '--reference' fasta file to run"
    }

    //
    // Create reference meta-map for downstream steps
    //
    Channel
        .value(file(reference, type: 'file', checkIfExists: true))
        .map { file ->
            def firstLine = file.withReader { it.readLine() }
            def meta = [id: firstLine.split()[0].replace('>', '')]
            [meta, file]
        }
        .set { ch_reference }

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
    if ( ch_primer_bed ) {
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

        // Format amplicon bed files channel to be [ val(name), path(bed_file) ]
        //  Mostly for clair3 to run on specific pools
        SPLIT_AMPLICON_REGION.out.bed
                .flatten()
                .map{ bedF -> [ bedF.baseName, file(bedF) ] }
                .set { ch_split_amp_pools_bed }
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

    // Get the genotype from the nextclade output as a value channel
    NEXTCLADE_RUN_REFERENCE.out.csv
        .map{ meta, csv -> csv }
        .splitCsv(sep: ';', header: true, limit:1)
        .map{ row -> row.clade }
        .first()
        .set { ch_genotype }

    emit:
    reference           = ch_reference
    fai                 = GENERATE_REF_INTERMEDIATES.out.fai
    refstats            = GENERATE_REF_INTERMEDIATES.out.refstats
    genome_bed          = GENERATE_REF_INTERMEDIATES.out.genome_bed
    primer_bed          = ch_primer_bed
    amplicon_bed        = ch_amplicon_bed
    split_amp_pools_bed = ch_split_amp_pools_bed
    ref_n450            = ADJUST_FASTA_HEADER.out.consensus
    genotype              = ch_genotype
    versions            = ch_versions
}

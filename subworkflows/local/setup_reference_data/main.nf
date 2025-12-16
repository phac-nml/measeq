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
include { EXTRACT_GENOTYPE           } from '../../../modules/local/custom/extract_genotype/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SETUP_REFERENCE_DATA {

    take:
    ch_reference        // channel: [ meta[ref_id], fasta ]
    ch_primer_bed       // channel: [ meta[ref_id], primer_bed ?: null ]
    nextclade_dataset   // channel: [ dataset ]

    main:
    ch_versions            = Channel.empty()
    ch_amplicon_bed        = Channel.empty()
    ch_split_amp_pools_bed = Channel.empty()

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
    reference           = ch_reference                               // channel: [ meta(ref_id, genotype), fasta ]
    fai                 = GENERATE_REF_INTERMEDIATES.out.fai         // channel: [ meta(ref_id), *.fai ]
    refstats            = GENERATE_REF_INTERMEDIATES.out.refstats    // channel: [ meta(ref_id), refstats ]
    genome_bed          = GENERATE_REF_INTERMEDIATES.out.genome_bed  // channel: [ meta(ref_id), genome.bed ]
    amplicon_bed        = ch_amplicon_bed                            // channel: [ meta(ref_id), amplicon.bed ]
    split_amp_pools_bed = ch_split_amp_pools_bed                     // channel: [ meta(ref_id, pool), *.bed ]
    ref_n450            = ADJUST_FASTA_HEADER.out.consensus          // channel: [ meta(ref_id), consensus ]
    versions            = ch_versions                                // channel: [ path(versions.yml) ]
}

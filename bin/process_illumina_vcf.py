#!/usr/bin/env python
'''
Originally written by @jts from https://github.com/jts/ncov2019-artic-nf/blob/be26baedcc6876a798a599071bb25e0973261861/bin/process_gvcf.py

Adjustments made such to focus just on the variants and not the GVCF output due to a rare error in the gVCF that was seen when running measles genomes
    Full adjustments:
        - Added in genotype adjustments to more easily call IUPACs without splitting the output consensus VCF
        - Remvoed depth filtering here and just marking lower depth sites due to how freebayes handles reads and depth that
            differs from bcftools/samtools and causes some inconsistencies when masking low depth regions
        - Adjusting how Del+Snp complex sites are handled
        - TSV output for easier downstream visualiazation and reporting
        - Code adjustments and comments for why decisions were made and keeping track of variables better
'''

import argparse
import pysam
import csv
import re
from typing import Optional, Tuple, Generator

# Set for iupac assignment
iupac_map = {
    frozenset(['A']): 'A',
    frozenset(['T']): 'T',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['A', 'C']): 'M',
    frozenset(['A', 'G']): 'R',
    frozenset(['A', 'T']): 'W',
    frozenset(['C', 'G']): 'S',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}


def yield_alt_base(cigar: str, alt: str) -> Generator:
    """Expand CIGAR string and alt to yield alt the alt base for each reference position (handling indels)

    Params
    ------
        cigar (str): Variant CIGAR string
        alt (str): Alt allele string

    Returns
    -------
        Generator: alt base for each ref position
    """
    parts = re.findall(r'(\d+)([MXID])', cigar)
    expanded_cigar = ''.join(int(n) * op for n, op in parts)

    alt_pos = 0
    for cigar_str in expanded_cigar:
        if cigar_str == 'I':
            alt_pos += 1
        elif cigar_str == 'D':
            yield '-'
        else:
            yield alt[alt_pos]
            alt_pos += 1


def calculate_vafs(record: pysam.VariantRecord) -> list:
    """Calculate the variant allele fraction for each alt allele using freebayes' read/alt observation tags

    Params
    ------
        record (VariantRecord): Single variant record to parse

    Returns
    -------
        list: Containing VAFs for each alt allele
    """
    vafs = list()
    for i in range(0, len(record.alts)):
        alt_reads = int(record.info["AO"][i])
        vaf = round(float(alt_reads) / float(record.info["DP"]), 4)
        vafs.append(vaf)
    return vafs


def make_simple_record(vcf_header: pysam.VariantHeader, parent_record: pysam.VariantRecord,
                       position: int, ref: str, alt: str, vaf: float) -> pysam.VariantRecord:
    """Make a simple VCF record with the minimal information needed for the consensus vcf output

    Params
    ------
        vcf_header (VariantHeader): Simple VCF header for consensus generation
        parent_record (VariantRecord): Detailed parent record
        position (int): Genome position
        ref (str): Reference base
        alt (str): Alt base
        vaf (float): Variant allele fraction calculated

    Returns
    -------
        VariantRecord: Simple positional variant record for consensus vcf output
    """
    r = vcf_header.new_record() # blank header
    r.chrom = parent_record.chrom
    r.pos = position
    r.ref = ref
    r.alts = [ alt ]
    r.qual = parent_record.qual
    r.info["DP"] = parent_record.info["DP"]
    r.info["VAF"] = vaf
    return r


def base_max(vaf_by_base: dict, skip=None) -> Optional[str]:
    """Return the base with the highest value in VAF that is not skipped with optional skip param

    Params
    ------
        vaf_by_base (dict): Dict containing bases (keys) and their VAF for the position
        skip (str): Optional Bases to skip the check for, default None

    Returns
    -------
        str: The base with the highest VAF, defaults to None if there weren't any
    """
    max_vaf = 0.0
    max_b = None
    for b in ["A", "T", "G", "C"]:
        if b != skip and vaf_by_base[b] > max_vaf:
            max_vaf = vaf_by_base[b]
            max_b = b
    return max_b


def handle_sub(vcf_header: pysam.VariantHeader, record: pysam.VariantRecord) -> list:
    """Process *substitution* variants found by freebayes into a variant that
    can be applied to the final consensus sequence

    Params
    ------
        vcf_header (VariantHeader): Simple VCF header for consensus generation
        record (VariantRecord): Position variant record containing all needed information

    Returns
    -------
        list: List of tuples containing the final VariantRecords for the position and dict of their base_frequencies
    """
    output = list()

    # This can handle multi-allelic MNPs and the typical case of a biallelic SNP
    sub_length = len(record.ref)
    vafs = calculate_vafs(record)

    # Calculate the VAF of each base at each position of the MNP
    base_frequency = list()
    for i in range(0, sub_length):
        base_frequency.append({ "A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})

    for alt, vaf, cigar in zip(record.alts, vafs, record.info['CIGAR']):
        for i, base in enumerate(yield_alt_base(cigar, alt)):
            if base == '-':
                continue
            base_frequency[i][base] += vaf

    # Construct output records
    for i in range(0, sub_length):
        # Choose base with highest frequency, skipping the reference
        #  That way we can identify mixed sites where the reference is the most represented base
        #  and that are above our minimum frequency
        max_b = base_max(base_frequency[i], record.ref[i])
        if max_b is None:
            continue
        r = make_simple_record(vcf_header, record, record.pos + i, record.ref[i], max_b, base_frequency[i][max_b])

        # Make sure the reference base is accounted for if we need to assign an IUPAC for the site
        total_vafs = sum(vafs)
        if total_vafs < 1:
            ref_freq = round(1 - total_vafs, 4)
            base_frequency[i][record.ref[i]] += ref_freq

        output.append((r, base_frequency[i]))

    return output


def handle_indel(vcf_header: pysam.VariantHeader, record: pysam.VariantRecord, min_indel_threshold: float) -> list:
    """Process *indel* variants found by freebayes into a variant that can be applied to the final consensus sequence

    Params
    ------
        vcf_header (VariantHeader): Simple VCF header for consensus generation
        record (VariantRecord): Position variant record containing all needed information

    Returns
    -------
        list: List of tuples containing the final VariantRecords for the position and a blank dict to match the handle_snp output
    """
    output = list()
    vafs = calculate_vafs(record)

    # Special case, if we have evidence for multiple possible indels (eg CTTT -> C, CTTT -> CT)
    #  we decide whether to apply an indel based on the summed VAF across all alt alleles, then
    #  apply the most frequent ALT. This is because there is evidence for /an/ indel but it is
    #  ambiguous which one. We can't represent ambiguous indels in a consensus fasta so this
    #  is the best we can do.
    # Also note that we can get indels called below the threshold if there is an alt allele SNP and an indel at the same point
    #  Probably better as if we just went with reference base we'd be completely inaccurate
    #  But it may be a revisit to mask the site
    if sum(vafs) <= min_indel_threshold:
        return output

    # Argmax without bringing in numpy
    best_vaf_idx = None
    max_vaf = 0.0
    for idx, value in enumerate(vafs):
        if value > max_vaf:
            max_vaf = value
            best_vaf_idx = idx

    # Have to add the sub evaluation if the prevalent mutation is a SNP instead of an indel
    #  Otherwise we have incorrect results
    if len(record.ref) == len(record.alts[best_vaf_idx]):
        output = handle_sub(vcf_header, record)
    else:
        r = make_simple_record(vcf_header, record, record.pos, record.ref, record.alts[best_vaf_idx], [ max_vaf ])

        # Have to add in the Genotype for bcftools 1.20 to apply the variant
        #  So as were only filtering to 1 indel genotype use that
        #  SNPS end up slightly different to account for iupac codes
        r.samples[0]['GT'] = (1,)

        output.append((r, {}))
    return output


def get_base_code(base_frequencies: dict, consensus_ambiguity_threshold: float) -> Tuple[str, frozenset]:
    """Determine which base(s) are above the 1 - consensus ambiguity threshold and return single IUPAC code

    Params
    ------
        base_frequencies (dict): Dict containing the float frequencies of each base for the position
        ambiguity_threshold (float): Ambiguity threshold to include base as significant for position

    Returns
    -------
        Tuple: String IUPAC code along with the frozenset of significant bases identified
    """
    # Find bases with frequencies above consensus ambiguity threshold
    threshold = 1 - consensus_ambiguity_threshold
    significant_bases = {k for k, v in base_frequencies.items() if v >= threshold}

    # Look up the Consensus or IUPAC code for the set of significant bases
    return iupac_map.get(frozenset(significant_bases), 'N'), frozenset(significant_bases)


def init_parser() -> argparse.ArgumentParser:
    """Init Parser and parse CL args"""
    description = 'Process an output freebayes .vcf file to create summary VCF and TSV files'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-v', '--variants-output', required=True,
            help="The output file name for the filtered variants with no MNP splitting and added VAF calculations (VCF)")

    parser.add_argument('-c', '--consensus-sites-output', required=True,
            help="The output file name for split MNP variants and IUPAC variants that will be applied to generate the consensus sequence (VCF)")

    parser.add_argument('-t', '--tsv-sites-output', required=True,
            help="The output file name for consensus variants ormatted as a TSV file for later reporting")

    parser.add_argument('-d', '--mark-depth', type=int, default=10,
            help="Minimum depth to mark a variant as low depth. Depth filtering should be done with the Freebayes argument")

    parser.add_argument('-l', '--lower-ambiguity-frequency', type=float, default=0.25,
            help="Variants with frequency less than -l will be discarded entirely")

    parser.add_argument('-u', '--upper-ambiguity-frequency', type=float, default=0.75,
            help="Substitution variants with frequency less than -u will be encoded with IUPAC ambiguity codes based on underlying variation")

    parser.add_argument('-m', '--minimum-indel-threshold', type=float, default=0.60,
            help="Indel variants with frequency less than the -m threshold will be skipped")

    parser.add_argument('-q', '--min-quality', type=int, default=20,
            help="Minimum quality to call a variant")

    parser.add_argument('-n', '--no-frameshifts', action="store_true",
            help="Skip indel mutations that are not divisible by 3")

    parser.add_argument('invcf', action='store', nargs=1)

    return parser


def main() -> None:
    """Main entry point"""
    parser = init_parser()
    args = parser.parse_args()

    # Load input VCF file and get header info
    vcf = pysam.VariantFile(open(args.invcf[0], 'r'))
    out_header = vcf.header

    # Open the filtered only output VCF file to write unmodified accepted variants along with the addition of their VAFs
    out_header.info.add("VAF", number="A", type='Float', description="Variant allele fraction, called from observed reference/alt reads")
    filtered_variants_out = pysam.VariantFile(args.variants_output, 'w', header=out_header)

    # Open the consensus VCF that will have split sites, add the genotypes to call IUPACs and contain minimal information for the consensus sequence generation
    #  This includes an additional 2 headers in the VCF file
    out_header.info.add("ConsensusTag", number=1, type='String', description="The type of base included in the consensus sequence (ambiguous/consensus/indel)")
    out_header.info.add("ConsensusBase", number=1, type='String', description="The base included in the consensus sequence (ATGC or IUPAC)")
    consensus_variants_out = pysam.VariantFile(args.consensus_sites_output, 'w', header=out_header)

    # Setup list of dicts for later TSV creation and reporting
    tsv_data_list = []

    # Parsing VCF records to assign final consensus variants and filter out poor variant calls for the filtered VCF
    for base_record in vcf:
        # Determine if any allele in the variant is an indel to handle it properly
        has_indel = False
        for i in range(0, len(base_record.alts)):
            has_indel = has_indel or len(base_record.ref) != len(base_record.alts[i])

        # Get the variant depth just to mark if it is lower than our threshold
        #  The pipeline freebayes minimum depth is set to match the other minimum depth parts
        #  Freebayes calculates depth in a way that isn't easy to reproduce so we are going to flag
        #  these sites but not drop them here and instead mask them downstream if they do have too few reads
        # And then the True/False is for tracking
        depth = base_record.info["DP"]
        low_depth = False
        if depth < args.mark_depth:
            low_depth = True

        # Process the input variant record to handle multi-allelic variants and MNPs
        out_records = list()
        if has_indel:
            # indels need to be handled specially as we can't apply ambiguity codes to them
            out_records = handle_indel(out_header, base_record, args.minimum_indel_threshold)
        else:
            out_records = handle_sub(out_header, base_record)

        # Classify variants using VAF cutoffs for IUPAC ambiguity codes, quality cutoffs, etc.
        #  For out_tuple, its formatted as: ( record, base_frequency dict )
        accept_variant = False
        for out_tuple in out_records:
            final_record = out_tuple[0]

            # Should have resolved MNP sites to just be 1
            assert(len(final_record.alts) == 1)
            is_indel = len(final_record.ref) != len(final_record.alts[0])

            ### Filters for Poor Variants ###
            # Add final VAF value
            vaf = round(final_record.info["VAF"][0], 4)

            # Discard low frequency variants lower than our minimum ambiguity threshold
            if vaf < args.lower_ambiguity_frequency:
                continue

            # Discard low quality sites as recommended by freebayes
            if base_record.qual < args.min_quality:
                continue

            # Discard fs indels if provided the arg
            if is_indel:
                base_change = len(final_record.alts[0]) - len(final_record.ref)
                if base_change % 3 != 0 and args.no_frameshifts:
                    continue

            ### Setup basic TSV/VCF output requirements ###
            tsv_tag = "PASS"
            consensus_tag = "None"
            consensus_base = final_record.alts[0]
            genotype = (1,)

            # High-frequency subs and indels are always applied without ambiguity
            #  We don't have to do an indel VAF check here as it is dealt with in handle_indel
            if vaf > args.upper_ambiguity_frequency or is_indel:
                consensus_tag = "consensus"
                if is_indel:
                    consensus_tag = "indel"
            # Otherwise data is between the upper and lower frequency so setup and apply IUPAC
            else:
                tsv_tag = "AMBIGUOUS"
                consensus_tag = "ambiguous"
                iupac_base, fzset = get_base_code(out_tuple[1], args.upper_ambiguity_frequency)

                # Set genotype for bcftools with the `-I` arg to properly use
                #  `-I` will apply an IUPAC based on the given genotype (ex. (0,1) will give an IUPAC based on the ref and alt)
                if final_record.ref in fzset:
                    genotype = tuple(range(0, len(fzset)))
                    # Can't manip fzset inplace so rewrite
                    #  Need to take the reference out as we add the fzset bases back into the alt
                    tmp_set = fzset - {final_record.ref}
                    fzset = tmp_set
                else:
                    genotype = tuple(range(1, len(fzset) + 1))

                # Setting up to write to the consensus VCF file and to track what the base is
                final_record.alts = fzset
                vaf = sum([v for k, v in out_tuple[1].items() if k in fzset])
                consensus_base = iupac_base


            # Output for consensus generation and reporting
            final_record.info["ConsensusTag"] = consensus_tag
            final_record.info["ConsensusBase"] = consensus_base
            final_record.samples[0]['GT'] = genotype
            consensus_variants_out.write(final_record)

            # Setting up for TSV output for later reporting
            #  Format: chrom, pos, ref, alt, qual, depth, vaf, tag
            if low_depth:
                tsv_tag = "LOWDEPTH"
            tsv_data_list.append([
                final_record.chrom,
                final_record.pos,
                final_record.ref,
                ','.join(final_record.alts),
                int(final_record.qual),
                depth,
                vaf,
                tsv_tag
            ])

            # Only 1 has to be accepted to add the full original variant to the final filtered VCF
            #  But only that 1 good variant will be in the consensus
            accept_variant = True

        # Write the base record with the added VAF calculation to the final filtered VCF file
        if accept_variant:
            base_record.info["VAF"] = calculate_vafs(base_record)
            filtered_variants_out.write(base_record)

    # Write to TSV at the end
    headers = ['chrom', 'pos', 'ref', 'alt', 'qual', 'depth', 'vaf', 'tag']
    with open(args.tsv_sites_output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(tsv_data_list)

if __name__ == "__main__":
    main()

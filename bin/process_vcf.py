#!/usr/bin/env python
# Written by @jts from https://github.com/jts/ncov2019-artic-nf/blob/be26baedcc6876a798a599071bb25e0973261861/bin/process_gvcf.py

# Adjustments made such that were focused on just the mutations and not the GVCF info
#  Along with that, added in genotype to allow new versions of bcftools consensus to work
#  And adjusting how Del+Snp complex sites are handled

import argparse
import pysam
import csv

# Set for iupac assignment
iupac_map = {
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

def calculate_vafs(record):
    '''Calculate the variant allele fraction for each alt allele using freebayes' read/alt observation tags'''
    vafs = list()
    for i in range(0, len(record.alts)):
        alt_reads = int(record.info["AO"][i])
        vaf = round(float(alt_reads) / float(record.info["DP"]), 4)
        vafs.append(vaf)
    return vafs

def make_simple_record(vcf_header, parent_record, position, ref, alt, vaf):
    '''Make a simple VCF record with the minimal information needed to make the consensus sequence'''
    r = vcf_header.new_record()
    r.chrom = parent_record.chrom
    r.pos = position
    r.ref = ref
    r.alts = [ alt ]
    r.qual = parent_record.qual
    r.info["DP"] = parent_record.info["DP"]
    r.info["VAF"] = vaf
    return r

def base_max(vaf_by_base, skip=None):
    '''
    Return the base with the highest value in vaf_by_base
    Optionally skipping a character (eg. the reference base)
    '''
    max_vaf = 0.0
    max_b = None
    for b in "ACGT":
        if b != skip and vaf_by_base[b] > max_vaf:
            max_vaf = vaf_by_base[b]
            max_b = b
    return max_b

def handle_sub(vcf_header, record):
    '''
    Process substitution variants found by freebayes into a variant that can be applied to the
    final consensus sequence
    '''
    output = list()

    # this code is general enough to handle multi-allelic MNPs
    # and the typical case of a biallelic SNP
    sub_length = len(record.ref)

    vafs = calculate_vafs(record)

    # calculate the VAF of each base at each position of the MNP
    base_frequency = list()
    indel_vafs = 0
    for i in range(0, sub_length):
        base_frequency.append({ "A":0.0, "C":0.0, "G":0.0, "T":0.0 })

    for alt, vaf in zip(record.alts, vafs):
        # For rare case when a SNP is the most prevalent but there is also an idel
        if len(alt) != sub_length:
            indel_vafs += vaf
            continue
        for i,b in enumerate(alt):
            base_frequency[i][b] += vaf

    # construct output records
    for i in range(0, sub_length):

        # choose base with highest frequency, skipping the reference
        max_b = base_max(base_frequency[i], record.ref[i])
        if max_b is None:
            continue
        r = make_simple_record(vcf_header, record, record.pos + i, record.ref[i], max_b, base_frequency[i][max_b])

        # add Reference % post record for IUPAC tracking
        alt_freq_sum = sum(base_frequency[i].values()) + indel_vafs
        if alt_freq_sum < 1:
            ref_freq = round(1 - alt_freq_sum, 4)
            base_frequency[i][record.ref[i]] += ref_freq # For if it doesn't add to 1, don't want to overwrite the ref

        output.append((r, base_frequency[i]))
    return output

def handle_indel(vcf_header, record, min_indel_threshold):
    '''
    Process indel variants found by freebayes into a variant that should be
    applied to the consensus sequence
    '''
    output = list()
    vafs = calculate_vafs(record)

    # special case, if we have evidence for multiple possible indels (eg CTTT -> C, CTTT -> CT)
    # we decide whether to apply an indel based on the summed VAF across all alt alleles, then
    # apply the most frequent ALT. This is because there is evidence for /an/ indel but it is
    # ambiguous which one. We can't represent ambiguous indels in a consensus fasta so this
    # is the best we can do.
    if sum(vafs) <= min_indel_threshold:
        return output

    # argmax without bringing in numpy
    max_idx = None
    max_vaf = 0.0
    for idx, value in enumerate(vafs):
        if value > max_vaf:
            max_vaf = value
            max_idx = idx

    # Have to add the sub evaluation if the prevalent mutation is a SNP instead of an indel
    #  Otherwise we have incorrect results
    if len(record.ref) == len(record.alts[max_idx]):
        output = handle_sub(vcf_header, record)
    else:
        r = make_simple_record(vcf_header, record, record.pos, record.ref, record.alts[max_idx], [ max_vaf ])

        # Have to add in atleast the Genotype for bcftools 1.20 to apply the variant
        #  So as were only filtering to 1 genotype use that
        #  SNPS slightly different to account for iupac codes
        r.samples[0]['GT'] = (1,)

        output.append((r, {}))
    return output

def get_base_code(base_dict, upper_ambiguity):
    # Filter bases with value above 1 - threshold
    threshold = 1 - upper_ambiguity
    significant_bases = {k for k, v in base_dict.items() if v >= threshold}
    print(significant_bases, base_dict, upper_ambiguity)

    # Consensus
    if len(significant_bases) == 1:
        return significant_bases.pop()

    # Look up the IUPAC code for the set of significant bases
    return iupac_map.get(frozenset(significant_bases), 'N')

def main():
    '''Main entry point'''
    description = 'Process a .gvcf file to create a file of consensus variants, low-frequency variants and a coverage mask'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-v', '--variants-output', required=True,
            help=f"The output file name for variants (VCF records)\n")

    parser.add_argument('-c', '--consensus-sites-output', required=True,
            help=f"The output file name for variants that will be applied to generate the consensus sequence\n")

    parser.add_argument('-t', '--tsv-sites-output', required=True,
            help=f"The output file name for consensus variants that is formatted as a TSV file for reporting later\n")

    parser.add_argument('-d', '--min-depth', type=int, default=10,
            help=f"Minimum depth to call a variant")

    parser.add_argument('-l', '--lower-ambiguity-frequency', type=float, default=0.25,
            help=f"Variants with frequency less than -l will be discarded")

    parser.add_argument('-u', '--upper-ambiguity-frequency', type=float, default=0.75,
            help=f"Substitution variants with frequency less than -u will be encoded with IUPAC ambiguity codes")

    parser.add_argument('-m', '--minimum-indel-threshold', type=float, default=0.60,
            help=f"Indel variants with frequency less than the -m threshold will be skipped")

    parser.add_argument('-q', '--min-quality', type=int, default=20,
            help=f"Minimum quality to call a variant")

    parser.add_argument('-n', '--no-frameshifts', action="store_true",
            help=f"Skip indel mutations that are not divisible by 3")

    parser.add_argument('file', action='store', nargs=1)

    args = parser.parse_args()

    # Load VCF
    vcf = pysam.VariantFile(open(args.file[0],'r'))

    out_header = vcf.header

    # Open the output file with the filtered variant sites
    out_header.info.add("VAF", number="A", type='Float', description="Variant allele fraction, called from observed reference/alt reads")
    variants_out = pysam.VariantFile(args.variants_output, 'w', header=out_header)

    # Open the output file with the changes to apply to the consensus fasta
    # This includes an additional 2 headers in the VCF file
    out_header.info.add("ConsensusTag", number=1, type='String', description="The type of base to be included in the consensus sequence (ambiguous or consensus)")
    # Set it to 1 as the Multi-allelic sites should be resolved
    out_header.info.add("ConsensusBase", number=1, type='String', description="The base included in the consensus sequence (to track IUPACs mostly)")
    consensus_sites_out = pysam.VariantFile(args.consensus_sites_output, 'w', header=out_header)

    # Setup TSV data list for later reporting
    tsv_data_list = []

    for record in vcf:

        # Determine if any allele in the variant is an indel
        has_indel = False
        for i in range(0, len(record.alts)):
            has_indel = has_indel or len(record.ref) != len(record.alts[i])

        # Get the variant depth
        depth = record.info["DP"]
        low_depth = False

        # Ignore records that don't meet our minimum depth other than for multi-allelic sites
        if depth < args.min_depth:
            # For multi-allelic variants only, the freebayes depth may not be exactly positional
            #  So let them pass and mask them later if the variant position is below the given minimum depth
            #  This helps resolve on the cut-off depths especially in the MF-NCR
            if not any([x for x in record.alts if len(x) > 1]):
                continue
            low_depth = True

        # process the input variant record to handle multi-allelic variants and MNPs
        out_records = list()
        if has_indel:
            # indels need to be handle specially as we can't apply ambiguity codes
            out_records = handle_indel(out_header, record, args.minimum_indel_threshold)
        else:
            out_records = handle_sub(out_header, record)

        # Classify variants using VAF cutoffs for IUPAC ambiguity codes, etc
        #  For out_tuple, its record, base_frequency dict (for IUPAC)
        accept_variant = False
        for out_tuple in out_records:
            out_r = out_tuple[0]
            # at this point we should have resolved multi-allelic variants
            assert(len(out_r.alts) == 1)

            # Add final VAF
            vaf = out_r.info["VAF"][0]
            is_indel = len(out_r.ref) != len(out_r.alts[0])

            # Recheck that an indel is not low-depth if the site was multi-allelic initially
            #  This should help avoid errant indels in the MF-NCR
            if (is_indel) and (depth < args.min_depth):
                continue

            # Discard low frequency variants
            if vaf < args.lower_ambiguity_frequency:
                continue

            # Discard fs indels if provided the arg
            if is_indel and len(out_r.alts[0]) % 3 != 0 and args.no_frameshifts:
                continue

            # Discard low quality sites as recommended by freebayes
            #  Based on the data nothing really lower than the default min_qual unless it is low low quality
            if record.qual < args.min_quality:
                # Tracking for TSV only
                tsv_data_list.append([
                out_r.chrom,
                out_r.pos,
                out_r.ref,
                out_r.alts[0],
                int(out_r.qual),
                depth,
                vaf,
                'Fail'
                ])
                continue

            # Write a tag describing what to do with the variant
            consensus_tag = "None"
            consensus_base = out_r.alts[0]
            genotype = (1,)

            # High-frequency subs and indels are always applied without ambiguity
            # we don't have to do an indel VAF check here as it is dealt with in handle_indel
            if vaf > args.upper_ambiguity_frequency or is_indel:
                # always apply these to the consensus
                consensus_tag = "consensus"
                tsv_tag = "Passing Alt Allele Fraction"
                if vaf >= 0.90:
                    tsv_tag = "High Alt Allele Fraction"
                elif is_indel:
                    tsv_tag = "Passing Indel"
            else:
                # To capture IUPACs in reports easier have a separate column
                iupac_base = get_base_code(out_tuple[1], args.upper_ambiguity_frequency)

                # Record ambiguous SNPs in the consensus sequence with IUPAC codes
                consensus_tag = "ambiguous"
                # Genotype needs to be mixed to get an iupac if that is what the base should be
                if iupac_base not in ['A', 'T', 'G', 'C']:
                    genotype = (0,1)
                    tsv_tag = f"Ambiguous Base - {iupac_base}"
                else:
                    tsv_tag = "Passing Low Alt Allele Fraction"

            # Output for consensus generation and reporting
            out_r.info["ConsensusTag"] = consensus_tag
            out_r.info["ConsensusBase"] = consensus_base
            out_r.samples[0]['GT'] = genotype
            consensus_sites_out.write(out_r)

            # Setting up for TSV output for later reporting
            #  Format: chrom, pos, ref, alt, qual, depth, vaf, tag
            if low_depth:
                tsv_tag = "Low Depth"
            tsv_data_list.append([
                out_r.chrom,
                out_r.pos,
                out_r.ref,
                consensus_base,
                int(out_r.qual),
                depth,
                vaf,
                tsv_tag
            ])

            accept_variant = True

        if accept_variant:
            record.info["VAF"] = calculate_vafs(record)
            variants_out.write(record)

    # Write to TSV at the end
    headers = ['chrom', 'pos', 'ref', 'alt', 'qual', 'depth', 'vaf', 'tag']
    with open(args.tsv_sites_output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(tsv_data_list)

if __name__ == "__main__":
    main()

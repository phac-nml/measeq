#!/usr/bin/env python
# Written by @jts from https://github.com/jts/ncov2019-artic-nf/blob/be26baedcc6876a798a599071bb25e0973261861/bin/process_gvcf.py
# May want to adjust/update it but for now it is unchanged

import argparse
import pysam

# calculate the variant allele fraction for each alt allele using freebayes' read/alt observation tags
def calculate_vafs(record):
    vafs = list()
    total_depth = float(record.info["DP"])
    for i in range(0, len(record.alts)):
        alt_reads = int(record.info["AO"][i])
        vaf = float(alt_reads) / float(record.info["DP"])
        vafs.append(vaf)
    return vafs

# make a simple VCF record with the minimal information needed to make the consensus sequence
def make_simple_record(vcf_header, parent_record, position, ref, alt, vaf):
    r = vcf_header.new_record()
    r.chrom = parent_record.chrom
    r.pos = position
    r.ref = ref
    r.alts = [ alt ]
    r.qual = parent_record.qual
    r.info["DP"] = parent_record.info["DP"]
    r.info["VAF"] = vaf
    return r

# process indel variants found by freebayes into a variant that should be
# applied to the consensus sequence
def handle_indel(vcf_header, record):
    output = list()
    vafs = calculate_vafs(record)

    # special case, if we have evidence for multiple possible indels (eg CTTT -> C, CTTT -> CT)
    # we decide whether to apply an indel based on the summed VAF across all alt alleles, then
    # apply the most frequent ALT. This is because there is evidence for /an/ indel but it is
    # ambiguous which one. We can't represent ambiguous indels in a consensus fasta so this
    # is the best we can do.
    if sum(vafs) < 0.5:
        return output

    # argmax without bringing in numpy
    max_idx = None
    max_vaf = 0.0
    for idx, value in enumerate(vafs):
        if value > max_vaf:
            max_vaf = value
            max_idx = idx
    
    r = make_simple_record(vcf_header, record, record.pos, record.ref, record.alts[max_idx], [ max_vaf ])

    # Have to add in atleast the Genotype for bcftools 1.20 to apply the variant
    #  So as were only filtering to 1 genotype use that
    #  SNPS slightly different to account for iupac codes
    r.samples[0]['GT'] = (1,)

    output.append(r)
    return output

# return the base with the highest value in vaf_by_base,
# optionally skipping a character (eg. the reference base)
def base_max(vaf_by_base, skip=None):
    max_vaf = 0.0
    max_b = None
    for b in "ACGT":
        if b != skip and vaf_by_base[b] > max_vaf:
            max_vaf = vaf_by_base[b]
            max_b = b
    return max_b

def handle_sub(vcf_header, record):
    output = list()

    # this code is general enough to handle multi-allelic MNPs
    # and the typical case of a biallelic SNP
    sub_length = len(record.ref)

    vafs = calculate_vafs(record)

    # calculate the VAF of each base at each position of the MNP
    base_frequency = list()
    for i in range(0, sub_length):
        base_frequency.append({ "A":0.0, "C":0.0, "G":0.0, "T":0.0 })

    for alt, vaf in zip(record.alts, vafs):
        assert(len(alt) == sub_length)
        for i,b in enumerate(alt):
            base_frequency[i][b] += vaf

    # construct output records
    for i in range(0, sub_length):

        # choose base with highest frequency, skipping the reference
        max_b = base_max(base_frequency[i], record.ref[i])
        if max_b is None:
            continue
        r = make_simple_record(vcf_header, record, record.pos + i, record.ref[i], max_b, base_frequency[i][max_b])

        output.append(r)
    return output

def main():

    description = 'Process a .gvcf file to create a file of consensus variants, low-frequency variants and a coverage mask'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-v', '--variants-output', required=True,
            help=f"The output file name for variants (VCF records)\n")

    parser.add_argument('-c', '--consensus-sites-output', required=True,
            help=f"The output file name for variants that will be applied to generate the consensus sequence\n")

    parser.add_argument('-d', '--min-depth', type=int, default=10,
            help=f"Minimum depth to call a variant")

    parser.add_argument('-l', '--lower-ambiguity-frequency', type=float, default=0.25,
            help=f"Variants with frequency less than -l will be discarded")

    parser.add_argument('-u', '--upper-ambiguity-frequency', type=float, default=0.75,
            help=f"Substitution variants with frequency less than -u will be encoded with IUPAC ambiguity codes")
    
    parser.add_argument('-q', '--min-quality', type=int, default=20,
            help=f"Minimum quality to call a variant")
    
    parser.add_argument('-n', '--no-frameshifts', action="store_true",
            help=f"Skip indel mutations that are not divisible by 3")

    parser.add_argument('file', action='store', nargs=1)

    args = parser.parse_args()

    # Load VCF
    vcf = pysam.VariantFile(open(args.file[0],'r'))

    # Initalize depth mask to all zeros for all contigs
    for r in vcf.header.records:
        if r.type == "CONTIG":
            genome_length = int(r['length'])

    out_header = vcf.header

    # Open the output file with the filtered variant sites
    out_header.info.add("VAF", number="A", type='Float', description="Variant allele fraction, called from observed reference/alt reads")
    variants_out = pysam.VariantFile(args.variants_output, 'w', header=out_header)

    # Open the output file with the changes to apply to the consensus fasta
    # This includes an additional tag in the VCF file
    out_header.info.add("ConsensusTag", number=1, type='String', description="The type of base to be included in the consensus sequence (ambiguous or fixed)")
    consensus_sites_out = pysam.VariantFile(args.consensus_sites_output, 'w', header=out_header)

    for record in vcf:

        # Set depth for this part of the genome
        # this works for both gVCF blocks and regular variants
        # because pos/stop are set appropriately
        v_start = record.pos
        v_end = record.stop
        depth = record.info["DP"]

        # Do nothing with Records that don't meet our minimum depth
        if depth < args.min_depth:
            continue

        # Determine if any allele in the variant is an indel
        has_indel = False
        for i in range(0, len(record.alts)):
            has_indel = has_indel or len(record.ref) != len(record.alts[i])

        # process the input variant record to handle multi-allelic variants and MNPs
        out_records = list()
        if has_indel:
            # indels need to be handle specially as we can't apply ambiguity codes
            out_records = handle_indel(out_header, record)
        else:
            out_records = handle_sub(out_header, record)

        # classify variants using VAF cutoffs for IUPAC ambiguity codes, etc
        accept_variant = False
        for out_r in out_records:

            # at this point we should have resolved multi-allelic variants
            assert(len(out_r.alts) == 1)

            vaf = out_r.info["VAF"][0]
            is_indel = len(out_r.ref) != len(out_r.alts[0])

            # Discard low frequency variants
            if vaf < args.lower_ambiguity_frequency:
                continue

            # Discard fs indels if provided the arg
            if is_indel and len(out_r.alts[0]) % 3 != 0 and args.no_frameshifts:
                continue

            # Discard low quality sites as recommended by freebayes
            #  Might need to add a proper calculation here for it based on depth but based on data
            #  nothing really is this low unless its very mixed or low low depth
            if record.qual < args.min_quality:
                print(out_r.qual)
                continue

            # Write a tag describing what to do with the variant
            consensus_tag = "None"
            genotype = (1,)

            # high-frequency subs and indels are always applied without ambiguity
            # we don't have to do an indel VAF check here as it is dealt with in handle_indel
            if vaf > args.upper_ambiguity_frequency or is_indel:
                # always apply these to the consensus
                consensus_tag = "fixed"
            else:
                # record ambiguous SNPs in the consensus sequence with IUPAC codes
                consensus_tag = "ambiguous"
                # Genotype needs to be mixed to get an iupac
                genotype = (0,1)
            out_r.info["ConsensusTag"] = consensus_tag
            out_r.samples[0]['GT'] = genotype
            consensus_sites_out.write(out_r)
            accept_variant = True

        if accept_variant:
            record.info["VAF"] = calculate_vafs(record)
            variants_out.write(record)

if __name__ == "__main__":
    main()

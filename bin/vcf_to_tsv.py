#!/usr/bin/env python3
'''
Simple script to parse clair3 VCF file to create a tsv file
'''

import argparse
import vcf

def init_parser() -> argparse.ArgumentParser:
    '''
    Specify command line arguments
    Returns command line parser with inputs
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        '--sample',
        required=True,
        type=str,
        help='Sample name'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        type=str,
        help='Path to artic Passing VCF file'
    )
    parser.add_argument(
        '-f',
        '--fail_vcf',
        type=str,
        help='Path to artic Failing VCF'
    )
    parser.add_argument(
        '-o',
        '--outfile',
        required=False,
        type=str,
        help='Name of outfile. Default is "SAMPLE.consensus.tsv"'
    )
    parser.add_argument(
        '--annotated',
        required=False,
        default=False,
        action='store_true',
        help='Specify if the vcf is annotated with snpeff or not'
    )
    return parser

def process_variants_details(var, annotated: bool, variant_type: str, variants_analyzed=[]) -> dict:
    """
    Process variant and create dictionary for entry
    Remove duplicated variants
    """
    # Check if seen before - if not add it to list and continue on
    variant_str = f'{var.REF}{var.POS}{var.ALT[0]}' # Only first alt allele, shouldn't have more than one with the process currently
    if variant_str in variants_analyzed:
        return {}
    variants_analyzed.append(variant_str)

    info = 'Full Alignment'
    if 'P' in var.INFO:
        info = 'Pileup'

    # Non-annotated we just need the info in TSV format
    # Tag is for final report, moreso for illumina but needed here
    if not annotated:
        out = {
            'chrom': str(var.CHROM),
            'pos': int(var.POS),
            'ref': str(var.REF),
            'alt': str(var.ALT[0]),
            'qual': str(var.QUAL),
            'info': info,
            'depth': int(var.samples[0].data.DP),
            'vaf': round(var.samples[0].data.AF[0], 4),
            'tag': variant_type
        }
        return out

    # VCF ANN is formatted as a "," list to start so if we find it, we take the first entry and split on |
    #  Also want to make sure that we don't have duplicates by checking if our variant str has been seen already
    var_ann = var.INFO.get('ANN', '')
    if not var_ann:
        out = {
            'chrom': str(var.CHROM),
            'pos': int(var.POS),
            'variant': '',
            'consequence': '',
            'ref': str(var.REF),
            'alt': str(var.ALT[0]),
            'qual': str(var.QUAL),
            'info': info,
            'depth': int(var.samples[0].data.DP),
            'vaf': round(var.samples[0].data.AF[0]),
            'tag': variant_type
        }
        return out
    ann_list = var_ann[0].split('|')

    # Create output detail dict based on VCFannotation format position
    #  https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
    impact = ann_list[2]
    gene = ann_list[3]
    prot = ann_list[10]
    pro_var = f'{gene}'
    if prot:
        pro_var = f'{gene} {prot}'
    consequence = ann_list[1]
    out = {
        'chrom': str(var.CHROM),
        'pos': int(var.POS),
        'variant': pro_var,
        'consequence': consequence,
        'ref': str(var.REF),
        'alt': str(var.ALT[0]),
        'qual': str(var.QUAL),
        'info': info,
        'depth': int(var.samples[0].data.DP),
        'vaf': round(var.samples[0].data.AF[0]),
        'tag': variant_type
    }
    return out

def process_variant_file(vcfin: str, annotated: bool, var_type: str) -> list:
    """Parse given VCF to return a list of variant dicts"""
    vcf_reader = vcf.Reader(filename=vcfin)
    variants = []

    # To fix out of range issue from pyvcf on empty files
    try:
        for var in vcf_reader:
            var_dict = process_variants_details(var, annotated, var_type)
            if var_dict != {}:
                variants.append(var_dict)
    except IndexError:
        pass
    return variants

def write_outfile(outfile: str, variants: list) -> None:
    """Write output file to given location"""
    columns = variants[0].keys()
    header = '\t'.join(columns)
    sorted_variants = sorted(variants, key=lambda x: x['pos'])
    with open(outfile, 'w') as f:
        f.write(header)
        f.write('\n')
        for var_dict in sorted_variants:
            f.write('\t'.join([str(var_dict[col]) for col in columns]))
            f.write('\n')

def main() -> None:
    """Main process"""
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Load in VCF and parse
    variants = process_variant_file(args.vcf, args.annotated, 'Consensus')
    failing_variants = process_variant_file(args.fail_vcf, args.annotated, 'Masked')

    variants.extend(failing_variants)

    # Write outfile
    outfile = f'{args.sample}.consensus.tsv'
    if args.outfile:
        outfile = args.outfile

    # Create empty file with correct headers if no variants
    if variants == []:
        base_cols = [
            'chrom',
            'pos',
            'variant',
            'consequence',
            'ref',
            'alt',
            'qual',
            'info',
            'depth',
            'vaf',
            'tag',
        ]
        if not args.annotated:
            base_cols = [col for col in base_cols if col not in ['variant', 'consequence']]

        variants = [{key: '' for key in base_cols}]

    write_outfile(outfile, variants)

if __name__ == '__main__':
    main()

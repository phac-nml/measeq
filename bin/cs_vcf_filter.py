#!/usr/bin/env python3
'''
Filter VCF Variants originally from https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcf_filter.py

Notable Changes:
  - Added a custom parameter to adjust the QUAL threshold with it set at 8 by default.
  - Added a filter for RefCall bases.
  - Depth based kick out for overlapping areas where there is a variant that doesn't make sense
'''

from cyvcf2 import VCF, Writer
from collections import defaultdict
import subprocess

def create_depth_map(bam):
    """
    Creating a depth map with samtools structured based on Genomic position as {Pos: Depth}
    """
    depth_map = {}
    with subprocess.Popen(
        ['samtools', 'depth', '-a', bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        bufsize=1
    ) as process:
        for line in process.stdout:
            # chrom - pos - depth
            pos_data = line.strip().split('\t')
            depth_map[int(pos_data[1])] = int(pos_data[2])
    return depth_map

def in_frame(v):
    if len(v.ALT) > 1:
        print("This code does not support multiple genotypes!")
        raise SystemExit
    ref = v.REF
    alt = v.ALT[0]
    bases = len(alt) - len(ref)
    if not bases:
        return True
    if bases % 3 == 0:
        return True
    return False

class Clair3Filter:
    def __init__(self, no_frameshifts, min_depth, min_variant_quality):
        self.no_frameshifts = no_frameshifts
        self.min_depth = min_depth
        self.min_variant_quality = min_variant_quality
        self.min_frameshift_quality = 30
        self.min_allele_frequency = 0.60

    def check_filter(self, v):
        qual = v.QUAL

        if qual < self.min_variant_quality:
            return False

        # Filter out low allele frequency variants
        try:
            allele_freq = v.format("AF")[0][0]
        except Exception:
            print(
                f"ERROR: Could not find AF for variant at {v.CHROM}:{v.POS}, cannot filter on allele frequency"
            )
            raise SystemExit(1)

        if allele_freq < self.min_allele_frequency:
            return False

        if not in_frame(v):
            if self.no_frameshifts:
                return False
            # require a higher quality for frameshifting indels, they're far more likely to be errors
            if qual < self.min_frameshift_quality:
                return False

        try:
            # We don't really care about the depth here, just skip it if it isn't there
            depth = v.INFO["DP"]
        except KeyError:
            depth = v.format("DP")[0][0]

        if depth < self.min_depth:
            return False

        return True

def go(args):
    vcf_reader = VCF(args.inputvcf)
    vcf_writer = Writer(args.output_pass_vcf, vcf_reader, "w")
    vcf_writer.write_header()
    vcf_writer_filtered = Writer(args.output_fail_vcf, vcf_reader, "w")
    vcf_writer_filtered.write_header()
    filter = Clair3Filter(args.no_frameshifts, args.min_depth, args.min_qual)

    depth_map = create_depth_map(args.inputbam)

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = "%s-%s" % (v.CHROM, v.POS)
        group_variants[indx].append(v)

    for v in variants:

        # Quick pre-filter to remove rubbish that we don't want adding to the mask
        try:
            if v.format("DP")[0][0] <= 1:
                print(f"Suppress variant {v.POS} due to low depth")
                continue

            # If there is a massive depth difference in an amplicon overlap area ignore
            #  We'd expect ~50/50 for an overlap with normalization so 10% would be very low
            #  And a real variant would be in both pools so if one was lower than the other that'd add to it here
            if v.format("DP")[0][0] < args.min_threshold_depth * depth_map[v.POS]:
                print(f"Supress variant at {v.POS} due to depth difference of {v.format('DP')[0][0]} to {depth_map[v.POS]}")
                continue

        except KeyError:
            pass

        # Completely skip RefCalls
        if v.ALT == []:
            print(f"skipping RefCall at {v.POS}")
            continue

        # now apply the filter to send variants to PASS or FAIL file
        if filter.check_filter(v):
            vcf_writer.write_record(v)
        else:
            variant_passes = False

            indx = "%s-%s" % (v.CHROM, v.POS)
            if len(group_variants[indx]) > 1:
                for check_variant in group_variants[indx]:
                    if filter.check_filter(check_variant):
                        variant_passes = True

            if not variant_passes:
                vcf_writer_filtered.write_record(v)

            else:
                print("Suppress variant %s" % (v.POS))


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-frameshifts", action="store_true")
    parser.add_argument("--min-depth", type=int)
    parser.add_argument("--min-qual", type=int, default=8)
    parser.add_argument("--min-threshold-depth", type=float, default=0.05)
    parser.add_argument("inputvcf")
    parser.add_argument("inputbam")
    parser.add_argument("output_pass_vcf")
    parser.add_argument("output_fail_vcf")

    args = parser.parse_args()

    go(args)


if __name__ == "__main__":
    main()

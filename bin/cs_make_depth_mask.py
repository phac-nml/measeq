#!/usr/bin/env python3
'''
Script to create a genomic depth mask for consensus generation
    Adapted from https://github.com/artic-network/fieldbioinformatics/blob/master/artic/make_depth_mask.py
    to remove the need for RG tags steps for non-amplicon data
'''

from Bio import SeqIO
import argparse
import itertools
import os
import pysam

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=20,
                        help="Max depth at which a position will be masked. Default: 20")
    parser.add_argument('--ignore-deletions', action="store_true", default=False,
                        help="if set, positional depth counts will ignore reads with reference deletions\
                        (i.e. evaluates positional depths on ref matches, not read span")
    parser.add_argument('reference')
    parser.add_argument('bamfile')
    parser.add_argument('outfile')
    return parser

# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    """Create list of intervals with sequential numbers"""
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def collect_depths(bamfile: str, refName: str, minDepth: int, ignoreDeletions: bool) -> list:
    """Collect read depth of coverage per reference position in a BAM file.

    Parameters
    ----------
    bamfile : string
        The BAM file that needs processing

    refName : string
        The name of the reference sequence to collect the depths for

    minDepth : int
        The minimum depth to report coverage for (0 will be reported if coverage < minDepth at a given position)

    ignoreDeletions : bool
        If true, positional depth counts will ignore reads with reference deletions

    Returns
    -------
    list
        Index is the reference position, value is the corresponding coverage depth
    """
    # check the BAM file exists
    if not os.path.exists(bamfile):
        raise Exception("bamfile doesn't exist (%s)" % bamfile)

    # open the BAM file
    bamFile = pysam.AlignmentFile(bamfile, 'rb')

    # get the TID for the reference
    tid = bamFile.get_tid(refName)
    if tid == -1:
        raise Exception(
            "bamfile does not contain specified reference (%s)" % refName)

    # create a depth vector to hold the depths at each reference position
    depths = [0] * bamFile.get_reference_length(refName)

    # generate the pileup
    for pileupcolumn in bamFile.pileup(refName, max_depth=10000, truncate=False, min_base_quality=0):

        # process the pileup column
        for pileupread in pileupcolumn.pileups:
            # process the pileup read
            if pileupread.is_refskip:
                continue
            if pileupread.is_del:
                if not ignoreDeletions:
                    depths[pileupcolumn.pos] += 1
            elif not pileupread.is_del:
                depths[pileupcolumn.pos] += 1
            else:
                raise Exception("unhandled pileup read encountered")

        # if final depth for pileup column < minDepth, report 0 and update the mask_vector
        if depths[pileupcolumn.pos] < minDepth:
            depths[pileupcolumn.pos] = 0

    # close file and return depth vector
    bamFile.close()
    return depths

def main():
    """Main entry point"""
    # Args
    parser = init_parser()
    args = parser.parse_args()

    # Open the reference sequence and collect the sequence header and sequence length of the first record
    record = list(SeqIO.parse(args.reference, "fasta"))[0]
    seqID = record.id
    seqLength = len(record.seq)

    # Collect the depths from the pileup, replacing any depth<minDepth with 0
    try:
        depths = collect_depths(args.bamfile, seqID,
                                    args.depth, args.ignore_deletions)
    except Exception as e:
        print(e)
        raise SystemExit(1)

    # Check the number of positions in the reported depths matches the reference sequence
    if len(depths) != seqLength:
        print("pileup length did not match expected reference sequence length")

    # Create a mask_vector that records reference positions where depth < minDepth
    mask_vector = []
    for pos, depth in enumerate(depths):
        if depth == 0:
            mask_vector.append(pos)

    # Get the intervals from the mask_vector
    intervals = list(intervals_extract(mask_vector))

    # Create the mask outfile
    maskfh = open(args.outfile, 'w')
    for i in intervals:
        maskfh.write("%s\t%s\t%s\n" % (seqID, i[0]+1, i[1]+1))
    maskfh.close()

if __name__ == "__main__":
    main()

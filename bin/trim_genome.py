#!/usr/bin/env python3
"""
Trim input FASTA file by aligning to an UTR-trimmed genome with minimap2 and calculating
the coordinates where it aligns to pass to samtools faidx

Example diagram:

    input_seq        (target): 1 |-----------------------------| 15894
    trim_guide        (query):     1 |---------------------| 15678
    alignment  (query coords):        30 |-------------| 15400
    alignment (target coords):       130 |-------------| 15500
    unaligned (for placement):       |---|             |---|
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


@dataclass(frozen=True)
class PafAlignment:
    """Dataclass to hold the Minimap2 PAF output
    Query being the trimmed reference, target the input
    """
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapping_quality: int

    @classmethod
    def from_line(cls, line: str):
        """Parse and validate a PAF alignment line"""
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 12:
            raise ValueError("PAF records must contain at least 12 tab-separated fields")

        try:
            alignment = cls(
                query_name=fields[0],
                query_length=int(fields[1]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target_name=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                matches=int(fields[9]),
                block_length=int(fields[10]),
                mapping_quality=int(fields[11])
            )

        except (ValueError, IndexError) as e:
            raise ValueError(f"Could not parse PAF record: {line}") from e

        return alignment


def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser(
        description="Calculate and extract a final trimmed genome based on a minimap2 alignment"
    )

    parser.add_argument(
        "-t",
        "--trim_guide",
        type=Path,
        required=True,
        help="Already trimmed FASTA file to act as a guide for trimming input genome"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Input FASTA genome to trim based on alignment with the trim guide"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=Path,
        default=Path("trimmed_out.fasta"),
        help="Path to output trimmed FASTA file (default: trimmed_out.fasta)"
    )
    parser.add_argument(
        "--preset",
        default="asm10",
        choices=[
            "asm5",
            "asm10",
            "asm20",
            "lr:hq"
        ],
        help="minimap2 alignment preset to use for alignment (default: asm10)"
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=60,
        help="Number of bases per output FASTA line (default: 60)"
    )

    return parser


def run_command(command: list, tool: list) -> str:
    """Run given command with subprocess.run and return string stdout"""
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout

    except subprocess.CalledProcessError as e:
        details = (e.stderr or "").strip()
        message = f"{tool} failed with exit code {e.returncode}."

        if details:
            message += f"\n{details}"
        raise RuntimeError(message) from e


def find_best_alignment(alignment_paf: str) -> PafAlignment:
    """Select best minimap2 alignment using the length.
    There should only be 1 unless the input fasta had large issues
    """
    alignments = []

    # Each line is a tab-separated PAF output so use that to populate
    #  a list of PafAlignment dataclasses
    for line in alignment_paf.splitlines():
        if not line.strip():
            continue

        alignments.append(PafAlignment.from_line(line))

    if not alignments:
        print("Minimap2 produced no alignments. Exiting with no output")
        sys.exit(0)

    # Just in case there are more than 1 alignments, select the largest
    best = max(
        alignments,
        key=lambda alignment: (alignment.query_end - alignment.query_start)
    )

    return best


def calculate_trim_coordinates(alignment: PafAlignment) -> Tuple[str, int, int]:
    """Calculate where the input (target) should be trimmed too based off of
    the trim_guide sequence boundry (query in minimap2 PAF output).
    Query is mapped to the target to get placement, then we calculate
    the unaligned start and end to get the whole trimmed region we need
    """
    # Determine number of unaligned bases at the start and end of the genome to properly place region
    missing_aln_start = alignment.query_start
    missing_aln_end = alignment.query_length - alignment.query_end

    # Positive strand means left and right are start and end respectively
    #  Adding in 1 to get to samtools 1-index based input from PAF 0-index output
    if alignment.strand == "+":
        raw_query_start = alignment.target_start - missing_aln_start + 1
        raw_query_end = alignment.target_end + missing_aln_end
    # We shouldn't get negative strand but if we do we swap the overhangs
    else:
        raw_query_start = alignment.target_start - missing_aln_end + 1
        raw_query_end = alignment.target_end + missing_aln_start

    # Determine the boundry but make sure we stay inside the target coordinates
    #  So smallest is starting at 1 and longest is the length of the input query sequence
    start = max(1, raw_query_start)
    end = min(alignment.target_length, raw_query_end)

    if start > end:
        raise ValueError(f"Calculated interval {alignment.target_name}:{start}-{end} is invalid as the start is after the end")

    # Region argument for samtools to extract
    region = f"{alignment.target_name}:{start}-{end}"

    return region, start, end


def extract_region(input_fasta, region, expected_length, line_width):
    """Using samtools faidx and the calculated region, extract a FASTA formatted string of
    the region from the input_fasta file and validate that it matches the target length
    """
    cmd = [
        "samtools",
        "faidx",
        "-n",
        str(line_width),
        str(input_fasta),
        region,
    ]
    fasta_text = run_command(cmd, "samtools faidx")

    # Validate we have a nt sequence output and for measles that it starts with ATG
    #  as we are trimming the UTRs at the start and end of the genome
    nt_sequence = "".join(
        line.strip()
        for line in fasta_text.splitlines()
        if line and not line.startswith(">")
    )

    if len(nt_sequence) != expected_length:
        raise ValueError(
            f"Expected {expected_length} bp from {region}, but samtools returned {len(nt_sequence)} bp"
        )

    # What we can start with, its just the start codon with increasing Ns incase the region has no data
    valid_starts = {"ATG", "NTG", "NNG", "NNN"}

    if nt_sequence[:3] not in valid_starts:
        raise ValueError(
            f"Expected nt sequence to start with a start codon or Ns, got {nt_sequence[:3]}"
        )

    if not nt_sequence:
        raise ValueError(f"No sequence was extracted for {region}.")

    return fasta_text


def main():
    """Main entry point to the script"""
    parser = init_parser()

    args = parser.parse_args()

    for path in (args.trim_guide, args.input):
        if not path.is_file():
            raise FileNotFoundError(f"Input FASTA not found: {path}")

    # Find the best alignment to extract the region
    minimap2_result = run_command(
        [
            "minimap2",
            "-x",
            args.preset,
            str(args.input),
            str(args.trim_guide)
        ],
        "minimap2"
    )
    alignment = find_best_alignment(minimap2_result)

    # Calculate trimming coordinates
    region, start, end = calculate_trim_coordinates(alignment)
    expected_length = end - start + 1

    # Print quickly for log tracking if there is ever an issue
    print(alignment)
    print(f"Expected trimmed length: {expected_length}")

    # Create final fasta file
    fasta_text = extract_region(
        args.input,
        region,
        expected_length,
        args.line_width,
    )

    with open(args.outfile, "w") as f:
        f.write(fasta_text)


if __name__ == "__main__":
    main()

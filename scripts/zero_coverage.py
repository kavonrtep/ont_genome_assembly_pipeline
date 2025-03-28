#!/usr/bin/env python
"""
This script extracts regions with zero coverage from an input BAM file and writes them to a BED file.
For each reference contig, it iterates over pileup columns (which are only produced for positions with at least one read)
and then infers the gaps (i.e. positions with zero coverage) as the intervals between consecutive covered positions.
If a contig has no coverage, the entire contig is output.

Usage:
    python zero_coverage.py --bam input.bam --out zero_coverage.bed [--min_base_quality 0]
"""

import pysam
import argparse

def find_zero_coverage_regions(bam_path, bed_out, min_base_quality=0):
    bam = pysam.AlignmentFile(bam_path, "rb")
    with open(bed_out, "w") as bedfh:
        # Iterate over each reference contig.
        for contig in bam.references:
            contig_length = bam.get_reference_length(contig)
            prev = -1  # previous covered position
            found_coverage = False
            # Iterate over pileup columns for this contig.
            for column in bam.pileup(contig, min_base_quality=min_base_quality, truncate=True):
                found_coverage = True
                pos = column.pos  # 0-indexed
                # If there's a gap between previous covered position and current column, output the gap.
                if pos > prev + 1:
                    start = prev + 1
                    end = pos
                    bedfh.write(f"{contig}\t{start}\t{end}\n")
                prev = pos
            # If no coverage was found on the contig, output the entire contig.
            if not found_coverage:
                bedfh.write(f"{contig}\t0\t{contig_length}\n")
            # Check for a gap at the end of the contig.
            if found_coverage and prev < contig_length - 1:
                start = prev + 1
                end = contig_length
                bedfh.write(f"{contig}\t{start}\t{end}\n")
    bam.close()

def main():
    parser = argparse.ArgumentParser(
        description="Extract zero coverage regions from a BAM file and output as a BED file."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file (must be indexed)")
    parser.add_argument("--out", required=True, help="Output BED file for zero coverage regions")
    parser.add_argument("--min_base_quality", type=int, default=0,
                        help="Minimum base quality for pileup (default: 0)")
    args = parser.parse_args()
    find_zero_coverage_regions(args.bam, args.out, args.min_base_quality)
    print(f"Zero coverage regions written to {args.out}")

if __name__ == "__main__":
    main()

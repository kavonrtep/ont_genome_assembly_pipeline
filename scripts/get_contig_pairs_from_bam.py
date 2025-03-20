#!/usr/bin/env python
"""
Script to evaluate all possible connections between contigs based on read mappings.

For each read in the BAM file, the script:
  - Iterates over all alignments (using pysam).
  - Filters out alignments with query_alignment_length below a specified threshold.
  - Checks if the alignment falls within a window (default 100kb) at either the 5' or 3' end
    of the contig. (An alignment may qualify for both.)
  - Groups valid alignments by read.
  - For each read with at least two valid alignments on different contigs, it outputs every
    pairing with:
       read_id, contig1, contig2, connection_type (e.g. "3->5"), aln_size1, aln_size2.

Usage:
  python contig_connection.py --bam myfile.bam --min_aln_length 50 --window_size 100000 --out connections.tsv
"""

import pysam
import argparse
import os

def get_contig_lengths(bamfile):
    """Extract contig lengths from the BAM header."""
    lengths = {}
    for sq in bamfile.header.get("SQ", []):
        lengths[sq["SN"]] = int(sq["LN"])
    return lengths

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate contig connections via read mappings at contig ends."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file (indexed)")
    parser.add_argument("--window_size", type=int, default=100000,
                        help="Window size (bp) at contig ends to consider (default: 100kb)")
    parser.add_argument("--min_aln_length", type=int, default=10000,
                        help="Minimum alignment length to consider (default: 10kb)")
    parser.add_argument("--out", default="contig_connections.tsv",
                        help="Output file name (default: contig_connections.tsv)")
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.bam, "rb")
    print(f"Reading contig lengths from {args.bam}")
    contig_lengths = get_contig_lengths(bamfile)

    # Dictionary: read_id -> list of valid alignment entries.
    # Each entry is a dict with: "contig", "end" (either "5" or "3"),
    # "aln_length", and "qstart" (read coordinate for sorting).
    read_alignments = {}

    # Iterate over all alignments in the BAM file.
    # Using fetch(until_eof=True) ensures we process the entire file.
    print("Getting read alignments...")
    for aln in bamfile.fetch(until_eof=True):
        if aln.is_unmapped:
            continue
        if aln.query_alignment_length < args.min_aln_length:
            continue

        contig = aln.reference_name
        if contig not in contig_lengths:
            continue
        contig_length = contig_lengths[contig]

        valid_end_types = []
        # Check if alignment falls in the 5' window (0-based start < window_size)
        if aln.reference_start < args.window_size:
            valid_end_types.append("5")
        # Check if alignment falls in the 3' window:
        # Since pysam returns reference_end as 0-based exclusive,
        # if contig_length - aln.reference_end < window_size then it is in the 3' window.
        if (contig_length - aln.reference_end) < args.window_size:
            valid_end_types.append("3")
        if not valid_end_types:
            continue

        read_id = aln.query_name
        # Create an entry for each end type the alignment qualifies for.
        for end_type in valid_end_types:
            entry = {
                "contig": contig,
                "end": end_type,
                "aln_length": aln.query_alignment_length,
                "qstart": aln.query_alignment_start
            }
            read_alignments.setdefault(read_id, []).append(entry)

    bamfile.close()
    print("Writing contig connections to", args.out)
    # Open output file and write header.
    with open(args.out, "w") as out_f:
        out_f.write("read_id\tcontig1\tcontig2\tconnection_type\taln_size_1\taln_size_2\n")
        # For each read, if it has at least two valid alignments from different contigs,
        # sort by the read's alignment start coordinate (qstart) and output each pair.
        for read_id, alns in read_alignments.items():
            # Only consider reads with at least 2 alignments.
            if len(alns) < 2:
                continue
            # Sort alignments by query alignment start.
            alns.sort(key=lambda x: x["qstart"])
            n = len(alns)
            for i in range(n):
                for j in range(i + 1, n):
                    # Report only if the two alignments are from different contigs.
                    if alns[i]["contig"] == alns[j]["contig"]:
                        continue
                    connection_type = f"{alns[i]['end']}->{alns[j]['end']}"
                    out_f.write(f"{read_id}\t{alns[i]['contig']}\t{alns[j]['contig']}\t"
                                f"{connection_type}\t{alns[i]['aln_length']}\t{alns[j]['aln_length']}\n")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import subprocess
import csv
from Bio import SeqIO


def merge_intervals(intervals):
    """Merge overlapping intervals.
    intervals: list of (start, end) tuples.
    Returns a list of merged intervals.
    """
    if not intervals:
        return []
    # Sort intervals by start position
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev = merged[-1]
        # If current interval overlaps or touches previous, merge them.
        # (Using <= ensures that contiguous intervals are merged)
        if current[0] <= prev[1] + 1:
            merged[-1] = (prev[0], max(prev[1], current[1]))
        else:
            merged.append(current)
    return merged


def run_blast(query_fasta, blast_db, num_threads, word_size, raw_output_file):
    """Run BLASTn in tabular mode with  saving the output to raw_output_file.
    The BLAST output format is 6 with columns:
    qseqid, pident, length, qstart, qend, qlen.
    """
    blast_cmd = ["blastn", "-query", query_fasta, "-db", blast_db, "-out",
        raw_output_file, "-outfmt", "6 qseqid pident length qstart qend qlen",
        "-num_threads", str(num_threads), "-word_size", str(word_size)]
    print("Running BLAST command: " + " ".join(blast_cmd))
    result = subprocess.run(
        blast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        raise RuntimeError("BLAST failed: " + result.stderr)
    print("BLAST results saved to", raw_output_file)


def parse_blast_file(raw_output_file, min_coverage, min_identity):
    """
    Parse BLAST tabular output from file.
    For each query, merge non-overlapping HSP intervals, cap total aligned length to qlen,
    and compute weighted average identity.
    Returns a list of tuples:
       (read_id, qlen, merged_aligned_length, coverage%, mean_identity%)
    """
    # Dictionary: read_id -> {'qlen': int, 'intervals': list of (start, end),
    #                          'weighted_identity_sum': float, 'total_aln_length': int}
    results = {}
    with open(raw_output_file, "r") as infile:
        for line in infile:
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            qseqid, pident, aln_length, qstart, qend, qlen = fields
            pident = float(pident)
            # Convert positions to integers
            aln_length = int(aln_length)
            qstart = int(qstart)
            qend = int(qend)
            qlen = int(qlen)
            # Initialize entry for read if needed.
            if qseqid not in results:
                results[qseqid] = {'qlen': qlen, 'intervals': [],
                    'weighted_identity_sum': 0.0, 'total_aln_length': 0}
            # Get start and end (ensure start <= end)
            start = min(qstart, qend)
            end = max(qstart, qend)
            # Save the interval (assuming inclusive coordinates)
            results[qseqid]['intervals'].append((start, end))
            # Accumulate identity weighted by this interval's length
            interval_length = end - start + 1
            results[qseqid]['weighted_identity_sum'] += pident * interval_length
            results[qseqid]['total_aln_length'] += interval_length

    final_results = []
    for read_id, data in results.items():
        merged_intervals = merge_intervals(data['intervals'])
        total_aligned = sum(end - start + 1 for start, end in merged_intervals)
        # Cap total aligned length at the query length to avoid overcounting
        if total_aligned > data['qlen']:
            total_aligned = data['qlen']
        coverage = (total_aligned / data['qlen']) * 100
        # Calculate weighted mean identity using total_aln_length from BLAST (not the merged one)
        # This may overestimate if overlaps exist, but it's acceptable if we filter by coverage.
        mean_identity = (data['weighted_identity_sum'] / data['total_aln_length']) if \
        data['total_aln_length'] > 0 else 0.0
        if coverage >= min_coverage and mean_identity >= min_identity:
            final_results.append(
                (read_id, data['qlen'], total_aligned, round(coverage, 2),
                 round(mean_identity, 2))
                )
    return final_results


def main():
    parser = argparse.ArgumentParser(
        description="Filter Nanopore reads for plastid DNA contamination using BLAST."
        )
    parser.add_argument(
        "-i", "--input", required=True, help="Input FASTA file of Nanopore reads"
        )
    parser.add_argument(
        "-d", "--db", required=True,
        help="BLAST database containing plastid DNA sequences"
        )
    parser.add_argument(
        "-o", "--output", required=True, help="Output tab-delimited file"
        )
    parser.add_argument(
        "--min_coverage", type=float, default=90.0,
        help="Minimum coverage of the read required (in percentage)"
        )
    parser.add_argument(
        "--min_identity", type=float, default=90.0,
        help="Minimum weighted average identity (in percentage)"
        )
    parser.add_argument(
        "--num_cpu", type=int, default=1, help="Number of CPU threads for BLAST"
        )
    parser.add_argument(
        "--word_size", type=int, default=11, help="Word size parameter for BLAST"
        )
    parser.add_argument(
        "--raw_blast", default="output.raw_blast",
        help="Filename to save raw BLAST output (tabular)"
        )
    args = parser.parse_args()

    # Run BLAST and save output to file
    run_blast(args.input, args.db, args.num_cpu, args.word_size, args.raw_blast)

    # Parse the saved BLAST output file
    print("Parsing BLAST output from", args.raw_blast)
    filtered_results = parse_blast_file(
        args.raw_blast, args.min_coverage, args.min_identity
        )

    # Write filtered results to output file (tab-delimited)
    with open(args.output, "w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(
            ["ReadID", "QueryLength", "AlignedLength", "Coverage(%)", "MeanIdentity(%)"]
            )
        for row in filtered_results:
            writer.writerow(row)
    print("Done. Filtered results written to", args.output)


if __name__ == "__main__":
    main()

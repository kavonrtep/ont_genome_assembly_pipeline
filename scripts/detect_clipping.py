#!/usr/bin/env python
"""
This script scans an input BAM file for soft-clipping events that might indicate assembly problems.
For each read (primary or supplementary) that has a soft clip at either end with a length at least a threshold
(default: 10kb), the script records an event with:
  - contig name (seqname)
  - clipping boundary (position: for left clip, read.reference_start; for right clip, read.reference_end)
  - clipping side ("left" or "right")
  - read ID
Optionally, it can compute the per-base coverage at that position (--add_coverage).

The script then groups (merges) nearby clipping events (within a tolerance, default: 10 bp) to produce BED intervals.
Finally, for each merged interval that has at least 2 clipped reads, the script defines an extended window
(+/- 10kb by default) and counts the number of reads spanning (i.e. covering) the original interval.
The final BED files (one for left and one for right clipping) have columns:
  seqname, start, end, score, name
where "score" is the number of clipped reads in the merged interval and "name" is the number of spanning reads.
"""

import pysam
import argparse
from collections import defaultdict
import re


def get_clip_events(bam_path, min_clip):
    """
    Iterate over all alignments and collect soft-clipping events with clip length >= min_clip.
    Returns a list of events (dicts) with keys: 'contig', 'position', 'side', and 'read_id'.
    """
    events = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                continue
            # do not use alignment with low mapping quality (MAPQ == 0)
            if read.mapping_quality == 0:
                continue
            # alignment length must be at least min_clip
            if read.query_alignment_length < min_clip:
                continue
            # Check left-end soft clipping (first CIGAR operation; op code 4)
            op, length = read.cigartuples[0]
            if op == 4 and length >= min_clip:
                event = {
                    "contig": read.reference_name,
                    "position": read.reference_start,
                    "side": "left",
                    "read_id": read.query_name
                }
                events.append(event)
            # Check right-end soft clipping (last CIGAR operation)
            op, length = read.cigartuples[-1]
            if op == 4 and length >= min_clip:
                event = {
                    "contig": read.reference_name,
                    "position": read.reference_end,  # one past the last aligned base
                    "side": "right",
                    "read_id": read.query_name
                }
                events.append(event)
    return events

def write_events(events, output_path):
    """
    Write detailed clipping events to a TSV file with columns:
    seqname, position, clipping_side, coverage, read_id.
    """
    with open(output_path, "w") as outfh:
        outfh.write("seqname\tposition\tclipping_side\tcoverage\tread_id\n")
        for event in events:
            cov_val = event.get("coverage", "NA")
            outfh.write(f'{event["contig"]}\t{event["position"]}\t{event["side"]}\t{cov_val}\t{event["read_id"]}\n')


def get_clip_events(bam_path, min_clip):
    """
    Iterate over alignments and collect soft-clipping events with clip length >= min_clip.
    For each left- or right-clipped read, after detecting the clip, we parse its SA tag
    (if present) to collect a list of contigs that might represent an extension.

    For a left clip:
      - Use the primary alignment's start (read.reference_start) as the boundary.
      - Then, for each SA alignment, if it is on the same contig and its end (computed from its CIGAR)
        is <= primary alignment's start, consider it an extension.
      - If the SA alignment is on a different contig, add that contig.

    For a right clip:
      - Use the primary alignment's end (read.reference_end) as the boundary.
      - Then, for each SA alignment, if it is on the same contig and its start is >= primary alignment's end,
        consider it an extension.
      - If the SA alignment is on a different contig, add that contig.

    Returns a list of events (dicts) with keys: 'contig', 'position', 'side', 'read_id', and 'extension_contigs'.
    """
    events = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                continue

            # Process left soft clipping.
            left_op, left_len = read.cigartuples[0]
            if left_op == 4 and left_len >= min_clip:
                primary_start = read.reference_start
                ext_contigs = []
                if read.has_tag("SA"):
                    sa_tag = read.get_tag(
                        "SA"
                        )  # SA tag: "rname,pos,strand,CIGAR,mapQ,NM;" repeated.
                    for sa in sa_tag.strip(";").split(";"):
                        fields = sa.split(",")
                        if len(fields) < 6:
                            continue
                        sa_rname = fields[0]
                        try:
                            sa_pos = int(fields[1]) - 1  # Convert 1-based to 0-based.
                        except ValueError:
                            continue
                        sa_cigar = fields[3]
                        # Calculate reference length consumed by SA alignment.
                        sa_ref_len = sum(
                            int(n) for n, op in
                            re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar) if op in "MDN=X"
                            )
                        sa_start = sa_pos
                        sa_end = sa_pos + sa_ref_len
                        if sa_rname != read.reference_name:
                            ext_contigs.append(sa_rname)
                        else:
                            if sa_end <= primary_start:
                                ext_contigs.append(sa_rname)
                event = {"contig": read.reference_name, "position": primary_start,
                    "side": "left", "read_id": read.query_name,
                    "extension_contigs": ",".join(sorted(set(ext_contigs)))}
                events.append(event)

            # Process right soft clipping.
            right_op, right_len = read.cigartuples[-1]
            if right_op == 4 and right_len >= min_clip:
                primary_end = read.reference_end
                ext_contigs = []
                if read.has_tag("SA"):
                    sa_tag = read.get_tag("SA")
                    for sa in sa_tag.strip(";").split(";"):
                        fields = sa.split(",")
                        if len(fields) < 6:
                            continue
                        sa_rname = fields[0]
                        try:
                            sa_pos = int(fields[1]) - 1
                        except ValueError:
                            continue
                        sa_cigar = fields[3]
                        sa_ref_len = sum(
                            int(n) for n, op in
                            re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar) if op in "MDN=X"
                            )
                        sa_start = sa_pos
                        sa_end = sa_pos + sa_ref_len
                        if sa_rname != read.reference_name:
                            ext_contigs.append(sa_rname)
                        else:
                            if sa_start >= primary_end:
                                ext_contigs.append(sa_rname)
                event = {"contig": read.reference_name, "position": primary_end,
                    "side": "right", "read_id": read.query_name,
                    "extension_contigs": ",".join(sorted(set(ext_contigs)))}
                events.append(event)
    return events


def write_events(events, output_path):
    """
    Write detailed clipping events to a TSV file with columns:
    seqname, position, clipping_side, coverage, read_id, extension_contigs.
    If coverage was not computed, outputs "NA".
    """
    with open(output_path, "w") as outfh:
        outfh.write(
            "seqname\tposition\tclipping_side\tcoverage\tread_id\textension_contigs\n"
            )
        for event in events:
            cov_val = event.get("coverage", "NA")
            ext = event.get("extension_contigs", "")
            outfh.write(
                f'{event["contig"]}\t{event["position"]}\t{event["side"]}\t{cov_val}\t{event["read_id"]}\t{ext}\n'
            )

def add_coverage_to_events(events, bam_path):
    """
    For each event, compute the per-base coverage using pysam.count_coverage.
    If the event position is beyond the contig length (e.g. for right-end events), adjust it.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for event in events:
            contig = event["contig"]
            pos = event["position"]
            ref_len = bam.get_reference_length(contig)
            if pos >= ref_len:
                pos = ref_len - 1
            if pos < 0:
                pos = 0
            cov = bam.count_coverage(contig, pos, pos+1)
            total_cov = sum(c[0] for c in cov)
            event["coverage"] = total_cov
            event["position"] = pos
    return events


def merge_events_by_side(events, side, tolerance):
    """
    Group events (for a given side: 'left' or 'right') by contig and merge events that are within
    the specified tolerance (in bp). Returns a list of merged intervals (dicts) with keys:
    'contig', 'start', 'end', and 'count' (number of clipped reads merged).
    """
    grouped = defaultdict(list)
    for ev in events:
        if ev["side"] == side:
            grouped[ev["contig"]].append(ev["position"])
    merged_intervals = []
    for contig, positions in grouped.items():
        positions.sort()
        start = positions[0]
        end = positions[0]
        count = 1
        for pos in positions[1:]:
            if pos - end <= tolerance:
                end = pos
                count += 1
            else:
                merged_intervals.append({
                    "contig": contig,
                    "start": start,
                    "end": end + 1,  # BED format: half-open interval [start, end)
                    "count": count
                })
                start = pos
                end = pos
                count = 1
        merged_intervals.append({
            "contig": contig,
            "start": start,
            "end": end + 1,
            "count": count
        })
    return merged_intervals

def write_bedgraph(intervals, output_path, extra_col=False):
    """
    Write intervals to a BED file. If extra_col is True, assume intervals have a 'spanning'
    field and output columns: contig, start, end, score (clipped count), name (spanning count).
    Otherwise, output contig, start, end, score.
    """
    with open(output_path, "w") as outfh:
        for iv in intervals:
            if extra_col:
                outfh.write(f'{iv["contig"]}\t{iv["start"]}\t{iv["end"]}\t{iv["count"]}\t{iv["spanning"]}\n')
            else:
                outfh.write(f'{iv["contig"]}\t{iv["start"]}\t{iv["end"]}\t{iv["count"]}\n')

def add_spanning_info(intervals, bam_path, span_window, side):
    """
    For each merged interval in intervals, define an extended window of ±span_window bp around the interval,
    and count the number of reads in the BAM file that span the original interval and are not clipped
    at the same side. The count is stored in the 'spanning' field.
    """
    updated_intervals = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for iv in intervals:
            if iv["count"] < 2:
                continue
            # Define extended window
            window_start = max(0, iv["start"] - span_window)
            window_end = iv["end"] + span_window
            spanning_count = 0
            for read in bam.fetch(iv["contig"], window_start, window_end):
                if read.is_unmapped:
                    continue
                # Must span the original interval
                if read.reference_start > iv["start"] or read.reference_end < iv["end"]:
                    continue
                # Exclude reads that are clipped at the problematic end.
                if side == "left":
                    if read.cigartuples and read.cigartuples[0][0] == 4:
                        continue
                elif side == "right":
                    if read.cigartuples and read.cigartuples[-1][0] == 4:
                        continue
                spanning_count += 1
            iv["spanning"] = spanning_count
            updated_intervals.append(iv)
    return updated_intervals

def main():
    parser = argparse.ArgumentParser(
        description="Detect soft-clipping events from a BAM file, merge nearby events, and summarize spanning reads."
    )

    parser.add_argument("--bam", required=True, help="Input BAM file (filtered and indexed)")
    parser.add_argument("--out", required=True, help="Output file prefix")
    parser.add_argument("--min_clip", type=int, default=10000,
                        help="Minimum soft clipping length to report (default: 10000 "
                             "bp). Minimal alignment length must be at least this value.")
    parser.add_argument("--add_coverage", action="store_true", default=False,
                        help="Compute coverage at each clipping position")
    parser.add_argument("--merge_tolerance", type=int, default=10,
                        help="Tolerance (bp) for merging nearby clipping events (default: 10 bp)")
    parser.add_argument("--span_window", type=int, default=10000,
                        help="Window size (bp) for counting spanning reads around an interval (default: 10000 bp)")
    args = parser.parse_args()

    # Collect detailed clipping events.
    events = get_clip_events(args.bam, args.min_clip)
    if args.add_coverage:
        events = add_coverage_to_events(events, args.bam)
    detailed_out = f"{args.out}_clip_info.tsv"
    write_events(events, detailed_out)
    print(f"Written detailed clipping events to {detailed_out}")

    # Merge events separately for left and right.
    left_intervals = merge_events_by_side(events, "left", args.merge_tolerance)
    right_intervals = merge_events_by_side(events, "right", args.merge_tolerance)

    left_bedgraph = f"{args.out}_left_clip.bedgraph"
    right_bedgraph = f"{args.out}_right_clip.bedgraph"
    write_bedgraph(left_intervals, left_bedgraph)
    write_bedgraph(right_intervals, right_bedgraph)
    print(f"Written merged clipping events to {left_bedgraph} and {right_bedgraph}")

    # For each merged interval (with at least 2 clipped reads), count spanning reads in the extended window.
    left_intervals_spanning = add_spanning_info(left_intervals, args.bam,
                                                args.span_window, side="left")
    right_intervals_spanning = add_spanning_info(right_intervals, args.bam,
                                                 args.span_window, side="right")

    left_span_bed = f"{args.out}_left_span.bed"
    right_span_bed = f"{args.out}_right_span.bed"
    # Write BED files with 5 columns: seqname, start, end, score (clipped count), name (spanning count)
    write_bedgraph(left_intervals_spanning, left_span_bed, extra_col=True)
    write_bedgraph(right_intervals_spanning, right_span_bed, extra_col=True)
    print(f"Written spanning information to {left_span_bed} and {right_span_bed}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Iterative Greedy Merging of PAF intervals

This script reads a PAF file (assumed sorted by qname and qstart), groups
the alignments by query name, target name, and strand, and then performs
an iterative greedy merging based on a gap (offset) criterion. Two intervals
within a group are merged if:
  - They share the same qname, tname, and strand.
  - They appear in order (i.e. candidate.qstart >= current.qend).
  - The gap in query and target are similar, i.e.:
       Let delta_query = candidate.qstart - current.qend
       For plus strand: delta_target = candidate.tstart - current.tend
                      For minus strand: delta_target = current.tstart - candidate.tend
    Then if both deltas are ≥ 0 and the absolute difference between them is
    less than or equal to (mean_delta * tolerance_fraction), we merge.
The merging is done iteratively (all-to-all within each group) until no more
merges can be performed.

Chains with fewer intervals than `min_intervals` are discarded.

Usage:
    python iterative_merge_paf.py input.paf output.paf -t 20 -m 2

  where -t is the tolerance percentage (default 0 means exact gap match)
  and -m is the minimum number of intervals for a chain to be kept.
"""

import sys
import math
import argparse
from collections import defaultdict


def can_merge(interval1, interval2, tolerance_fraction):
    """
    Decide whether two intervals (assumed in order) can be merged.
    Both intervals must share the same qname, tname, and strand.

    For plus strand:
        delta_query = interval2['qstart'] - interval1['qend']
        delta_target = interval2['tstart'] - interval1['tend']
    For minus strand:
        delta_query = interval2['qstart'] - interval1['qend']
        delta_target = interval1['tstart'] - interval2['tend']
    They are mergeable if both deltas are >= 0 and:
        abs(|delta_query| - |delta_target|) <= (mean_delta * tolerance_fraction)
    """
    # Ensure order in query
    if interval2['qstart'] < interval1['qend']:
        return False

    delta_query = interval2['qstart'] - interval1['qend']
    strand = interval1['strand']
    if strand == '+':
        delta_target = interval2['tstart'] - interval1['tend']
    elif strand == '-':
        delta_target = interval1['tstart'] - interval2['tend']
    else:
        sys.stderr.write(f"Warning: Unrecognized strand '{strand}'.\n")
        return False

    if delta_query < 0 or delta_target < 0:
        return False

    # Compute mean delta and allowed difference
    mean_delta = (abs(delta_query) + abs(delta_target)) / 2.0
    if mean_delta > 1000000:
        return False
    # intervals width
    q_width = abs(interval1['qend'] - interval1['qstart'])
    t_width = abs(interval1['tend'] - interval1['tstart'])
    # merge only if mean delta in max 10x of the width of the intervals
    if mean_delta > 10 * max(q_width, t_width):
        return False
    diff = abs(abs(delta_query) - abs(delta_target))
    if mean_delta == 0 or diff <= mean_delta * tolerance_fraction:
        return True
    return False


def merge_intervals(interval1, interval2):
    """
    Merge two intervals (which are assumed mergeable) into a single interval.
    For plus strand, the union is taken:
        new_qstart = min(qstart1, qstart2)
        new_qend   = max(qend1, qend2)
        new_tstart = min(tstart1, tstart2)
        new_tend   = max(tend1, tend2)
    For minus strand, we follow the convention used in previous chaining:
        new_qstart = min(qstart1, qstart2)
        new_qend   = max(qend1, qend2)
        new_tstart = max(tstart1, tstart2)
        new_tend   = min(tend1, tend2)
    Other fields (nmatch, alen) are summed and mapq becomes the minimum.
    The interval count is summed.
    """
    new_interval = {}
    new_interval['qname'] = interval1['qname']
    new_interval['qlen'] = interval1['qlen']
    new_interval['tname'] = interval1['tname']
    new_interval['tlen'] = interval1['tlen']
    new_interval['strand'] = interval1['strand']

    new_interval['qstart'] = min(interval1['qstart'], interval2['qstart'])
    new_interval['qend'] = max(interval1['qend'], interval2['qend'])

    if new_interval['strand'] == '+':
        new_interval['tstart'] = min(interval1['tstart'], interval2['tstart'])
        new_interval['tend'] = max(interval1['tend'], interval2['tend'])
    else:  # for minus strand, follow previous chaining convention
        new_interval['tstart'] = max(interval1['tstart'], interval2['tstart'])
        new_interval['tend'] = min(interval1['tend'], interval2['tend'])

    new_interval['nmatch'] = interval1['nmatch'] + interval2['nmatch']
    new_interval['alen'] = interval1['alen'] + interval2['alen']
    new_interval['mapq'] = min(interval1['mapq'], interval2['mapq'])
    new_interval['intervals'] = interval1.get('intervals', 1) + interval2.get(
        'intervals', 1
        )

    return new_interval


def iterative_greedy_merge(intervals, tolerance_fraction):
    """
    Given a list of intervals (all with the same qname, tname, and strand),
    iteratively merge any pair that are mergeable until no further merging is possible.
    """
    merged = True
    iterations = 0
    while merged:
        iterations += 1
        print("Merging, iteration", iterations)
        merged = False
        new_intervals = []
        used = [False] * len(intervals)
        # We'll try a simple pairwise merging:
        for i in range(len(intervals)):
            if used[i]:
                continue
            current = intervals[i]
            for j in range(i + 1, len(intervals)):
                if used[j]:
                    continue
                candidate = intervals[j]
                if can_merge(current, candidate, tolerance_fraction):
                    current = merge_intervals(current, candidate)
                    used[j] = True
                    merged = True
            new_intervals.append(current)
        intervals = new_intervals
    return intervals


def iterative_greedy_merge(intervals, tolerance_fraction, max_scan=10):
    """
    Given a list of intervals (all with the same qname, tname, and strand),
    iteratively merge any pair that are mergeable until no further merging is possible.
    This version only scans at most max_scan candidate intervals for each interval.

    Parameters:
      intervals: list of interval dictionaries.
      tolerance_fraction: allowed percentage difference (as a fraction) for merging.
      max_scan: maximum number of subsequent intervals to consider for merging.

    Returns:
      A list of merged intervals.
    """
    merged = True
    iterations = 0
    while merged:
        iterations += 1
        print("Merging, iteration", iterations)
        merged = False
        new_intervals = []
        used = [False] * len(intervals)
        # We'll try a pairwise merge, but for each interval only scan the next max_scan intervals.
        for i in range(len(intervals)):
            if used[i]:
                continue
            current = intervals[i]
            # Limit the inner loop to only max_scan candidates
            scan_limit = min(len(intervals), i + 1 + max_scan)
            for j in range(i + 1, scan_limit):
                if used[j]:
                    continue
                candidate = intervals[j]
                if can_merge(current, candidate, tolerance_fraction):
                    current = merge_intervals(current, candidate)
                    used[j] = True
                    merged = True
            new_intervals.append(current)
        intervals = new_intervals
        # should be sorted by qstart, but just in case
        intervals.sort(key=lambda r: r['qstart'])
    return intervals


def merge_paf_intervals(input_file, output_file, tolerance_percent, min_intervals):
    """
    Read a PAF file (sorted by qname and qstart), group records by (qname, tname, strand),
    then perform iterative greedy merging (using an absolute tolerance on the gap offsets)
    within each group. Chains with fewer intervals than min_intervals are discarded.
    Merged chains are written to the output file.
    """
    tolerance_fraction = tolerance_percent / 100.0
    records = []
    print("Reading PAF records...")
    with open(input_file, 'r') as infile:
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 12:
                sys.stderr.write(
                    f"Warning: Line {line_num} has fewer than 12 fields. Skipping.\n"
                    )
                continue
            try:
                qname = fields[0]
                qlen = int(fields[1])
                qstart = int(fields[2])
                qend = int(fields[3])
                strand = fields[4]
                tname = fields[5]
                tlen = int(fields[6])
                tstart = int(fields[7])
                tend = int(fields[8])
                nmatch = int(fields[9])
                alen = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                sys.stderr.write(
                    f"Warning: Line {line_num} has invalid numeric fields. Skipping.\n"
                    )
                continue

            # Ensure coordinates are increasing.
            qstart, qend = min(qstart, qend), max(qstart, qend)
            tstart, tend = min(tstart, tend), max(tstart, tend)

            rec = {'qname': qname, 'qlen': qlen, 'qstart': qstart, 'qend': qend,
                'strand': strand, 'tname': tname, 'tlen': tlen, 'tstart': tstart,
                'tend': tend, 'nmatch': nmatch, 'alen': alen, 'mapq': mapq,
                'intervals': 1}
            records.append(rec)

    # Group records by (qname, tname, strand)
    groups = defaultdict(list)
    print("Creating groups and merging intervals...")
    for rec in records:
        key = (rec['qname'], rec['tname'], rec['strand'])
        groups[key].append(rec)

    # For each group, sort by qstart and apply iterative greedy merging.
    merged_records = []
    print("iterating over groups")
    for key, recs in groups.items():
        recs.sort(key=lambda r: r['qstart'])
        print("Merging group", key)
        merged_group = iterative_greedy_merge(recs, tolerance_fraction)
        for merged in merged_group:
            if merged['intervals'] >= min_intervals:
                merged_records.append(merged)

    # Sort final merged records by qname and qstart
    merged_records.sort(key=lambda r: (r['qname'], r['qstart']))

    with open(output_file, 'w') as outfile:
        for rec in merged_records:
            out_fields = [rec['qname'], str(rec['qlen']), str(rec['qstart']),
                str(rec['qend']), rec['strand'], rec['tname'], str(rec['tlen']),
                str(rec['tstart']), str(rec['tend']), str(rec['nmatch']),
                str(rec['alen']), str(rec['mapq'])]
            outfile.write('\t'.join(out_fields) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="""
        Iterative Greedy Merging of PAF alignments.
        The script groups alignments by qname, tname, and strand, then iteratively merges intervals
        if the gap (offset) between them is within the allowed tolerance (as a percentage of the mean gap).
        Alignments with fewer intervals than min_intervals are discarded.
        The PAF file should be sorted by qname and qstart before use.
        """
        )
    parser.add_argument("input_paf", help="Input PAF file (sorted by qname,qstart)")
    parser.add_argument(
        "output_paf", help="Output PAF file with merged/chained alignments"
        )
    parser.add_argument(
        "-t", "--tolerance", type=float, default=0.0,
        help="Tolerance percentage for gap differences (default: 0.0)"
        )
    parser.add_argument(
        "-m", "--min-intervals", type=int, default=1,
        help="Minimum number of merged intervals required to keep a chain (default: 1)"
        )
    args = parser.parse_args()
    print("Processing", args.input_paf)
    merge_paf_intervals(
        args.input_paf, args.output_paf, args.tolerance, args.min_intervals
        )


if __name__ == "__main__":
    main()

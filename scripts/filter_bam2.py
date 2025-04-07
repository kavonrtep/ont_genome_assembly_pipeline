#!/usr/bin/env python

import pysam
import argparse
import collections
import itertools
import subprocess
import tempfile
import os

def passes_filters(aln, max_mismatch_prop, min_mapq, min_aln_length):
    """
    Return True if the alignment passes filtering criteria.
    """
    if aln.is_unmapped:
        return False

    aln_length = aln.query_alignment_length
    if not aln_length or aln_length == 0:
        return False

    if min_aln_length is not None and aln_length < min_aln_length:
        return False

    if min_mapq is not None and aln.mapping_quality < min_mapq:
        return False

    if max_mismatch_prop is not None:
        try:
            nm = aln.get_tag("NM")
        except KeyError:
            return False  # if no NM tag, skip
        mismatch_prop = nm / float(aln_length)
        if mismatch_prop > max_mismatch_prop:
            return False

    return True

def update_read_group(aln_list):
    """
    For a given list of alignments for a read (all that passed filtering),
    update flags and tags so that if any non-supplementary alignments exist,
    one is designated as primary (the best by MAPQ) and the others are marked as secondary.
    Also update the NH (number of hits) and HI (hit index) tags accordingly.

    Returns the updated list of alignments (supplementary alignments are left unchanged).
    """
    non_supp = [aln for aln in aln_list if not aln.is_supplementary]

    if non_supp:
        primaries = [aln for aln in non_supp if not aln.is_secondary]
        if primaries:
            primary_aln = max(primaries, key=lambda aln: aln.mapping_quality)
        else:
            primary_aln = max(non_supp, key=lambda aln: aln.mapping_quality)
            primary_aln.flag &= ~0x100  # clear secondary flag

        for aln in non_supp:
            if aln != primary_aln:
                aln.flag |= 0x100

        new_nh = len(non_supp)
        for idx, aln in enumerate(non_supp, start=1):
            if aln.has_tag("NH"):
                aln.set_tag("NH", new_nh, value_type='i')
            if aln.has_tag("HI"):
                aln.set_tag("HI", idx, value_type='i')
    return aln_list

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline to filter BAM alignments. The input BAM is first sorted by read name "
                    "using samtools (with specified number of threads), then filtered, and the resulting "
                    "filtered BAM is coordinate sorted and indexed. All intermediate files are removed."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--out", required=True, help="Output BAM file (coordinate sorted filtered BAM)")
    parser.add_argument("--max_mismatch_prop", type=float, default=None,
                        help="Maximum allowed mismatch proportion (NM / alignment length)")
    parser.add_argument("--min_mapq", type=int, default=None, help="Minimum MAPQ value")
    parser.add_argument("--min_aln_length", type=int, default=None,
                        help="Minimum alignment length (query alignment length)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for samtools sort")
    args = parser.parse_args()

    if args.max_mismatch_prop is None and args.min_mapq is None and args.min_aln_length is None:
        parser.error("At least one filtering parameter (--max_mismatch_prop, --min_mapq, --min_aln_length) must be set.")

    # Create a temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp(prefix="temp_bam_sort_")
    namesorted_bam = os.path.join(temp_dir, "namesorted.bam")
    filtered_unsorted_bam = os.path.join(temp_dir, "filtered_unsorted.bam")

    # Step 1: Sort input BAM by read name using samtools sort -n
    sort_cmd = ["samtools", "sort", "-n", "-@", str(args.threads), "-o", namesorted_bam, args.bam]
    subprocess.run(sort_cmd, check=True)

    # Step 2: Filter the namesorted BAM
    bam_in = pysam.AlignmentFile(namesorted_bam, "rb")
    out_bam = pysam.AlignmentFile(filtered_unsorted_bam, "wb", header=bam_in.header)

    for read_name, group in itertools.groupby(bam_in.fetch(until_eof=True), key=lambda aln: aln.query_name):
        group_list = list(group)
        passing = [aln for aln in group_list if passes_filters(aln, args.max_mismatch_prop, args.min_mapq, args.min_aln_length)]
        if not passing:
            continue
        updated = update_read_group(passing)
        for aln in updated:
            out_bam.write(aln)

    bam_in.close()
    out_bam.close()

    # Step 3: Coordinate sort the filtered BAM using samtools sort
    sort_final_cmd = ["samtools", "sort", "-@", str(args.threads), "-o", args.out, filtered_unsorted_bam]
    subprocess.run(sort_final_cmd, check=True)

    # Step 4: Create BAM index using samtools index
    index_cmd = ["samtools", "index", args.out]
    subprocess.run(index_cmd, check=True)

    # Cleanup: Delete intermediate files and temporary directory
    os.remove(namesorted_bam)
    os.remove(filtered_unsorted_bam)
    os.rmdir(temp_dir)

if __name__ == "__main__":
    main()

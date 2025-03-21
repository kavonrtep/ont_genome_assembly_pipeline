#!/usr/bin/env python3
import argparse
from Bio import SeqIO


def read_removal_ids(removal_files):
    """Read one or more files containing read IDs (first column) to remove."""
    removal_set = set()
    for filename in removal_files:
        with open(filename, "r") as fh:
            for line in fh:
                parts = line.strip().split()
                if parts:
                    removal_set.add(parts[0])
    return removal_set


def extract_consensus_parts(input_file, file_format):
    """
    First pass: Scan the input FASTA/Q file and extract a set of individual read IDs
    that appear as parts of consensus reads (i.e. IDs that come from splitting a record ID
    that contains a semicolon).
    """
    consensus_parts = set()
    with open(input_file, "r") as fh:
        for record in SeqIO.parse(fh, file_format):
            if ";" in record.id:
                # Split on semicolon and add each trimmed part.
                parts = record.id.split(";")
                for part in parts:
                    consensus_parts.add(part.strip())
    return consensus_parts


def filter_records(
        input_file, file_format, removal_set, consensus_parts, output_file, trim_start,
        trim_end
        ):
    """
    Second pass: Filter records from input_file.
      - Skip any record whose ID is in removal_set.
      - If a record's ID contains a semicolon (i.e. a consensus record), output it with its
        new ID (first part + "_d") after trimming.
      - Otherwise, if a record's ID is in consensus_parts (i.e. it is an individual read that is
        part of a consensus), skip it.
      - Else, output the record as-is after trimming.

    Trimming: For each record that is output, remove the first 'trim_start' bases and the last
    'trim_end' bases. If the remaining sequence would be empty or negative in length,
    the record is skipped.

    For FASTQ records, quality scores are trimmed accordingly.
    """
    with open(input_file, "r") as infh, open(output_file, "w") as outfh:
        for record in SeqIO.parse(infh, file_format):
            # If record.id is in the external removal list, skip it.
            if record.id in removal_set:
                continue

            # For consensus records (IDs containing ";"), rename and process.
            if ";" in record.id:
                parts = record.id.split(";")
                new_id = parts[0].strip() + "_d"
                record.id = new_id
                record.name = new_id
                record.description = new_id
            else:
                # If the record is an individual read that appears as part of a consensus, skip it.
                if record.id in consensus_parts:
                    continue

            # Check if the record is long enough for trimming.
            if len(record.seq) <= (trim_start + trim_end):
                continue  # Skip records that would become empty after trimming.

            orig_length = len(record.seq)
            # For FASTQ, capture the trimmed quality scores first.
            if file_format == "fastq" and "phred_quality" in record.letter_annotations:
                new_qual = record.letter_annotations["phred_quality"][
                           trim_start: orig_length - trim_end]
            else:
                new_qual = None

            # Compute the trimmed sequence.
            trimmed_seq = record.seq[trim_start: orig_length - trim_end]

            # Clear letter annotations to allow reassigning record.seq.
            record.letter_annotations = {}
            record.seq = trimmed_seq

            # Reassign trimmed quality scores if applicable.
            if new_qual is not None:
                record.letter_annotations["phred_quality"] = new_qual

            # Write the filtered (and trimmed) record.
            SeqIO.write(record, outfh, file_format)


def main():
    parser = argparse.ArgumentParser(
        description="Filter a FASTA/FASTQ file by removing records whose IDs are in removal files or that represent individual reads of consensus triplets. "
                    "For consensus reads (IDs like 'AA;BB'), output only the consensus record (renamed to 'AA_d'). "
                    "Additionally, trim a specified number of bases from the beginning and end of each record."
        )
    parser.add_argument(
        "-i", "--input", required=True, help="Input FASTA or FASTQ file"
        )
    parser.add_argument(
        "-r", "--removal_files", nargs="+", required=True,
        help="File(s) with list of read IDs to remove (ID in first column)"
        )
    parser.add_argument(
        "-o", "--output", required=True, help="Output file (in same format as input)"
        )
    parser.add_argument(
        "--format", choices=["fasta", "fastq"], default="fasta",
        help="Input file format (default: fasta)"
        )
    parser.add_argument(
        "--trim_start", type=int, default=100,
        help="Number of bases to remove from the beginning of each read (default: 100)"
        )
    parser.add_argument(
        "--trim_end", type=int, default=10,
        help="Number of bases to remove from the end of each read (default: 10)"
        )
    args = parser.parse_args()

    print("Extracting consensus read components from input file...")
    consensus_parts = extract_consensus_parts(args.input, args.format)
    print(f"Found {len(consensus_parts)} consensus component IDs.")

    removal_set = read_removal_ids(args.removal_files)
    print(f"Loaded {len(removal_set)} removal IDs from removal file(s).")

    print("Filtering and trimming records (second pass)...")
    filter_records(
        args.input, args.format, removal_set, consensus_parts, args.output,
        args.trim_start, args.trim_end
        )
    print("Filtering complete. Filtered records written to", args.output)


if __name__ == "__main__":
    main()


def filter_fasta_fastq(
        input_file, removal_file, output_file, file_format, trim_start=100, trim_end=10
        ):
    """
    Filter a FASTA or FASTQ file by removing records whose IDs are in a removal file,
    applying consensus filtering and trimming each record.
    """
    consensus_parts = extract_consensus_parts(input_file, file_format)
    removal_set = read_removal_ids([removal_file])
    filter_records(
        input_file, file_format, removal_set, consensus_parts, output_file, trim_start,
        trim_end
        )

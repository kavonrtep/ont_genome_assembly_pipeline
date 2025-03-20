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

def filter_records(input_file, file_format, removal_set, consensus_parts, output_file):
    """
    Second pass: Filter records from input_file.
      - Skip any record whose ID is in removal_set.
      - If a record's ID contains a semicolon (i.e. a consensus record), output it with its
        new ID (first part + "_d").
      - Otherwise, if a record's ID is in consensus_parts (i.e. it is an individual read that is
        part of a consensus), skip it.
      - Else, output the record as-is.
    """
    with open(input_file, "r") as infh, open(output_file, "w") as outfh:
        # Use a writer that writes one record at a time.
        for record in SeqIO.parse(infh, file_format):
            # If record.id is in the external removal list, skip it.
            if record.id in removal_set:
                continue
            # If this record is a consensus record (its id contains ";")
            if ";" in record.id:
                parts = record.id.split(";")
                new_id = parts[0].strip() + "_d"
                record.id = new_id
                record.name = new_id
                record.description = new_id
                SeqIO.write(record, outfh, file_format)
            else:
                # If the record is an individual read that appears as part of a consensus, skip it.
                if record.id in consensus_parts:
                    continue
                # Otherwise, write the record as is.
                SeqIO.write(record, outfh, file_format)

def main():
    parser = argparse.ArgumentParser(
        description="Filter a FASTA/FASTQ file by removing records whose IDs are in removal files or that represent individual reads of consensus triplets. "
                    "For consensus reads (IDs like 'AA;BB'), output only the consensus record (renamed to 'AA_d')."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA or FASTQ file")
    parser.add_argument("-r", "--removal_files", nargs="+", required=True,
                        help="File(s) with list of read IDs to remove (ID in first column)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file (in same format as input)")
    parser.add_argument("--format", choices=["fasta", "fastq"], default="fasta",
                        help="Input file format (default: fasta)")
    args = parser.parse_args()

    print("Extracting consensus read components from input file...")
    consensus_parts = extract_consensus_parts(args.input, args.format)
    print(f"Found {len(consensus_parts)} consensus component IDs.")

    removal_set = read_removal_ids(args.removal_files)
    print(f"Loaded {len(removal_set)} removal IDs from removal file(s).")

    print("Filtering records (second pass)...")
    filter_records(args.input, args.format, removal_set, consensus_parts, args.output)
    print("Filtering complete. Filtered records written to", args.output)

if __name__ == "__main__":
    main()



def filter_fasta_fastq(input_file, removal_file, output_file, file_format):
    """
    Filter a FASTA or FASTQ file by removing records whose IDs are in a removal file.
    """
    consensus_parts = extract_consensus_parts(input_file, file_format)
    removal_set = read_removal_ids([removal_file])
    filter_records(input_file, file_format, removal_set, consensus_parts, output_file)


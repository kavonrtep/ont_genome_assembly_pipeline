#!/usr/bin/env python
"""
Script to read a table of contig connection pairs and run extract_contig_ends2.py
for each unique pair.

The input table should have at least 4 columns:
  read_id, contig1, contig2, connection_type, ...
where connection_type is in the format "X->Y" (e.g. "5->5", "5->3", "3->5", "3->3").

The script determines the appropriate contig order and orientation for
extract_contig_ends2.py such that the final combined FASTA always has:
  - contig1: extracted in 3' orientation
  - contig2: extracted in 5' orientation

For each unique pair, an output directory is created with the name:
  contig1_<orient1>_contig2_<orient2>
and the extract_contig_ends2.py script is called with the parameters.

Usage:
  python run_extract_for_contig_pairs.py --pairs_file connections.tsv --window_size 100000 --bam my.bam --reads_fasta reads.fasta --out_dir base_output --genome_fasta genome.fasta --min_aln_length 50
"""

import argparse
import csv
import os
import subprocess
import concurrent.futures


def process_pair(pair, args):
    # Unpack the tuple for the pair.
    contigA, contigB, connection_type = pair
    # Determine the orientation and choose contigs.
    chosen_contig1, orient1, chosen_contig2, orient2 = choose_orientation(contigA, contigB, connection_type)

    # Create output directory name in the form: contig1_<orient1>_contig2_<orient2>
    pair_out_dir = os.path.join(args.out_dir, f"{chosen_contig1}_{orient1}_{chosen_contig2}_{orient2}")
    os.makedirs(pair_out_dir, exist_ok=True)

    # Build the command to call extract_contig_ends2.py
    cmd = [
        "python", args.extract_script,
        "--contig1", chosen_contig1,
        "--orient1", orient1,
        "--contig2", chosen_contig2,
        "--orient2", orient2,
        "--window_size", str(args.window_size),
        "--bam", args.bam,
        "--reads_fasta", args.reads_fasta,
        "--out_dir", pair_out_dir,
        "--genome_fasta", args.genome_fasta,
        "--min_aln_length", str(args.min_aln_length)
    ]

    print("Running command:", " ".join(cmd))

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        # Log errors to a file.
        with open(os.path.join(pair_out_dir, "error.txt"), "a") as f:
            f.write(" ".join(cmd) + "\n")
            f.write("Error running the command\n")
    return

def choose_orientation(contigA, contigB, connection_type):
    """
    Given a pair (contigA, contigB) and a connection type "X->Y" (X and Y are "3" or "5"),
    decide on the order and the orientations so that the final call to extract_contig_ends2.py
    yields:
      - contig1 that provides the 3' end
      - contig2 that provides the 5' end.

    The extraction script works as follows:
      - For contig1: F means extract original 3' region, RC means extract original 5' region (and then RC).
      - For contig2: F means extract original 5' region, RC means extract original 3' region (and then RC).

    Option 1: Use order as given:
      contig1 = A, contig2 = B.
      Then:
        orient1 = F if X=="3", else RC if X=="5"
        orient2 = F if Y=="5", else RC if Y=="3"

    Option 2: Swap the order:
      contig1 = B, contig2 = A.
      Then:
        orient1 = F if Y=="3", else RC if Y=="5"
        orient2 = F if X=="5", else RC if X=="3"

    We then choose the option with fewer RC flags (i.e. preferring F). If equal, choose Option 1.
    """
    X, _, Y = connection_type.partition("->")
    # Option 1 (no swap)
    orient1_opt1 = "F" if X == "3" else "RC"
    orient2_opt1 = "F" if Y == "5" else "RC"
    count1 = (1 if orient1_opt1 == "RC" else 0) + (1 if orient2_opt1 == "RC" else 0)

    # Option 2 (swap order)
    orient1_opt2 = "F" if Y == "3" else "RC"
    orient2_opt2 = "F" if X == "5" else "RC"
    count2 = (1 if orient1_opt2 == "RC" else 0) + (1 if orient2_opt2 == "RC" else 0)

    if count2 < count1:
        # Use swapped order
        return contigB, orient1_opt2, contigA, orient2_opt2
    else:
        return contigA, orient1_opt1, contigB, orient2_opt1


def main():
    parser = argparse.ArgumentParser(
        description="Run extract_contig_ends2.py for each unique contig pair from a connections table."
        )
    parser.add_argument(
        "--pairs_file", required=True,
        help="Input table file with columns: read_id, contig1, contig2, connection_type, ..."
        )
    parser.add_argument(
        "--window_size", required=False, type=int, default=100000,
        help="Window size (number of bases) to extract"
        )
    parser.add_argument("--bam", required=True, help="Input BAM file (indexed)")
    parser.add_argument(
        "--reads_fasta", required=True, help="FASTA file with full read sequences"
        )
    parser.add_argument("--out_dir", required=True, help="Base output directory")
    parser.add_argument(
        "--genome_fasta", required=True,
        help="Genome FASTA file (indexed for seqkit faidx)"
        )
    parser.add_argument(
        "--min_aln_length", type=int, default=10000, required=False,
        help="Minimum alignment length required for a read to be considered"
        )
    parser.add_argument(
        "--extract_script", required=True,
        help="Path to the extract_contig_ends2.py script"
        )
    parser.add_argument(
        "--delimiter", default="\t",
        help="Delimiter used in the pairs file (default: tab)",
        required=False
        )
    parser.add_argument(
        "--cpu", type=int, default=1,
        help="Number of CPU threads to use (default: 1)"
        )

    args = parser.parse_args()

    # Build a set of unique (contig1, contig2, connection_type) tuples.
    unique_pairs = {}
    with open(args.pairs_file, "r") as infile:
        reader = csv.reader(infile, delimiter=args.delimiter)
        # Skip header
        next(reader)
        for row in reader:
            if len(row) < 4:
                continue
            # row: read_id, contig1, contig2, connection_type, ...
            contigA = row[1].strip()
            contigB = row[2].strip()
            connection_type = row[3].strip()
            key = (contigA, contigB, connection_type)
            unique_pairs[key] = True

    print(f"Found {len(unique_pairs)} unique contig pairs and connection types.")


    # Use concurrent.futures to run multiple processes in parallel.
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpu) as executor:
        futures = {executor.submit(process_pair, pair, args) for pair in unique_pairs}
        for future in concurrent.futures.as_completed(futures):
            future.result()




if __name__ == "__main__":
    main()

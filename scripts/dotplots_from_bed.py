#!/usr/bin/env python
"""
Script to extract genomic regions for BED intervals and generate self dotplots.
For each interval, it extracts the region defined by the BED interval ±150k bases,
skips intervals too close to the ends, and creates a self dotplot using Gepard.
Each interval gets its own directory named "chrom_bedStart_bedEnd" with the extracted
sequence (genomic_sequence.fasta) and a "dotplots" subdirectory containing the dotplot image.

External tools required: seqkit, Gepard.
"""

import argparse
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import concurrent.futures


def get_contig_length(genome_fasta, contig):
    """Return the length of a contig from the genome FASTA using Biopython."""
    with open(genome_fasta, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            if record.id == contig:
                return len(record.seq)
    raise ValueError(f"Contig {contig} not found in {genome_fasta}")


def extract_genomic_region(genome_fasta, contig, start, end, out_file):
    """
    Extract a genomic region using seqkit faidx.
    Region format: contig:start-end.
    Writes the FASTA output to out_file and returns (sequence, region_string).
    """
    region = f"{contig}:{start}-{end}"
    cmd = ["seqkit", "faidx", genome_fasta, region]
    with open(out_file, "w") as fh:
        subprocess.check_call(cmd, stdout=fh)
    # Concatenate all sequence lines (skip header)
    seq = ""
    with open(out_file, "r") as fh:
        for line in fh:
            if not line.startswith(">"):
                seq += line.strip()
    return seq, region


def run_gepard(out_dir, seq1, seq2, word=None):
    """
    Run Gepard dotplot generation.

    seq1: path to first sequence file.
    seq2: path to second sequence file.
    The output image is stored in out_dir.
    """
    gepard_dir = os.getenv("GEPARD", "/mnt/raid/opt/gepard-1.30")
    seq1 = os.path.abspath(seq1)
    seq2 = os.path.abspath(seq2)
    out_dir = os.path.abspath(out_dir)
    seq1_file = os.path.basename(seq1)
    seq1_path = os.path.dirname(seq1)
    seq2_file = os.path.basename(seq2)
    seq2_path = os.path.dirname(seq2)
    outfile = os.path.join(out_dir, f"{seq1_file}_x_{seq2_file}.png")
    command = f"cd {gepard_dir}; ./gepardcmd.sh -seq1 {seq1_path}/{seq1_file} -seq2 {seq2_path}/{seq2_file} -matrix matrices/edna.mat"
    if word is not None:
        command += f" -word {word}"
    command += f" -outfile {outfile}; cd {out_dir}"
    ret = subprocess.call(command, shell=True)
    if ret != 0:
        raise RuntimeError(f"Error executing command: {command}")


def process_interval_self(bed_line, genome_fasta, extension, out_dir):
    """
    Process one BED interval:
      - Create a directory named "chrom_bedStart_bedEnd"
      - Skip the interval if it is too close to the 5' or 3' end
        (i.e., if there are not at least 'extension' bases available on either side)
      - Otherwise, extract the genomic sequence for (BED interval ± extension)
      - Generate a self dotplot (comparing the extracted sequence to itself)
        and save it in a 'dotplots' subdirectory.
    """
    fields = bed_line.strip().split()
    if len(fields) < 3:
        print("Skipping invalid BED line:", bed_line)
        return

    chrom = fields[0]
    bed_start = int(fields[1])
    bed_end = int(fields[2])

    contig_len = get_contig_length(genome_fasta, chrom)

    # Check if there are at least 'extension' bases available on either side.
    if bed_start < extension:
        print(
            f"Skipping interval {chrom}_{bed_start}_{bed_end}: too close to 5' end (requires at least {extension} bases upstream)."
            )
        return
    if bed_end + extension > contig_len:
        print(
            f"Skipping interval {chrom}_{bed_start}_{bed_end}: too close to 3' end (requires at least {extension} bases downstream)."
            )
        return

    # Convert BED 0-based start to 1-based extraction coordinates.
    seq_start = bed_start + 1 - extension
    seq_end = bed_end + extension

    interval_name = f"{chrom}_{bed_start}_{bed_end}"
    interval_dir = os.path.join(out_dir, interval_name)
    os.makedirs(interval_dir, exist_ok=True)
    print(
        f"Processing interval {interval_name} with region: {chrom}:{seq_start}-{seq_end}"
        )

    # Extract genomic sequence.
    genomic_seq_file = os.path.join(interval_dir, interval_name + "_genomic_sequence.fasta")
    genomic_seq, region_str = extract_genomic_region(
        genome_fasta, chrom, seq_start, seq_end, genomic_seq_file
        )
    print(f"Extracted genomic sequence for {region_str} into {genomic_seq_file}")

    # Create a local dotplots subdirectory.
    local_dotplots_dir = interval_dir
    os.makedirs(local_dotplots_dir, exist_ok=True)

    # Generate self dotplot: compare genomic_seq_file against itself.
    try:
        run_gepard(local_dotplots_dir, genomic_seq_file, genomic_seq_file)
        print(f"Dotplot generated for interval {interval_name} in {local_dotplots_dir}")
    except Exception as e:
        print(f"Error generating dotplot for interval {interval_name}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract genomic regions (BED interval ± extension) and generate self dotplots."
        )
    parser.add_argument(
        "--bed", required=True, help="BED file (0-based, at least 3 columns)"
        )
    parser.add_argument(
        "--genome_fasta", required=True,
        help="Genome FASTA file (indexed for seqkit faidx)"
        )
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument(
        "--extension", type=int, default=150000,
        help="Extension on each side of the BED interval (default 150k)"
        )
    # cpu
    parser.add_argument(
        "--cpu", type=int, default=1,
        help="Number of CPU cores to use (default 1)"
        )

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Read and filter BED intervals.
    with open(args.bed, "r") as bed_fh:
        bed_lines = [line for line in bed_fh if line.strip() and not line.startswith("#")]

    # Process each interval in parallel.
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpu) as executor:
        futures = [executor.submit(
            process_interval_self, line, args.genome_fasta, args.extension, args.out_dir
            ) for line in bed_lines]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print("Error processing an interval:", e)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Extract full-length reads overlapping a genomic interval using external tools.
This script takes as input:
  --bam:        a BAM file (indexed)
  --interval:   a genomic interval in the format "contig:start-end" (commas allowed)
  --reads_fasta: FASTA file with full read sequences
  --genome_fasta: FASTA file with the genomic sequence
  --out_dir:    output directory

The script uses:
  - samtools view to extract read IDs overlapping the interval (separately for forward and reverse strands)
  - seqkit grep to extract those reads from the reads FASTA
  - For reverse reads, seqkit seq --reverse --complement is used to reorient them.
  - seqkit split is then used to split the combined FASTA into individual files.
  - Additionally, the genomic sequence for the input interval is extracted (via seqkit subseq).

Output:
  - A text file (read_ids.txt) listing the unique read IDs extracted.
  - A subdirectory (reads_split) containing individual FASTA files for each read.
  - A FASTA file (genomic_interval.fasta) with the genomic sequence for the interval.

Temporary files are used internally.
"""

import argparse
import os
import re
import subprocess
import tempfile


def parse_interval(interval_str):
    """Parse an interval string "contig:start-end" (commas allowed) and return (contig, start, end)."""
    interval_str = interval_str.replace(",", "")
    try:
        contig, pos = interval_str.split(":")
        start_str, end_str = pos.split("-")
        return contig.strip(), int(start_str), int(end_str)
    except Exception as e:
        raise ValueError(
            f"Interval '{interval_str}' is not in the expected format 'contig:start-end'"
            ) from e


def get_read_ids(bam, region, flag_filter):
    """
    Use samtools view with flag_filter (e.g. "-F 16" for forward or "-f 16" for reverse)
    to get read IDs overlapping the region.
    Returns a set of read IDs.
    """
    cmd = ["samtools", "view", flag_filter, bam, region]
    # Execute command; each line: read ID is the first field.
    output = subprocess.check_output(cmd, universal_newlines=True)
    read_ids = set()
    for line in output.strip().splitlines():
        fields = line.split("\t")
        if fields:
            read_ids.add(fields[0])
    return read_ids


def write_ids_to_file(ids, filepath):
    """Write a set of IDs to a file, one per line."""
    with open(filepath, "w") as fh:
        for r in ids:
            fh.write(r + "\n")


def run_seqkit_grep(fasta, id_file, out_file):
    """Run seqkit grep to extract sequences matching IDs from id_file."""
    cmd = ["seqkit", "grep", "-f", id_file, fasta]
    with open(out_file, "w") as outf:
        subprocess.check_call(cmd, stdout=outf)


def run_seqkit_faidx(fasta, id_list, out_file):
    """Run seqkit faidx to extract sequences matching IDs from id_file."""
    with open(out_file, "w") as outf:
        for ids in id_list:
            cmd = ["seqkit", "faidx", fasta, ids]
            subprocess.check_call(cmd, stdout=outf)


def run_seqkit_revcomp(in_file, out_file):
    """Run seqkit seq --reverse --complement on in_file, write to out_file."""
    cmd = ["seqkit", "seq", "--reverse", "--complement", in_file]
    tmp_out = out_file + ".tmp"
    with open(tmp_out, "w") as outf:
        subprocess.check_call(cmd, stdout=outf)
    # add _rc suffix to all sequence names
    with open(tmp_out) as inf, open(out_file, "w") as outf:
        for line in inf:
            if line.startswith(">"):
                line = line.strip() + "_rc\n"
            outf.write(line)


def run_seqkit_cat(files, out_file):
    """Concatenate multiple FASTA files into one output FASTA."""
    with open(out_file, "w") as outf:
        for f in files:
            with open(f) as inf:
                outf.write(inf.read())


def run_seqkit_split(in_file, out_dir):
    """Run seqkit split to split FASTA file into individual files in out_dir."""
    cmd = ["seqkit", "split", "-i", "-O", out_dir, in_file]
    subprocess.check_call(cmd)
    # rename split files


def extract_genomic_sequence(genome_fasta, interval_str, out_file):
    """Extract genomic sequence for interval using seqkit subseq."""
    cmd = ["seqkit", "faidx", genome_fasta, sanitize_region(interval_str)]
    with open(out_file, "w") as outf:
        subprocess.check_call(cmd, stdout=outf)

def sanitize_region(region):
    # remove commas, spaces
    region = re.sub(r"[,\s]", "", region)
    return region

def main():
    parser = argparse.ArgumentParser(
        description="Extract full-length reads overlapping a genomic interval and the genomic sequence using external seqkit commands."
        )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (indexed)")
    parser.add_argument(
        "--interval", required=True,
        help='Genomic interval in the format "contig:start-end" (e.g. ptg000010l:168,961,035-169,022,156)'
        )
    parser.add_argument(
        "-f", "--reads_fasta", required=True, help="FASTA file with full read sequences"
        )
    parser.add_argument(
        "-g","--genome_fasta", required=True, help="FASTA file with the genomic sequence"
        )
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory")
    args = parser.parse_args()

    # Parse interval.
    contig, start, end = parse_interval(args.interval)
    region_str = f"{contig}:{start}-{end}"

    os.makedirs(args.out_dir, exist_ok=True)

    # Use samtools view to get read IDs for forward and reverse separately.
    forward_ids = get_read_ids(args.bam, region_str, "-F 16")
    print(f"Found {len(forward_ids)} forward read IDs.")
    reverse_ids = get_read_ids(args.bam, region_str, "-f 16")
    print(f"Found {len(reverse_ids)} reverse read IDs.")

    all_ids = forward_ids.union(reverse_ids)
    # Write all unique read IDs to a text file.
    ids_txt = os.path.join(args.out_dir, "read_ids.txt")
    write_ids_to_file(all_ids, ids_txt)
    print(f"Unique read IDs written to {ids_txt}")

    # Create temporary directory for intermediate FASTA files.
    with tempfile.TemporaryDirectory() as tmpdir:
        forward_fasta = os.path.join(tmpdir, "forward_reads.fasta")
        reverse_fasta = os.path.join(tmpdir, "reverse_reads.fasta")
        reverse_rc_fasta = os.path.join(tmpdir, "reverse_reads_rc.fasta")
        combined_fasta = os.path.join(tmpdir, "reads.fasta")

        # Use seqkit grep to extract forward and reverse reads.
        run_seqkit_faidx(args.reads_fasta, forward_ids, forward_fasta)
        run_seqkit_faidx(args.reads_fasta, reverse_ids, reverse_fasta)

        # For reverse reads, get reverse complement.
        run_seqkit_revcomp(reverse_fasta, reverse_rc_fasta)

        # Concatenate forward and reverse (reverse-complemented) reads.
        run_seqkit_cat([forward_fasta, reverse_rc_fasta], combined_fasta)

        # Use seqkit split to split combined FASTA into individual FASTA files.
        reads_split_dir = os.path.join(args.out_dir, "reads_split")
        os.makedirs(reads_split_dir, exist_ok=True)
        run_seqkit_split(combined_fasta, reads_split_dir)

        # Also, copy the combined FASTA to the output directory if desired.
        combined_out = os.path.join(args.out_dir, "reads.fasta")
        subprocess.check_call(["cp", combined_fasta, combined_out])


    # Extract genomic sequence for the interval using seqkit subseq.
    genomic_out = os.path.join(
        args.out_dir,
        F"{sanitize_region(args.interval)}.fasta")
    print("Extracting genomic sequence for the interval using seqkit subseq...")
    extract_genomic_sequence(args.genome_fasta, args.interval, genomic_out)

    # run custom script on extracted fasta files:
    # script /mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/scripts/run_gepard_any_two_seq.pl
    # Params are
    # -d output_dir
    # -a seq1
    # -b seq2 (or seq1 for self-comparison)
    # create directory - dotplots
    # make read self-comparison dotplot
    # make read-genome aganst all reads dotplots (individually)

    def run_gepard_petr(out_dir, seq1, seq2):
        cmd = [
            "/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/scripts/run_gepard_any_two_seq.pl",
            "-d", out_dir,
            "-a", seq1,
            "-b", seq2
        ]
        subprocess.check_call(cmd)

    def run_gepard_single(out_dir, seq1, seq2, word=None):
        print("seq1:", seq1)
        print("seq2:", seq2)
        """
        Runs gepard using the provided sequence files and output directory.

        Parameters:
          seq1   : Path to the first sequence file.
          seq2   : Path to the second sequence file.
          out_dir: Output directory for the gepard result.
          word   : (Optional) Value for the '-word' option (default is None).

        The function constructs and executes the command:
            cd /mnt/raid/opt/gepard-1.30;
            ./gepardcmd.sh -seq1 <abs_seq1_dir>/<seq1_file>
                            -seq2 <abs_seq2_dir>/<seq2_file>
                            -matrix matrices/edna.mat
                            [ -word <word> ]
                            -outfile <out_dir>/<seq1_file>_x_<seq2_file>.png;
            cd <out_dir>
        If the command fails, a RuntimeError is raised.
        """
        # set gepard directory to environment variable GEPARD if exists
        gepard_dir = os.getenv("GEPARD")
        if not gepard_dir:
            gepard_dir = "/mnt/raid/opt/gepard-1.30"

        # Convert input paths to absolute paths
        seq1 = os.path.abspath(seq1)
        seq2 = os.path.abspath(seq2)
        out_dir = os.path.abspath(out_dir)

        # Extract file names and their directories
        seq1_file = os.path.basename(seq1)
        seq1_path = os.path.dirname(seq1)

        seq2_file = os.path.basename(seq2)
        seq2_path = os.path.dirname(seq2)

        # Build the output filename (png)
        outfile = os.path.join(out_dir, f"{seq1_file}_x_{seq2_file}.png")

        # Build the command string step by step.
        # Change directory to the gepard installation directory.
        command = f"cd {gepard_dir}; ./gepardcmd.sh -seq1 {seq1_path}/{seq1_file} -seq2 {seq2_path}/{seq2_file} -matrix matrices/edna.mat"

        # Add the -word option if provided.
        if word is not None:
            command += f" -word {word}"

        # Add outfile parameter and then change directory to out_dir.
        command += f" -outfile {outfile}; cd {out_dir}"

        # Execute the command.
        ret = subprocess.call(command, shell=True)
        if ret != 0:
            raise RuntimeError(f"Error when executing command: {command}")

    # run gepard on individual read sequences
    def run_gepard_on_read_dir(reads_split_dir, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for root, dirs, files in os.walk(reads_split_dir):
            for file in files:
                if file.endswith(".fasta"):
                    seq1 = os.path.join(root, file)
                    run_gepard_single(outdir, seq1, seq1)
    def run_gepard_on_genomic_read(reads_split_dir, genomic_out, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for root, dirs, files in os.walk(reads_split_dir):
            for file in files:
                if file.endswith(".fasta"):
                    seq1 = os.path.join(root, file)
                    run_gepard_single(outdir, genomic_out, seq1)

    run_gepard_on_read_dir(reads_split_dir,
                           os.path.join(args.out_dir,"dotplots_reads"))
    run_gepard_on_genomic_read(reads_split_dir, genomic_out,
                               os.path.join(args.out_dir, "dotplots_genomic_reads"))


if __name__ == "__main__":
    main()

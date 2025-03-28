#!/usr/bin/env python
"""
Modified script to extract reads overlapping a specified window from two contigs,
adjust orientations based on user-specified orientation (F or RC), and combine the two
contig parts (with a gap) into a FASTA file for dotplotting.

Inputs:
  --contig1      : Name of contig1
  --orient1      : Orientation of contig1 (F for forward, RC for reverse complement)
  --contig2      : Name of contig2
  --orient2      : Orientation of contig2 (F for forward, RC for reverse complement)
  --window_size  : Number of bases to extract from the contig end.
                  Extraction is done using original coordinates as follows:
                    - For contig1: if orient1==F, extract the 3' end;
                                 if orient1==RC, extract the 5' end.
                    - For contig2: if orient2==F, extract the 5' end;
                                 if orient2==RC, extract the 3' end.
                  After extraction, if orientation is RC, the sequence is reverse complemented.
  --bam          : BAM file (indexed)
  --reads_fasta  : FASTA file with full read sequences
  --out_dir      : Output directory
  --genome_fasta : Genome FASTA file (indexed for seqkit faidx)

Tools used:
  - samtools view for read extraction from BAM.
  - seqkit (faidx, grep, seq, split, cat) for FASTA processing.
  - Biopython for extracting contig lengths.
"""

import argparse
import os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
import pysam

def get_contig_length(genome_fasta, contig):
    """Return the length of contig from genome_fasta using Biopython."""
    with open(genome_fasta, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            if record.id == contig:
                return len(record.seq)
    raise ValueError(f"Contig {contig} not found in {genome_fasta}")

def extract_contig_region(genome_fasta, contig, which_end, window_size, out_file):
    """
    Extract a window from the contig using seqkit faidx.
    For the 5' end, extract bases 1-window_size.
    For the 3' end, extract bases (length-window_size+1)-length.
    Returns the extracted sequence (as a string) and the region string.
    """
    contig_len = get_contig_length(genome_fasta, contig)

    if which_end == "5":
        start = 1
        end = window_size
    elif which_end == "3":
        start = contig_len - window_size + 1
        end = contig_len
    else:
        raise ValueError("which_end must be '5' or '3'")
    # what is end and start are out of bounds? - adjust them
    if start < 1:
        start = 1
    if end > contig_len:
        end = contig_len
    region = f"{contig}:{start}-{end}"
    cmd = ["seqkit", "faidx", genome_fasta, region]
    # extract the region to a file
    print(cmd)
    with open(out_file, "w") as fh:
        subprocess.check_call(cmd, stdout=fh)
    # Read the extracted sequence (concatenate all non-header lines)
    seq = ""
    with open(out_file, "r") as fh:
        for line in fh:
            if not line.startswith(">"):
                seq += line.strip()
    return seq, region

def reverse_complement(seq):
    """Return the reverse complement of the sequence."""
    return str(Seq(seq).reverse_complement())

def get_read_ids(bam, region, flag_filter):
    """
    Use samtools view with the given flag filter to get read IDs overlapping the region.
    flag_filter: e.g. "-F 16" for forward, "-f 16" for reverse.
    Returns a set of read IDs.
    """
    cmd = ["samtools", "view", flag_filter, bam, region]
    output = subprocess.check_output(cmd, universal_newlines=True)
    read_ids = set()
    for line in output.strip().splitlines():
        fields = line.split("\t")
        if fields:
            read_ids.add(fields[0])
    return read_ids

def write_ids_to_file(ids, filepath):
    """Write a set of IDs to a file (one per line)."""
    with open(filepath, "w") as fh:
        for r in ids:
            fh.write(r + "\n")


def run_seqkit_faidx_ids(fasta, id_list, out_file):
    """Run seqkit faidx to extract sequences matching IDs from id_file."""
    with open(out_file, "w") as outf:
        for ids in id_list:
            cmd = ["seqkit", "faidx", fasta, ids]
            subprocess.check_call(cmd, stdout=outf)

def run_seqkit_revcomp(in_file, out_file):
    """
    Run seqkit seq --reverse --complement on in_file and write to out_file.
    Also appends a suffix (_rc) to the sequence names.
    """
    cmd = ["seqkit", "seq", "--reverse", "--complement", in_file]
    tmp_out = out_file + ".tmp"
    with open(tmp_out, "w") as outf:
        subprocess.check_call(cmd, stdout=outf)
    with open(tmp_out, "r") as inf, open(out_file, "w") as outf:
        for line in inf:
            if line.startswith(">"):
                line = line.strip() + "_rc\n"
            outf.write(line)
    os.remove(tmp_out)

def run_seqkit_cat(files, out_file):
    """Concatenate multiple FASTA files into one output FASTA."""
    with open(out_file, "w") as outf:
        for f in files:
            with open(f, "r") as inf:
                outf.write(inf.read())

def run_seqkit_split(in_file, out_dir):
    """Split a FASTA file into individual files using seqkit split."""
    cmd = ["seqkit", "split", "-i", "-O", out_dir, in_file]
    subprocess.check_call(cmd)

def extract_reads_for_region(
        bam, region, reads_fasta, out_dir, reorient=False, min_aln_length=10000
        ):
    """
    Extract reads overlapping a given region, filtering out alignments with a
    query alignment length below min_aln_length. This version uses pysam.

    Standard procedure:
      - Use pysam to iterate over alignments overlapping the region.
      - Separate read IDs into forward and reverse based on the alignment orientation.
      - Use seqkit faidx (not grep) to extract reads from reads_fasta.
      - For one of the sets, run reverse complement.
    When reorient is False (default):
      - Forward reads are kept as is; reverse reads are reverse complemented.
    When reorient is True:
      - Forward reads are reverse complemented; reverse reads are kept as is.

    Parameters:
      bam            : BAM file (indexed)
      region         : Region in format "contig:start-end"
      reads_fasta    : FASTA file with full read sequences
      out_dir        : Output directory to store extracted reads and temporary files
      reorient       : Boolean flag for read reorientation
      min_aln_length : Minimum alignment length required for a read to be considered.
    """
    os.makedirs(out_dir, exist_ok=True)
    print("BAM file:", bam)
    print("Region:", region)

    # Parse the region string "contig:start-end"
    chrom, pos = region.split(":")
    start_str, end_str = pos.split("-")
    start = int(start_str)
    end = int(end_str)

    fwd_ids = set()
    rev_ids = set()
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        # pysam uses 0-based start indexing.
        for aln in bamfile.fetch(chrom, start - 1, end):
            # Use query_alignment_length to filter reads.
            if aln.query_alignment_length is None or aln.query_alignment_length < min_aln_length:
                continue
            # skip MAPQ == 0 reads
            if aln.mapping_quality == 0:
                continue
            if not aln.is_reverse:
                fwd_ids.add(aln.query_name)
            else:
                rev_ids.add(aln.query_name)

    print(f"Found {len(fwd_ids)} forward reads")
    print(f"Found {len(rev_ids)} reverse reads")
    all_ids = fwd_ids.union(rev_ids)
    ids_txt = os.path.join(out_dir, "read_ids.txt")
    write_ids_to_file(all_ids, ids_txt)

    with tempfile.TemporaryDirectory() as tmpdir:
        fwd_fasta = os.path.join(tmpdir, "forward_reads.fasta")
        rev_fasta = os.path.join(tmpdir, "reverse_reads.fasta")
        rev_rc_fasta = os.path.join(tmpdir, "reverse_reads_rc.fasta")
        # Extract reads using seqkit faidx (iterating over each read ID)
        run_seqkit_faidx_ids(reads_fasta, fwd_ids, fwd_fasta)
        run_seqkit_faidx_ids(reads_fasta, rev_ids, rev_fasta)

        if not reorient:
            # Standard: reverse complement the reverse reads.
            run_seqkit_revcomp(rev_fasta, rev_rc_fasta)
            reads_to_cat = [fwd_fasta, rev_rc_fasta]
        else:
            # Reorient: reverse complement the forward reads.
            run_seqkit_revcomp(fwd_fasta, fwd_fasta + ".rc")
            reads_to_cat = [fwd_fasta + ".rc", rev_fasta]

        combined_fasta = os.path.join(tmpdir, "reads_combined.fasta")
        run_seqkit_cat(reads_to_cat, combined_fasta)
        reads_split_dir = os.path.join(out_dir, "reads_split")
        os.makedirs(reads_split_dir, exist_ok=True)
        run_seqkit_split(combined_fasta, reads_split_dir)
        final_reads_fasta = os.path.join(out_dir, "reads.fasta")
        subprocess.check_call(["cp", combined_fasta, final_reads_fasta])

    print(f"Reads for region {region} extracted to {out_dir}")


def run_gepard(out_dir, seq1, seq2,  word=None):
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

def run_gepard_perl(out_dir, seq1, seq2):
    cmd = [
        "/mnt/ceph/454_data/MinION/analysis/centromere_assembly_Cameor/scripts/run_gepard_any_two_seq.pl",
        "-d", out_dir,
        "-a", seq1,
        "-b", seq2
    ]
    subprocess.check_call(cmd)

def run_gepard_on_read_dir(reads_split_dir, outdir, word=None):
    """Run gepard on individual read sequences (self dotplot)."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for root, dirs, files in os.walk(reads_split_dir):
        for file in files:
            if file.endswith(".fasta"):
                seq1 = os.path.join(root, file)
                run_gepard(outdir, seq1, seq1, word)

def run_gepard_on_genomic_read(reads_split_dir, genomic_out, outdir, word=None):
    """Run gepard for contigs against reads dotplots."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for root, dirs, files in os.walk(reads_split_dir):
        for file in files:
            if file.endswith(".fasta"):
                seq1 = os.path.join(root, file)
                run_gepard(outdir, genomic_out, seq1, word)




def main():
    parser = argparse.ArgumentParser(
        description="Extract reads overlapping specified contig regions and combine "
                    "contig parts for dotplot."
    )
    parser.add_argument("--contig1", required=True, help="Name of contig1")
    parser.add_argument("--orient1", required=True, choices=["F", "RC"],
                        help="Orientation of contig1 (F for forward, RC for reverse complement)")
    parser.add_argument("--contig2", required=True, help="Name of contig2")
    parser.add_argument("--orient2", required=True, choices=["F", "RC"],
                        help="Orientation of contig2 (F for forward, RC for reverse complement)")
    parser.add_argument("--window_size", required=False, type=int,
                        help="Window size (number of bases) to extract", default=100000)
    parser.add_argument("--bam", required=True, help="Input BAM file (indexed)")
    parser.add_argument("--reads_fasta", required=True,
                        help="FASTA file with full read sequences")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--genome_fasta", required=True,
                        help="Genome FASTA file (indexed for seqkit faidx)")
    parser.add_argument("--min_aln_length", type=int, default=10000,
                        help="Minimum alignment length required for a read to be considered")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Determine which region (original coordinates) to extract for each contig based on orientation.
    # For contig1:
    #   If orient1 == F: we want the final contig part to be the 3' end, so extract original 3' region.
    #   If orient1 == RC: we want the final contig part (after RC) to be the 3' end, which comes from the original 5' region.
    contig1_region_type = "3" if args.orient1 == "F" else "5"
    # For contig2:
    #   If orient2 == F: we want the final contig part to be the 5' end, so extract original 5' region.
    #   If orient2 == RC: we want the final contig part (after RC) to be the 5' end, which comes from the original 3' region.
    contig2_region_type = "5" if args.orient2 == "F" else "3"

    # Extract regions from the genome (using original coordinates).
    contig1_out = os.path.join(args.out_dir, f"{args.contig1}_{contig1_region_type}.fasta")
    contig1_seq, contig1_region = extract_contig_region(
        args.genome_fasta, args.contig1, contig1_region_type, args.window_size, contig1_out
    )
    print(f"Extracted {contig1_region} to {contig1_out}")
    contig2_out = os.path.join(args.out_dir, f"{args.contig2}_{contig2_region_type}.fasta")
    contig2_seq, contig2_region = extract_contig_region(
        args.genome_fasta, args.contig2, contig2_region_type, args.window_size, contig2_out
    )
    print(f"Extracted {contig2_region} to {contig2_out}")
    # For final combined FASTA, adjust sequences so that:
    #   - contig1 part is in 3' orientation.
    #   - contig2 part is in 5' orientation.
    # If contig1 orientation is RC, reverse complement the extracted 5' region.
    if args.orient1 == "RC":
        contig1_seq = reverse_complement(contig1_seq)
    # If contig2 orientation is RC, reverse complement the extracted 3' region.
    if args.orient2 == "RC":
        contig2_seq = reverse_complement(contig2_seq)

    # Extract reads overlapping the specified regions.
    # For contig1: if orient1 is RC, pass reorient=True so that forward reads are reversed.
    contig1_reads_dir = os.path.join(args.out_dir, f"{args.contig1}_reads")
    extract_reads_for_region(
        args.bam, contig1_region, args.reads_fasta, contig1_reads_dir,
        reorient=(args.orient1 == "RC"),
        min_aln_length=args.min_aln_length
    )
    # For contig2: if orient2 is RC, pass reorient=True.
    contig2_reads_dir = os.path.join(args.out_dir, f"{args.contig2}_reads")
    extract_reads_for_region(
        args.bam, contig2_region, args.reads_fasta, contig2_reads_dir,
        reorient=(args.orient2 == "RC"),
        min_aln_length=args.min_aln_length
    )

    # Create a combined contig FASTA file (for dotplotting).
    combined_fasta = os.path.join(args.out_dir, "combined_contigs.fasta")
    gap = "N" * 5000
    with open(combined_fasta, "w") as fh:
        fh.write(f">{args.contig1}_part\n{contig1_seq}\n")
        fh.write(">gap\n" + gap + "\n")
        fh.write(f">{args.contig2}_part\n{contig2_seq}\n")
    print(f"Combined contig FASTA written to {combined_fasta}")

    # Run gepard dotplots.
    contig1_reads_dir_split = os.path.join(contig1_reads_dir, "reads_split")
    contig2_reads_dir_split = os.path.join(contig2_reads_dir, "reads_split")
    dotplot_reads_dir = os.path.join(args.out_dir, "dotplot_reads")
    dotplot_genomic_dir = os.path.join(args.out_dir, "dotplot_genomic")
    dotplot_genomic_dir_w150 = os.path.join(args.out_dir, "dotplot_genomic_w150")
    run_gepard_on_read_dir(contig1_reads_dir_split, dotplot_reads_dir)
    run_gepard_on_read_dir(contig2_reads_dir_split, dotplot_reads_dir)
    run_gepard_on_genomic_read(contig1_reads_dir_split, combined_fasta, dotplot_genomic_dir)
    run_gepard_on_genomic_read(contig2_reads_dir_split, combined_fasta, dotplot_genomic_dir)
    run_gepard_on_genomic_read(contig1_reads_dir_split, combined_fasta, dotplot_genomic_dir_w150, word=150)
    run_gepard_on_genomic_read(contig2_reads_dir_split, combined_fasta, dotplot_genomic_dir_w150, word=150)


if __name__ == "__main__":
    main()

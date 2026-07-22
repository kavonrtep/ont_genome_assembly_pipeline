#!/usr/bin/env python3
"""Merge per-haplotype `yak triobin` output into one parental/subgenome table.

`yak triobin pat.yak mat.yak contigs.fa` classifies each contig by parental
k-mer content. This helper collects those per-haplotype classifications and
maps the parental call onto user-supplied subgenome labels (e.g. AA vs B for an
AAB triploid), producing a single tidy TSV.

Note: the exact yak triobin column layout can vary between yak versions. The
contig-name and classification columns are configurable (--name_col /
--class_col); defaults are the first and last whitespace-delimited fields.
"""

import argparse
import sys


def read_fasta(path):
    """Minimal FASTA reader -> {name: sequence}. Uses the first whitespace
    token of the header as the name (matching PAF/BAM contig naming)."""
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_fasta(path, records):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def classify_token(token):
    """Normalise a yak triobin type token to paternal/maternal/ambiguous."""
    token = token.strip().lower()
    if token in ("p", "pat", "paternal", "1"):
        return "paternal"
    if token in ("m", "mat", "maternal", "2"):
        return "maternal"
    return "ambiguous"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-t", "--triobin", action="append", default=[], metavar="HAP:FILE",
                        help="Haplotype label and its yak triobin output, e.g. hap1:hap1.triobin. "
                             "Repeatable.")
    parser.add_argument("-o", "--output", required=True, help="Merged TSV output path")
    parser.add_argument("--paternal_label", default="paternal",
                        help="Subgenome label for paternal contigs (default: paternal)")
    parser.add_argument("--maternal_label", default="maternal",
                        help="Subgenome label for maternal contigs (default: maternal)")
    parser.add_argument("--name_col", type=int, default=0,
                        help="0-based column holding the contig name (default: 0)")
    parser.add_argument("--class_col", type=int, default=-1,
                        help="0-based column holding the class token (default: -1, last)")
    parser.add_argument("--fasta", action="append", default=[], metavar="HAP:FILE",
                        help="Haplotype label and its FASTA, e.g. hap1:hap1.fa. Required "
                             "with --bin_dir to write per-subgenome FASTAs. Repeatable.")
    parser.add_argument("--bin_dir", default=None,
                        help="If set, write {bin_dir}/{hap}.{subgenome}.fa binned by call.")
    args = parser.parse_args()

    label_map = {
        "paternal": args.paternal_label,
        "maternal": args.maternal_label,
        "ambiguous": "ambiguous",
    }

    rows = []
    for spec in args.triobin:
        if ":" not in spec:
            sys.exit(f"--triobin expects HAP:FILE, got: {spec}")
        hap, path = spec.split(":", 1)
        with open(path) as handle:
            for line in handle:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) <= max(args.name_col, args.class_col % len(fields)):
                    continue
                contig = fields[args.name_col]
                parent = classify_token(fields[args.class_col])
                rows.append((hap, contig, parent, label_map[parent]))

    with open(args.output, "w") as out:
        out.write("haplotype\tcontig\tparent\tsubgenome\n")
        for hap, contig, parent, subgenome in rows:
            out.write(f"{hap}\t{contig}\t{parent}\t{subgenome}\n")

    if args.bin_dir:
        import os

        os.makedirs(args.bin_dir, exist_ok=True)
        fasta_by_hap = {}
        for spec in args.fasta:
            hap, path = spec.split(":", 1)
            fasta_by_hap[hap] = read_fasta(path)
        # contig -> subgenome, per haplotype
        assignment = {}
        for hap, contig, _parent, subgenome in rows:
            assignment.setdefault(hap, {})[contig] = subgenome
        for hap, seqs in fasta_by_hap.items():
            bins = {}
            for contig, seq in seqs.items():
                subgenome = assignment.get(hap, {}).get(contig, "ambiguous")
                bins.setdefault(subgenome, []).append((contig, seq))
            for subgenome, records in bins.items():
                write_fasta(os.path.join(args.bin_dir, f"{hap}.{subgenome}.fa"), records)


if __name__ == "__main__":
    main()

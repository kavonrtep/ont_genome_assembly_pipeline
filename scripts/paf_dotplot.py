#!/usr/bin/env python3
"""Render a whole-genome dotplot from a minimap2 PAF (assembly vs reference).

The query sequences (assembly haplotype) are laid out along the X axis and the
target sequences (reference) along the Y axis in a single concatenated
coordinate system, with alignments drawn as line segments coloured by strand.
Intended for the dotplot_vs_reference step of Snakefile_generic.
"""

import argparse
import sys


def parse_paf(path, min_len):
    """Yield alignment records from a PAF file, filtering by block length."""
    records = []
    query_len = {}
    target_len = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            strand = fields[4]
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            block = int(fields[10])
            query_len[qname] = qlen
            target_len[tname] = tlen
            if block < min_len:
                continue
            records.append((qname, qstart, qend, strand, tname, tstart, tend))
    return records, query_len, target_len


def cumulative_offsets(lengths):
    """Return {name: offset} and the total, ordered by descending length."""
    offsets = {}
    running = 0
    for name, length in sorted(lengths.items(), key=lambda kv: -kv[1]):
        offsets[name] = running
        running += length
    return offsets, running


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--paf", required=True, help="Input PAF file")
    parser.add_argument("-o", "--output", required=True, help="Output PNG path")
    parser.add_argument("--min_len", type=int, default=2000,
                        help="Minimum alignment block length to plot (default: 2000)")
    parser.add_argument("--query_name", default="assembly", help="X axis label")
    parser.add_argument("--target_name", default="reference", help="Y axis label")
    parser.add_argument("--title", default=None, help="Plot title")
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    records, query_len, target_len = parse_paf(args.paf, args.min_len)
    if not records:
        # Still emit an (empty) figure so the rule produces its declared output.
        sys.stderr.write("paf_dotplot.py: no alignments passed the length filter\n")

    q_off, q_total = cumulative_offsets(query_len)
    t_off, t_total = cumulative_offsets(target_len)

    fig, ax = plt.subplots(figsize=(10, 10))
    for qname, qstart, qend, strand, tname, tstart, tend in records:
        x0 = q_off[qname] + qstart
        x1 = q_off[qname] + qend
        if strand == "+":
            y0 = t_off[tname] + tstart
            y1 = t_off[tname] + tend
            colour = "#1f77b4"
        else:
            y0 = t_off[tname] + tend
            y1 = t_off[tname] + tstart
            colour = "#d62728"
        ax.plot([x0, x1], [y0, y1], color=colour, linewidth=0.6)

    for offset in list(q_off.values())[1:]:
        ax.axvline(offset, color="0.85", linewidth=0.4, zorder=0)
    for offset in list(t_off.values())[1:]:
        ax.axhline(offset, color="0.85", linewidth=0.4, zorder=0)

    ax.set_xlim(0, max(q_total, 1))
    ax.set_ylim(0, max(t_total, 1))
    ax.set_xlabel(args.query_name)
    ax.set_ylabel(args.target_name)
    ax.set_title(args.title or f"{args.query_name} vs {args.target_name}")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

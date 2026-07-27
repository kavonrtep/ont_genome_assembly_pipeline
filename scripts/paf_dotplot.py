#!/usr/bin/env python3
"""Whole-genome dotplot from a minimap2 PAF (assembly haplotype vs reference).

Query (assembly) contigs are laid along X, ordered by the reference sequence and
position they best match, so co-linear assemblies read as a diagonal; contigs
with no reference hit are grouped at the right under "unplaced". Reference
sequences are laid along Y. Both axes carry sequence-name ticks (every reference
sequence; the largest query contigs), so the plot is legible rather than an
anonymous concatenation. Alignments are drawn coloured by strand.
"""

import argparse
import sys
from collections import defaultdict


def parse_paf(path, min_len):
    records = []
    query_len = {}
    target_len = {}
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qlen, qstart, qend = f[0], int(f[1]), int(f[2]), int(f[3])
            strand = f[4]
            tname, tlen, tstart, tend = f[5], int(f[6]), int(f[7]), int(f[8])
            block = int(f[10])
            query_len[qname] = qlen
            target_len[tname] = tlen
            if block < min_len:
                continue
            records.append((qname, qstart, qend, strand, tname, tstart, tend))
    return records, query_len, target_len


def query_placement(records, query_len):
    """For each query: (best target by aligned bases, weighted-mean target pos)."""
    bases = defaultdict(lambda: defaultdict(int))
    positions = defaultdict(lambda: defaultdict(list))
    for qname, qs, qe, _strand, tname, ts, te in records:
        w = max(qe - qs, 1)
        bases[qname][tname] += w
        positions[qname][tname].append(((ts + te) / 2.0, w))
    placement = {}
    for q in query_len:
        if q in bases:
            best = max(bases[q], key=bases[q].get)
            pw = positions[q][best]
            wsum = sum(w for _, w in pw) or 1
            placement[q] = (best, sum(p * w for p, w in pw) / wsum)
        else:
            placement[q] = (None, 0.0)
    return placement


def offsets(order, lengths):
    off = {}
    running = 0
    for name in order:
        off[name] = running
        running += lengths[name]
    return off, running


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--paf", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--min_len", type=int, default=2000,
                        help="Minimum alignment block length to plot (default: 2000)")
    parser.add_argument("--top_n_labels", type=int, default=30,
                        help="Label this many of the largest query contigs on X (default: 30)")
    parser.add_argument("--query_name", default="assembly")
    parser.add_argument("--target_name", default="reference")
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    records, query_len, target_len = parse_paf(args.paf, args.min_len)
    if not records:
        sys.stderr.write("paf_dotplot.py: no alignments passed the length filter\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Reference axis: sequences ordered by descending length.
    target_order = sorted(target_len, key=lambda t: -target_len[t])
    t_rank = {t: i for i, t in enumerate(target_order)}
    t_off, t_total = offsets(target_order, target_len)

    # Query axis: ordered by (best reference, position); unplaced grouped last.
    placement = query_placement(records, query_len)

    def qkey(q):
        t, pos = placement[q]
        if t is None:
            return (len(target_order), 0.0, -query_len[q])
        return (t_rank[t], pos, 0.0)

    query_order = sorted(query_len, key=qkey)
    q_off, q_total = offsets(query_order, query_len)
    first_unplaced = next((q_off[q] for q in query_order if placement[q][0] is None), None)

    fig, ax = plt.subplots(figsize=(11, 10))
    for qname, qs, qe, strand, tname, ts, te in records:
        x0, x1 = q_off[qname] + qs, q_off[qname] + qe
        if strand == "+":
            y0, y1, colour = t_off[tname] + ts, t_off[tname] + te, "#1f77b4"
        else:
            y0, y1, colour = t_off[tname] + te, t_off[tname] + ts, "#d62728"
        ax.plot([x0, x1], [y0, y1], color=colour, linewidth=0.6)

    for off in list(q_off.values())[1:]:
        ax.axvline(off, color="0.9", linewidth=0.3, zorder=0)
    for off in list(t_off.values())[1:]:
        ax.axhline(off, color="0.9", linewidth=0.3, zorder=0)
    if first_unplaced is not None:
        ax.axvline(first_unplaced, color="0.5", linewidth=0.8, linestyle="--", zorder=1)

    # Reference (Y) ticks: every sequence at its midpoint.
    ax.set_yticks([t_off[t] + target_len[t] / 2.0 for t in target_order])
    ax.set_yticklabels(target_order, fontsize=7)
    # Query (X) ticks: the largest contigs only, to stay legible.
    labelled = sorted(query_len, key=lambda q: -query_len[q])[:args.top_n_labels]
    labelled = [q for q in query_order if q in set(labelled)]
    ax.set_xticks([q_off[q] + query_len[q] / 2.0 for q in labelled])
    ax.set_xticklabels(labelled, fontsize=6, rotation=90)

    ax.set_xlim(0, max(q_total, 1))
    ax.set_ylim(0, max(t_total, 1))
    ax.set_xlabel(f"{args.query_name}  (contigs ordered by reference; dashed = unplaced)")
    ax.set_ylabel(args.target_name)
    ax.set_title(args.title or f"{args.query_name} vs {args.target_name}")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Karyoplot of terminal telomeric repeats from a tidk search table.

Reads tidk's per-window repeat counts and draws each telomere-bearing contig as
an ideogram bar on a single shared bp scale, with telomeric repeat arrays shown
as bands at the contig termini (contiguous dense windows collapsed into one
band, not per-window counts). Only contigs with a telomere close to a terminus
are shown. Bands are coloured by strand (forward vs reverse repeat), which
encodes orientation; when a RagTag AGP is supplied, each contig is first
reoriented to its reference orientation so telomeres display consistently.

tidk search TSV columns: id, window, forward_repeat_number, reverse_repeat_number, telomeric_repeat
"""

import argparse
import sys
from collections import defaultdict

FWD_COLOUR = "#1f77b4"
REV_COLOUR = "#d62728"


def read_fai(path):
    lengths = {}
    with open(path) as handle:
        for line in handle:
            f = line.split("\t")
            if len(f) >= 2:
                lengths[f[0]] = int(f[1])
    return lengths


def read_agp_orientation(path):
    """component_id -> orientation (+/-) from AGP 'W' rows."""
    orient = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[4] == "W":
                orient[f[5]] = f[8]
    return orient


def read_windows(path):
    """id -> sorted list of (pos_start, pos_end, forward, reverse)."""
    by_contig = defaultdict(list)
    positions = defaultdict(list)
    with open(path) as handle:
        header = handle.readline()
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            cid, window = f[0], int(f[1])
            fwd, rev = int(f[2]), int(f[3])
            positions[cid].append(window)
            by_contig[cid].append([window, fwd, rev])
    # infer window size per contig from the smallest positive step
    result = {}
    for cid, rows in by_contig.items():
        steps = [b - a for a, b in zip(sorted(positions[cid]), sorted(positions[cid])[1:]) if b > a]
        wsize = min(steps) if steps else 10000
        windows = []
        for window, fwd, rev in rows:
            windows.append((max(window - wsize, 0), window, fwd, rev))
        windows.sort()
        result[cid] = (windows, wsize)
    return result


def terminal_bands(windows, wsize, length, min_repeats, term_bp, term_frac):
    """Collapse dense terminal windows into bands: list of (start, end, fwd, rev)."""
    zone = max(term_bp, int(term_frac * length))
    dense = [w for w in windows if (w[2] + w[3]) >= min_repeats
             and (w[1] <= zone or w[0] >= length - zone)]
    bands = []
    for start, end, fwd, rev in dense:
        if bands and start - bands[-1][1] <= wsize:
            bands[-1][1] = end
            bands[-1][2] += fwd
            bands[-1][3] += rev
        else:
            bands.append([start, end, fwd, rev])
    return bands


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tsv", required=True, help="tidk search *_telomeric_repeat_windows.tsv")
    parser.add_argument("--fai", required=True, help="assembly .fai (contig lengths)")
    parser.add_argument("-o", "--output", required=True, help="output PNG")
    parser.add_argument("--agp", default=None, help="RagTag scaffold AGP for reference reorientation")
    parser.add_argument("--motif", default="", help="telomere motif (for the title)")
    parser.add_argument("--min_repeats", type=int, default=25,
                        help="min forward+reverse repeats per window to count as telomeric (default 25)")
    parser.add_argument("--terminal_bp", type=int, default=20000)
    parser.add_argument("--terminal_frac", type=float, default=0.05)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    lengths = read_fai(args.fai)
    windows_by_contig = read_windows(args.tsv)
    orient = read_agp_orientation(args.agp) if args.agp else {}

    # Build the set of telomere-bearing contigs (terminal telomere present).
    contigs = []
    for cid, (windows, wsize) in windows_by_contig.items():
        length = lengths.get(cid)
        if not length:
            continue
        bands = terminal_bands(windows, wsize, length, args.min_repeats,
                               args.terminal_bp, args.terminal_frac)
        if not bands:
            continue
        flip = orient.get(cid) == "-"
        if flip:
            bands = [[length - e, length - s, rev, fwd] for s, e, fwd, rev in bands]
        contigs.append((cid, length, bands, flip))

    contigs.sort(key=lambda c: -c[1])  # longest first

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    orient_note = "reoriented to reference (AGP)" if orient else "orientation as assembled"
    title = args.title or f"Terminal telomeres{' ' + args.motif if args.motif else ''} — {orient_note}"

    if not contigs:
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(0.5, 0.5, "no terminal telomeres detected", ha="center", va="center")
        ax.axis("off")
        ax.set_title(title)
        fig.savefig(args.output, dpi=150)
        sys.stderr.write("telomere_karyoplot.py: no terminal telomeres detected\n")
        return

    max_len = max(c[1] for c in contigs)
    fig, ax = plt.subplots(figsize=(11, max(2.0, 0.4 * len(contigs) + 1)))
    for row, (cid, length, bands, flip) in enumerate(contigs):
        y = len(contigs) - row
        ax.add_patch(plt.Rectangle((0, y - 0.3), length, 0.6, facecolor="0.9",
                                   edgecolor="0.6", linewidth=0.5))
        for s, e, fwd, rev in bands:
            colour = FWD_COLOUR if fwd >= rev else REV_COLOUR
            ax.add_patch(plt.Rectangle((s, y - 0.3), max(e - s, max_len * 0.002), 0.6,
                                       facecolor=colour, edgecolor="none"))

    ax.set_ylim(0.3, len(contigs) + 0.9)
    ax.set_xlim(-max_len * 0.01, max_len * 1.01)
    ax.set_yticks([len(contigs) - r for r in range(len(contigs))])
    ax.set_yticklabels([f"{c[0]}{' (-)' if c[3] else ''}" for c in contigs], fontsize=7)
    ax.set_xlabel("position (bp, shared scale)")
    ax.set_title(title)
    ax.legend(handles=[Patch(color=FWD_COLOUR, label="forward repeat"),
                       Patch(color=REV_COLOUR, label="reverse repeat")],
              loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()

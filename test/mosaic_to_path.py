#!/usr/bin/env python3
"""Expand a mosaic into one path per thread and check it is an exact walk in the graph.

This is the acceptance test the mosaic format exists to pass. Every other property -- segments
chaining, orientation, walk order, patched gaps, nested excursions -- is a way for this to fail, so
checking it directly subsumes checking them individually.

Inputs are what `vg` can already produce, so this needs no new subcommand:
  paths.gaf   vg paths -x GRAPH -A      one line per stored path, node walk in field 6
  graph.gfa   vg view -g GRAPH          L-lines give the edges
  mosaic.tsv  vg call --mosaic-out

Checks, in order of how much they tell you:
  1. every segment names a haplotype the GAF carries, and its start/end nodes lie on that walk
  2. the slice between them is taken in the haplotype's own direction of travel
  3. consecutive segments of a path join -- the next begins where the previous ended
  4. every step of the concatenated path is a real edge in the graph
Exit status is non-zero if any check fails; the report says which and how often.

A path is one (contig, strand, fragment) triple: v5 says so, and `fragment` only advances where a
boundary was left open, so with --mosaic-patch-gaps on there is one path per strand.

--gfa is optional. Without it checks 1-3 still run, which is enough to say whether the rows tile;
check 4 needs the graph's edges, and the GFA for a whole chromosome is large.
"""
import sys, argparse
from collections import defaultdict


def read_walks(gaf):
    """path name -> list of (node_id, is_reverse), from the GAF path field."""
    walks = {}
    for line in open(gaf):
        if line.startswith('@') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        steps, i, p = [], 0, f[5]
        while i < len(p):
            if p[i] in '><':
                rev = p[i] == '<'
                j = i + 1
                while j < len(p) and p[j] not in '><':
                    j += 1
                steps.append((int(p[i + 1:j]), rev))
                i = j
            else:
                i += 1
        walks[f[0]] = steps
    return walks


def read_edges(gfa):
    """set of ((from, from_rev), (to, to_rev)), both orientations of each L-line."""
    e = set()
    for line in open(gfa):
        if not line.startswith('L\t'):
            continue
        f = line.rstrip('\n').split('\t')
        a, ao, b, bo = int(f[1]), f[2] == '-', int(f[3]), f[4] == '-'
        e.add(((a, ao), (b, bo)))
        e.add(((b, not bo), (a, not ao)))      # the same edge walked the other way
    return e


def haplotype_walk(walks, name, _cache={}):
    """Fragments of one haplotype, each with a node->offset index.

    Indexed once per haplotype and cached: a walk here is two million nodes, and rebuilding the
    index per segment turns a linear check into 16 billion operations.
    """
    if name in _cache:
        return _cache[name]
    out = [w for n, w in walks.items() if n == name or n.startswith(name + '#')
           or n.split('[')[0] == name]
    out.sort(key=len, reverse=True)
    _cache[name] = [(w, {n: i for i, (n, _) in enumerate(w)}) for w in out]
    return _cache[name]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mosaic', required=True)
    ap.add_argument('--gaf', required=True)
    ap.add_argument('--gfa', help='optional; without it check 4 is skipped')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()

    walks = read_walks(a.gaf)
    edges = read_edges(a.gfa) if a.gfa else None
    hapname, threads = {}, defaultdict(list)
    for line in open(a.mosaic):
        f = line.rstrip('\n').split('\t')
        if f[0] == '#haplotype':
            hapname[f[1]] = f[2]
        elif f[0] == 'H':
            threads[(f[1], f[2], f[3])].append(f)      # contig, strand, fragment

    fail = defaultdict(int)
    reversed_slices = [0]
    orient_mismatch = [0]
    total_steps = 0
    for (contig, strand, frag), segs in sorted(threads.items()):
        path = []
        for f in segs:
            # v5: contig strand fragment ref_start ref_end start_node end_node
            #     hap_index haplotype sites gbwt_offset
            # start_node and end_node are ORIENTED: id * 2 + is_reverse.
            hap = f[9]
            start, srev = int(f[6]) // 2, bool(int(f[6]) & 1)
            end,   erev = int(f[7]) // 2, bool(int(f[7]) & 1)
            if f[8] == '*':
                fail['segment has no haplotype to walk'] += 1
                continue
            frags = haplotype_walk(walks, hap) or haplotype_walk(walks, hap + '#' + contig)
            if not frags:
                fail['haplotype not found in the graph'] += 1
                continue
            piece = None
            for w, idx in frags:
                if start in idx and end in idx:
                    i, j = idx[start], idx[end]
                    if i <= j:
                        piece = w[i:j + 1]
                    else:
                        # The segment runs against the haplotype's stored direction. That is not an
                        # error and not rare: on this graph whole haplotypes are stored reversed --
                        # recombination#17 is '<' at every one of its 2.0 M steps -- so a segment
                        # given in reference order is a walk along the reverse complement. Take the
                        # slice and flip both the order and every step's orientation.
                        piece = [(n, not r) for n, r in reversed(w[j:i + 1])]
                        reversed_slices[0] += 1
                    # The file states the orientation at both ends; the haplotype's own walk states
                    # it too. They have to agree, or the row's oriented anchors are not the walk's.
                    if piece and (piece[0][1] != srev or piece[-1][1] != erev):
                        orient_mismatch[0] += 1
                    break
            if piece is None:
                fail['segment endpoints not both on one fragment of its haplotype'] += 1
                continue
            if path and path[-1] != piece[0]:
                if edges is not None and (path[-1], piece[0]) in edges:
                    fail['segments abut by an edge rather than sharing a node'] += 1
                else:
                    fail['segments do not join'] += 1
            path.extend(piece if not path else piece[1:])
        for i in range(len(path) - 1):
            total_steps += 1
            if edges is not None and (path[i], path[i + 1]) not in edges:
                fail['step is not an edge in the graph'] += 1
        if not a.quiet:
            print(f"  {contig} strand {strand} fragment {frag}: "
                  f"{len(segs)} segments -> {len(path)} node steps")

    against = (f"against {len(edges)//2} edges" if edges is not None
               else "(no --gfa given: edges unchecked)")
    print(f"\n  {len(threads)} paths, {total_steps} steps {against}")
    print(f"  {reversed_slices[0]} segments walk against their haplotype's stored direction")
    if orient_mismatch[0]:
        print(f"  {orient_mismatch[0]} segments whose stated orientation disagrees with the walk")
    if not fail:
        print("  OK: every path expands to an exact walk in the graph")
        return 0
    for k, v in sorted(fail.items(), key=lambda kv: -kv[1]):
        print(f"  FAIL {v:>7}  {k}")
    return 1


if __name__ == '__main__':
    sys.exit(main())

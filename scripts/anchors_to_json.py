#!/usr/bin/env python3
"""Filter a `vg call --anchors-out` file and reshape it for an anchor-graph assembler.

Emits the shape Shasta's `FromJson` anchor creation reads:

    [ [ "anchorName", [ ["readName", strand, begin, end], ... ] ], ... ]

An anchor here is a **zero-length pin**, so `begin == end`. What that position means depends on the
coordinate frame, and the two available ones differ, so it is explicit rather than assumed:

  --frame oriented   (default) positions are in the read after reverse-complementing when strand is
                     1, which is the convention Shasta's own reader documents.
  --frame sequenced  positions are in the read as sequenced, which is what vg writes and what the
                     read file contains. Use this if the consumer indexes the raw reads.

`vg` writes `offset` as the index of the last base BEFORE the pin, reading in the site's direction,
in the read as sequenced. The pin itself therefore sits at `offset + 1` on strand 0 and at `offset`
on strand 1; converting to the oriented frame mirrors that through the read length, which is why
--frame oriented needs --reads.

Filtering happens here rather than in `vg` so that one run's output can be cut several ways.
"""

import argparse
import gzip
import json
import sys


def open_maybe_gzip(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def read_lengths(path):
    lengths = {}
    with open_maybe_gzip(path) as handle:
        first = handle.readline()
        handle.seek(0)
        if first.startswith("@"):
            while True:
                name = handle.readline()
                if not name:
                    break
                seq = handle.readline().rstrip("\n")
                handle.readline()
                handle.readline()
                lengths[name[1:].split()[0]] = len(seq)
        else:
            name, length = None, 0
            for line in handle:
                if line.startswith(">"):
                    if name is not None:
                        lengths[name] = length
                    name, length = line[1:].split()[0], 0
                else:
                    length += len(line.strip())
            if name is not None:
                lengths[name] = length
    return lengths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--out", default="-")
    ap.add_argument("--reads", help="FASTA/FASTQ; required for --frame oriented")
    ap.add_argument("--frame", choices=["oriented", "sequenced"], default="oriented")
    ap.add_argument("--min-coverage", type=int, default=2, help="minimum reads per anchor [2]")
    ap.add_argument("--max-coverage", type=int, default=1000000)
    ap.add_argument("--min-score", type=float, default=0.0, help="per-read assignment phred")
    ap.add_argument("--min-gqn", type=float, default=0.0, help="site GQN; '.' never passes a >0 cut")
    ap.add_argument("--min-explained", type=float, default=0.0)
    ap.add_argument("--nodes", help="restrict to node IDs in this file, one per line")
    ap.add_argument("--keep-shared-names", action="store_true",
                    help="keep reads whose name is carried by more than one alignment. Paired-end "
                         "mates share a name, and a consumer keying on the name would merge two "
                         "different reads, so by default such reads are dropped entirely")
    args = ap.parse_args()

    lengths = read_lengths(args.reads) if args.reads else None
    if args.frame == "oriented" and lengths is None:
        sys.exit("error: --frame oriented needs --reads to mirror positions through the read length")

    # Names carried by more than one alignment, found the only way the file allows: a name in two
    # slots of one site must belong to two reads. Long reads are unpaired and this set is empty.
    shared_names = set()
    if not args.keep_shared_names:
        seen_slot = {}
        site = None
        slot = None
        in_anchor = set()
        with open(args.anchors) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] == "A":
                    site = (f[1], f[2])
                    slot = f[3]
                    in_anchor = set()
                elif f[0] == "R" and site is not None:
                    # Two proofs that a name belongs to two alignments: two slots of one site, or
                    # twice within one anchor. The second is not optional -- it is the half that
                    # catches mates landing on the same allele.
                    if f[1] in in_anchor:
                        shared_names.add(f[1])
                    in_anchor.add(f[1])
                    key = (site, f[1])
                    if key in seen_slot and seen_slot[key] != slot:
                        shared_names.add(f[1])
                    seen_slot[key] = slot
        del seen_slot

    # The interning table: R rows carry an id, and the consumer wants the name.
    id_name = {}
    with open(args.anchors) as handle:
        for line in handle:
            if line.startswith("#read\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) == 3:
                    id_name[f[1]] = f[2]
            elif not line.startswith("#"):
                break
    if not id_name:
        sys.exit(f"error: {args.anchors} has no #read table; expected an anchors-version 4 file")

    keep_nodes = None
    if args.nodes:
        with open(args.nodes) as handle:
            keep_nodes = {int(line.strip()) for line in handle if line.strip()}

    out = []
    current = None
    kept_reads = []
    dropped_site = 0
    dropped_cov = 0
    dropped_score = 0
    dropped_shared = 0

    def flush():
        nonlocal dropped_cov
        if current is None:
            return
        if not (args.min_coverage <= len(kept_reads) <= args.max_coverage):
            dropped_cov += 1
            return
        node, snarl, slot, allele = current
        # The id keeps slot, not allele: slot is what joins to the VCF's GT, and two slots of a
        # heterozygote can carry the same allele only at a site that should have collapsed.
        out.append([f"{snarl}_{node}_{slot}", kept_reads[:]])

    with open(args.anchors) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "A":
                flush()
                kept_reads = []
                node, snarl, slot = int(fields[1]), fields[2], int(fields[3])
                allele = int(fields[4])
                gqn = None if fields[5] == "." else float(fields[5])
                explained = float(fields[6])
                skip = False
                if keep_nodes is not None and node not in keep_nodes:
                    skip = True
                # A '.' GQN means the site offered no gap to normalise, which is not zero and must
                # not pass a positive threshold by accident.
                if args.min_gqn > 0.0 and (gqn is None or gqn < args.min_gqn):
                    skip = True
                if explained < args.min_explained:
                    skip = True
                if skip:
                    dropped_site += 1
                    current = None
                else:
                    current = (node, snarl, slot, allele)
            elif fields[0] == "R" and current is not None:
                read_id = fields[1]
                name = id_name.get(read_id)
                if name is None:
                    sys.exit(f"error: read id {read_id} is not in the #read table")
                strand, offset, score = int(fields[2]), int(fields[3]), float(fields[4])
                if score < args.min_score:
                    dropped_score += 1
                    continue
                if read_id in shared_names:
                    dropped_shared += 1
                    continue
                # The pin, as a position rather than as the base before it.
                pin = offset + 1 if strand == 0 else offset
                if args.frame == "oriented":
                    length = lengths.get(name)
                    if length is None:
                        sys.exit(f"error: read {name} is not in {args.reads}")
                    if strand == 1:
                        pin = length - pin
                kept_reads.append([name, strand, pin, pin])
    flush()

    stream = sys.stdout if args.out == "-" else open(args.out, "w")
    json.dump(out, stream)
    if stream is not sys.stdout:
        stream.close()
    print(f"{len(out)} anchors written; dropped {dropped_site} sites by filter, "
          f"{dropped_cov} by coverage, {dropped_score} read placements by score, "
          f"{dropped_shared} by shared read name ({len(shared_names)} such names)",
          file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())

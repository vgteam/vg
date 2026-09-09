#!/usr/bin/env python3
"""Validate a `vg call --anchors-out` file against the reads it indexes.

`vg` checks the pin invariant in process, against the graph, while the read is still live. This
checks the *written file* against the read file, which is the version that catches a mistake shared
between the collector and its own assertion -- and it is what to run before handing anchors to an
assembler, since the offsets index into the reads as sequenced.

Checks, in order of how badly a failure would mislead:

  1. no (read, strand, offset) appears in more than one anchor -- a pinned position is unique;
  2. every offset lies inside its read;
  3. A rows are ordered by node ID;
  4. the slots of one site have disjoint read sets;
  5. every R row belongs to an A row, its read id resolves in the #read table, and the header
     declares a version this script understands.

Checks 1 and 4 both assume the read NAME identifies one alignment. Paired-end mates share a name, so
on short-read data it does not, and both checks then fire on something `vg` cannot fix: two genuinely
different reads that happen to be called the same thing. Those are separated out and reported as a
warning rather than a failure, with the count, because the distinction matters -- a shared name is a
property of the read file, while a repeated position for one alignment is a bug in the pin geometry.
Long reads are unpaired, so this does not arise there.

**A shared name is only provable when the two alignments touch the same site**, though -- as a name in
two slots of one site, or twice within one anchor. Mates that land on different sites share a name and
are invisible to both proofs, so their position repeats look exactly like a geometry fault from inside
the file. `--paired` says the read file has mates sharing names, which downgrades those to warnings.
Leave it off for long reads, where a repeat really would be a bug.

Reads may be given as FASTA or FASTQ, plain or gzipped. Without them, checks 1, 3, 4 and 5 still run
-- which is most of the value -- and check 2 is skipped rather than silently passed.
"""

import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gzip(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def read_lengths(path):
    """Name -> length, from FASTA or FASTQ. Names are truncated at the first whitespace, which is
    what every aligner does when it writes the name into an alignment record."""
    lengths = {}
    with open_maybe_gzip(path) as handle:
        first = handle.readline()
        if not first:
            return lengths
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
        elif first.startswith(">"):
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
        else:
            sys.exit(f"error: {path} is neither FASTA nor FASTQ")
    return lengths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True, help="the --anchors-out TSV")
    ap.add_argument("--reads", help="FASTA/FASTQ the offsets index into; enables the bounds check")
    ap.add_argument("--quiet", action="store_true", help="print only failures")
    ap.add_argument("--paired", action="store_true",
                    help="the reads are paired-end, so mates share a name. Position repeats that "
                         "cannot be attributed to a detected shared name become warnings instead of "
                         "failures; see the module docstring for why they cannot always be detected")
    args = ap.parse_args()

    lengths = read_lengths(args.reads) if args.reads else None

    # The interning table. R rows carry an integer id; the table above them is the file's only copy
    # of each name, and everything below reports names rather than ids because an id means nothing
    # outside this one file.
    id_name = {}
    with open(args.anchors) as handle:
        for line in handle:
            if line.startswith("#read\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) == 3:
                    id_name[f[1]] = f[2]
            elif not line.startswith("#"):
                break

    # First pass: which names are provably carried by more than one alignment?
    #
    # Two independent proofs, and both are needed. A name in two different SLOTS of one site must
    # belong to two reads. And a name twice within ONE anchor must too: vg walks a site's reads once
    # and contributes at most one row per read per anchor, so a repeat there is a second alignment,
    # not a second visit. Checking only the first proof misses every mate pair that lands on the same
    # allele -- which is half of them, and exactly the half that then shows up as a position repeat.
    shared_names = set()
    seen_slot = {}
    in_anchor = set()
    current = None
    with open(args.anchors) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "A" and len(fields) == 7:
                current = (int(fields[1]), fields[2], int(fields[3]))
                in_anchor = set()
            elif fields[0] == "R" and current is not None and len(fields) == 5:
                if fields[1] in in_anchor:
                    shared_names.add(id_name.get(fields[1], fields[1]))
                in_anchor.add(fields[1])
                key = (current[0], current[1], fields[1])
                previous = seen_slot.get(key)
                if previous is not None and previous != current[2]:
                    shared_names.add(id_name.get(fields[1], fields[1]))
                seen_slot[key] = current[2]
    del seen_slot, in_anchor

    version = None
    anchors = 0
    rows = 0
    current = None
    last_node = None
    seen_position = {}
    site_reads = defaultdict(dict)     # (node, snarl) -> read -> slot
    failures = []
    shared_hits = 0
    undetected_hits = 0

    def fail(message):
        failures.append(message)

    with open(args.anchors) as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if line.startswith("#anchors-version"):
                version = line.split("\t")[1]
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "A":
                if len(fields) != 7:
                    fail(f"line {lineno}: an A row has {len(fields)} fields, expected 7")
                    continue
                node, snarl, slot = int(fields[1]), fields[2], int(fields[3])
                if last_node is not None and node < last_node:
                    fail(f"line {lineno}: node {node} follows {last_node}, so the file is not in "
                         "node order")
                last_node = node
                current = (node, snarl, slot)
                anchors += 1
            elif fields[0] == "R":
                if current is None:
                    fail(f"line {lineno}: an R row before any A row; R rows do not carry their own "
                         "anchor key, so it has nothing to belong to")
                    continue
                if len(fields) != 5:
                    fail(f"line {lineno}: an R row has {len(fields)} fields, expected 5")
                    continue
                name = id_name.get(fields[1])
                if name is None:
                    fail(f"line {lineno}: read id {fields[1]} is not in the #read table")
                    continue
                strand, offset = int(fields[2]), int(fields[3])
                rows += 1

                key = (name, strand, offset)
                if key in seen_position:
                    if name in shared_names:
                        shared_hits += 1
                    elif args.paired:
                        undetected_hits += 1
                    else:
                        fail(f"line {lineno}: {name} strand {strand} offset {offset} is pinned by "
                             f"{current} and already by {seen_position[key]}")
                else:
                    seen_position[key] = current

                if lengths is not None:
                    length = lengths.get(name)
                    if length is None:
                        fail(f"line {lineno}: read {name} is not in the read file")
                    elif not (0 <= offset < length):
                        fail(f"line {lineno}: offset {offset} is outside {name}, which is "
                             f"{length} bp")

                site = (current[0], current[1])
                previous = site_reads[site].get(name)
                if previous is not None and previous != current[2] and name not in shared_names:
                    fail(f"line {lineno}: {name} is in both slot {previous} and slot {current[2]} "
                         f"of {site}")
                site_reads[site][name] = current[2]
            else:
                fail(f"line {lineno}: unknown row type {fields[0]!r}")

    if version is None:
        fail("no #anchors-version header, so the format is unknown")
    elif version != "4":
        # v3 is refused rather than tolerated: it has these columns but its `slot` is in allele
        # order, so it reads cleanly and joins wrongly.
        fail(f"#anchors-version {version} is not the one this script understands (4)")
    if version == "4" and not id_name:
        fail("no #read table, but every version from 2 on interns every read name")

    if shared_names:
        print(f"WARN {len(shared_names)} read names are carried by more than one alignment "
              f"(paired-end mates share a name); {shared_hits} pinned positions repeat because of "
              f"it. Not a pin-geometry fault -- but a consumer that keys on the read name will "
              f"merge two different reads, so rename mates before using this file.", file=sys.stderr)

    if undetected_hits:
        print(f"WARN {undetected_hits} further pinned positions repeat for a read name whose sharing "
              f"could not be proved from the file -- mates landing on different sites never co-occur, "
              f"so neither proof reaches them. Reported because --paired was given; without it these "
              f"are failures.", file=sys.stderr)

    if failures:
        for message in failures[:50]:
            print(f"FAIL {message}", file=sys.stderr)
        if len(failures) > 50:
            print(f"... and {len(failures) - 50} more", file=sys.stderr)
        return 1

    if not args.quiet:
        checked = "offsets bounds-checked" if lengths is not None else "bounds check SKIPPED (no --reads)"
        print(f"OK  {anchors} anchors, {rows} read placements, "
              f"{len(seen_position)} distinct pinned positions, {checked}")
        if shared_names:
            print(f"    ({len(shared_names)} read names shared between alignments; see the warning "
                  f"above)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

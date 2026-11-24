#!/usr/bin/env python3
"""Keep the GAF records whose every node is in a subgraph.

    tabix <sorted.gaf.gz> '{}:FIRST-LAST' | extract_reads.py <node list> > reads.gaf

The whole-genome GAF also holds reads that leave the contig's subgraph, and vg call
aborts on a read whose path names a node its graph does not have. The input stays in
tabix's order, so the output is still sorted for indexing.
"""

import re
import sys

NODE = re.compile(r"[<>]([0-9]+)")


def main() -> None:
    with open(sys.argv[1]) as fh:
        nodes = {int(line) for line in fh}
    kept = dropped = 0
    out = sys.stdout
    for line in sys.stdin:
        path = line.split("\t", 6)[5]
        if all(int(n) in nodes for n in NODE.findall(path)):
            out.write(line)
            kept += 1
        else:
            dropped += 1
    print(f"kept {kept} records, dropped {dropped} that leave the subgraph", file=sys.stderr)


if __name__ == "__main__":
    main()

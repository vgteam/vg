#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 15


# Normal graph without repeated sequence and maximum node length 2418 bp.
# See 54_vg_haplotypes.t for details on the graph.
vg convert -g haplotype-sampling/micb-kir3dl1.gfa > base-graph.vg
vg dedup base-graph.vg -o dedup.vg
is $? 0 "dedup works without repeated sequence"
cmp base-graph.vg dedup.vg > /dev/null
is $? 1 "chopping changes the graph"

# With unique sequence in a long node, vg dedup puts the short piece first,
# while vg mod --chop puts it last. But if we unchop and convert to GFA,
# we will get identical graphs.
vg mod --chop 1024 base-graph.vg > chopped.vg
is "$(vg stats -z -l chopped.vg)" "$(vg stats -z -l dedup.vg)" "chopped graph has the same stats as with vg mod --chop"
vg mod --unchop dedup.vg | vg convert -f - > dedup-unchopped.gfa
vg mod --unchop chopped.vg | vg convert -f - > chopped-unchopped.gfa
cmp dedup-unchopped.gfa chopped-unchopped.gfa
is $? 0 "unchopped graphs are identical"

# The graph should simply pass through if there are no long nodes.
vg dedup -m 3000 base-graph.vg -o no-chopping.vg
is $? 0 "dedup works when there is no need for chopping"
cmp base-graph.vg no-chopping.vg
is $? 0 "graph passes through without changes"

# Cleanup
rm -f base-graph.vg dedup.vg chopped.vg dedup-unchopped.gfa chopped-unchopped.gfa no-chopping.vg


# A test graph with repeated sequence that can be deduplicated.
# See dedup/README.md for further details.
gunzip -c dedup/dedup.gfa.gz > base-graph.gfa
vg convert -g base-graph.gfa > base-graph.vg

# Just the graph output.
vg dedup base-graph.vg -o dedup.vg
is $? 0 "dedup works with repeated sequence"
cmp base-graph.vg dedup.vg > /dev/null
is $? 1 "chopping changes the graph"
rm -f dedup.vg

# GFA and FASTA output.
# TODO: Are the numbers stable across platforms?
vg dedup base-graph.vg -f dedup.fa -G dedup.gfa
is $? 0 "dedup works with GFA and FASTA output"
is "$(grep -c '^>' dedup.fa)" 1716 "1716 stored strings in FASTA"
is "$(grep -c '^S\t[0-9]*\t\*' dedup.gfa)" 12488 "12488 chopped pieces in GFA"
rm -f dedup.fa dedup.gfa

# All three outputs.
vg dedup base-graph.vg -o dedup.vg -f dedup.fa -G dedup.gfa
is $? 0 "dedup works with all three outputs"
cmp base-graph.vg dedup.vg > /dev/null
is $? 1 "chopping changes the graph"
is "$(grep -c '^>' dedup.fa)" 1716 "1716 stored strings in FASTA"
is "$(grep -c '^S\t[0-9]*\t\*' dedup.gfa)" 12488 "12488 chopped pieces in GFA"
rm -f dedup.vg dedup.fa dedup.gfa

# Cleanup.
rm -f base-graph.vg base-graph.gfa


# TODO: GBZ v4 from GFA + FASTA

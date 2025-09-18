#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 16

# See 54_vg_haplotypes.t for details on the graph.

# Build a graph and a distance index and find snarls.
vg gbwt -G haplotype-sampling/micb-kir3dl1.gfa --gbz-format -g graph.gbz
vg index -j graph.dist graph.gbz
vg snarls --include-trivial graph.gbz > graph.snarls

# Distance index, binary format
vg chains graph.gbz graph.dist > dist.stdout.binary
is $? 0 "chains from a distance index to stdout in binary format"
vg chains graph.gbz graph.dist -o dist.file.binary
is $? 0 "chains from a distance index to a file in binary format"
cmp dist.file.binary dist.stdout.binary
is $? 0 "the chain files are identical"

# Snarls, binary format
vg chains graph.gbz graph.snarls > snarls.stdout.binary
is $? 0 "chains from a snarls file to stdout in binary format"
vg chains graph.gbz graph.snarls -o snarls.file.binary
is $? 0 "chains from a snarls file to a file in binary format"
cmp snarls.file.binary snarls.stdout.binary
is $? 0 "the chain files are identical"

# Note that the files are not necessarily identical, if there are multiple
# possible start/end nodes for a chain (see chains_main.cpp).
cmp snarls.stdout.binary dist.stdout.binary
is $? 0 "chain files from different inputs are identical"
is $(md5sum snarls.stdout.binary | cut -f1 -d' ') cff1d755216134420c67115ee42160ed "binary output is correct" 

# Distance index, GFA format
vg chains graph.gbz graph.dist --gfa > dist.stdout.gfa
is $? 0 "chains from a distance index to stdout in GFA format"
vg chains graph.gbz graph.dist --gfa -o dist.file.gfa
is $? 0 "chains from a distance index to a file in GFA format"
cmp dist.file.gfa dist.stdout.gfa
is $? 0 "the chain files are identical"

# Snarls, GFA format
vg chains graph.gbz graph.snarls --gfa > snarls.stdout.gfa
is $? 0 "chains from a snarls file to stdout in GFA format"
vg chains graph.gbz graph.snarls --gfa -o snarls.file.gfa
is $? 0 "chains from a snarls file to a file in GFA format"
cmp snarls.file.gfa snarls.stdout.gfa
is $? 0 "the chain files are identical"

cmp snarls.stdout.gfa dist.stdout.gfa
is $? 0 "chain files from different inputs are identical"
is $(md5sum snarls.stdout.gfa | cut -f1 -d' ') a09ce9aa7125aff93e01cfeb76d026f2 "GFA output is correct"

# Cleanup
rm -f graph.gbz graph.dist graph.snarls
rm -f dist.stdout.binary dist.file.binary
rm -f snarls.stdout.binary snarls.file.binary
rm -f dist.stdout.gfa dist.file.gfa
rm -f snarls.stdout.gfa snarls.file.gfa

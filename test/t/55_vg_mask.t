#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 4

printf "GRCh38#0#chr1\t1\t4\n" > m.bed
vg gbwt -g m.gbz --gbz-format -G graphs/gfa_with_reference.gfa
vg mask -g -b m.bed m.gbz > mm.gbz

is $(vg view mm.gbz | grep "^S" | grep 4 | cut -f 3) "NNN" "A node can be masked in a GBZ"
is $(vg view mm.gbz | grep "^S" | grep -v 4 | grep N | wc -l | sed 's/^[[:space:]]*//') "0"  "Off target nodes are not masked in a GBZ"

vg convert -g graphs/gfa_with_reference.gfa > m.vg
vg mask -b m.bed m.vg > mm.vg

is $(vg view mm.vg | grep "^S" | grep 4 | cut -f 3) "NNN" "A node can be masked in a BDSG graph"
is $(vg view mm.vg | grep "^S" | grep -v 4 | grep N | wc -l | sed 's/^[[:space:]]*//') "0"  "Off target nodes are not masked in a BDSG graph"

rm m.bed m.gbz mm.gbz m.vg mm.vg

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 0
exit

# We have broken the index format that this was using
# I think it makes more sense to disable this and rewrite it using GCSA2 than continue maintaining the old index just for it


# Handmade-example where graphs share 3 kmers.

# Import the graphs
vg view -J -v compare/graph1.json > graph1.vg
vg view -J -v compare/graph2.json > graph2.vg

# Index 6mers
vg index -k 6 -d graph1.idx graph1.vg
vg index -k 6 -d graph2.idx graph2.vg

# compare
vg compare graph1.idx graph2.idx  > comparison.json

is $(jq --argfile a comparison.json --argfile b compare/truth.json -n '($a | (.. | arrays) |= sort) as $a | ($b | (.. | arrays) |= sort) as $b | $a == $b') true "vg compare produces the expected output for 6mer comparison test case."

rm -rf graph1.vg graph2.vg graph1.idx graph2.idx comparison.json


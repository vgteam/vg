#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

# Handmade-example where graphs share 3 kmers.

# Import the graphs
vg view -J -v compare/graph1.json > graph1.vg
vg view -J -v compare/graph2.json > graph2.vg

# Index 6mers
vg index -k 6 graph1.vg
vg index -k 6 graph2.vg

# compare
vg compare graph1.vg graph2.vg  > comparison.json

is $(jq --argfile a comparison.json --argfile b compare/truth.json -n '($a | (.. | arrays) |= sort) as $a | ($b | (.. | arrays) |= sort) as $b | $a == $b') true "vg compare produces the expected output for 6mer comparison test case."

rm -rf graph1.vg graph2.vg graph1.vg.index graph2.vg.index comparison.json


#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 5

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg

vg simplify --algorithm small x.vg > x.small.vg
is "${?}" "0" "vg simplify runs through when popping small bubbles"

is "$(vg mod --unchop x.small.vg | vg stats -N -)" "1" "simplification pops all the bubbles in a simple graph"

vg simplify --algorithm rare --min-count 2 -v small/x.vcf.gz x.vg > x.rare.vg
is "${?}" "0" "vg simplify runs through when removing rare variants"

vg validate x.rare.vg
is "${?}" "0" "the graph is valid after removing rare variants"

is "$(vg mod --unchop x.rare.vg | vg stats -N -)" "109" "simplification keeps only some variants"

rm -f x.vg x.small.vg x.rare.vg
 


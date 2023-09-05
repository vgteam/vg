#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 5

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg

# Make sure to discard all the warnings about having removed alt paths.
vg simplify --algorithm small x.vg > x.small.vg 2>/dev/null
is "${?}" "0" "vg simplify runs through when popping small bubbles"

is "$(vg paths -d -v x.small.vg | vg mod --unchop - | vg stats -N -)" "1" "simplification pops all the bubbles in a simple graph"

vg simplify --algorithm rare --min-count 2 -v small/x.vcf.gz x.vg > x.rare.vg
is "${?}" "0" "vg simplify runs through when removing rare variants"

vg validate x.rare.vg
is "${?}" "0" "the graph is valid after removing rare variants"

is "$(vg paths -d -v x.rare.vg | vg mod --unchop - | vg stats -N -)" "118" "simplification keeps only some variants"

rm -f x.vg x.small.vg x.rare.vg
 


#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 1

# Toy example of hand-made pileup to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg view -J -l call/pileup.json -L | vg call - -j > snps.gam.json
is $(jq --argfile a snps.gam.json --argfile b call/truth.json -n '($a | (.. | arrays) |= sort) as $a | ($b | (.. | arrays) |= sort) as $b | $a == $b') true "vg call produces the expected output for toy snp test case."
rm -f snps.gam.json 


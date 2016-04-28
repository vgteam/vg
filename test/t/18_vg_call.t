#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 2

# Toy example of hand-made pileup (and hand inspected truth) to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg view -J -l call/pileup.json -L | vg call tiny.vg - -j -s 10 -d 10 -q 10 -b 0.5 > calls.json
is $(jq --argfile a calls.json --argfile b call/truth.json -n '($a == $b)') true "vg call produces the expected output for toy snp test case."

vg view -J -l call/pileup.json -L | vg call tiny.vg - -j -s 10 -d 10 -q 10 -l 0.5 > calls_l.json
is $(jq --argfile a calls_l.json --argfile b call/truth_l.json -n '($a == $b)') true "vg call -l produces the expected output for toy snp test case."

rm -f calls_l.json calls.json tiny.vg



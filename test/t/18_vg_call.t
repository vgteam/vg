#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 2

# Toy example of hand-made pileup (and hand inspected truth) to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg

vg view -J -l call/pileup.json -L | vg call tiny.vg - -A calls_l.vg -s 10 -d 10 -q 10 -b 0.25 -f 0.2 -I > /dev/null 2> /dev/null
vg view -j calls_l.vg | jq . > calls_l.json
is $(jq --argfile a calls_l.json --argfile b call/truth_l.json -n '($a == $b)') true "vg call -l produces the expected output for toy snp test case."

rm -f calls_l.json calls_l.vg

vg view -J -l call/pileup.json -L > tiny.vgpu
vg call tiny.vg tiny.vgpu -A calls_l.vg -s 10 -d 10 -q 10 -b 0.25  > /dev/null 2> /dev/null
is $? "0" "vg call doesn't crash in multiallelic mode"

rm -f calls_l.json calls_l.vg tiny.vg tiny.vgpu



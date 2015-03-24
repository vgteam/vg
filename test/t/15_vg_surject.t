#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 2

vg construct -r small/x.fa >j.vg
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s -k 11 x.vg

is $(vg map -r <(vg sim -s 1337 -n 100 j.vg) x.vg | vg surject -p x -d x.vg.index - | vg view -a - | jq .score | grep -v 200 | wc -l) \
    0 "vg surject works perfectly for perfect reads derived from the reference"

is $(vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg surject -p x -d x.vg.index - | vg view -a - | wc -l) \
    100 "vg surject works for every read simulated from a dense graph"

rm -rf j.vg x.vg x.vg.index

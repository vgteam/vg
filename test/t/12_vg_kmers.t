#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 1

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | wc -l) \
    151216 \
    "correct numbers of kmers in the graph"

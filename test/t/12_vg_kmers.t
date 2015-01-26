#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 4

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | sort | uniq | wc -l) \
    7141 \
    "correct numbers of kmers in the graph"

is  $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | sort | uniq | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz -t 4 | vg kmers -k 11 - | wc -l) \
    "only unique kmers are produced"

vg construct -v small/x.vcf.gz -r small/x.fa >x.vg
vg index -s x.vg
is $(vg find -n 10 -c 1 x.vg | vg kmers -k 11 - | sort -n -k 2 -k 3 | tail -1 | cut -f 2) \
    14 \
    "kmers correctly generated from all nodes"
rm x.vg
rm -rf x.vg.index

is $(vg kmers -k 11 -e 7 jumble/j.vg | wc -l) \
    9300 \
    "edge-max correctly bounds the number of kmers in a complex graph"

#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 7

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | sort | uniq | wc -l) \
    7141 \
    "correct numbers of kmers in the graph"

is  $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | sort | uniq | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz -t 4 | vg kmers -k 11 - | wc -l) \
    "only unique kmers are produced"

vg construct -v small/x.vcf.gz -r small/x.fa >x.vg
vg index -s x.vg
is $(vg find -n 10 -c 1 x.vg | vg kmers -k 11 - | sort -n -k 2 -k 3 | tail -1 | cut -f 2) \
    12 \
    "kmers correctly generated from all nodes"

is $(vg kmers -g -k 11 -t 1 x.vg | wc -l) 2168 "GCSA2 output produces the expected number of lines"

#is $(vg kmers -g -k 11 -t 1 x.vg | cut -f 1 | sort | uniq | wc -l) $(vg kmers -k 11 x.vg | cut -f 1 | sort | uniq | wc -l) "GCSA2 produces output for all kmers"

is $(vg kmers -g -k 11 -t 1 x.vg | grep AATAAGGCTTG | md5sum | cut -f 1 -d\ ) "72e6700f7a6906d2d34e3bf12de78e9f" "GCSA2 output works when next position is multiple"

is $(vg kmers -g -k 11 -t 1 x.vg | grep CATATTAGCCA | md5sum | cut -f 1 -d\ ) "2bbb55f269882959418f9f55e651cd2a" "GCSA2 output works when previous characters are multiple"

rm x.vg
rm -rf x.vg.index

is $(vg kmers -k 11 -e 7 jumble/j.vg | wc -l) \
    9300 \
    "edge-max correctly bounds the number of kmers in a complex graph"

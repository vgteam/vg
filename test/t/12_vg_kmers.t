#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 12

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | sort | uniq | wc -l) \
    2133 \
    "correct numbers of kmers in the graph when reporting only start positions"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -n -k 11 - | sort | uniq | wc -l) \
    7141 \
    "correct numbers of kmers in the graph when including negative-offset entries"

is  $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 -d - | sort | uniq | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz -t 4 | vg kmers -k 11 -d - | wc -l) \
    "only unique kmers are produced"
    
is $(vg kmers -k 15 reversing/reversing_edge.vg | grep "CAAATAAGTGTAATC" | wc -l) 1 "to_end edges are handled correctly"

is $(vg kmers -k 15 reversing/reversing_edge.vg | grep "AAATAAGTGTAATCA" | wc -l) 1 "from_start edges are handled correctly"

vg construct -v small/x.vcf.gz -r small/x.fa >x.vg
vg index -s x.vg
is $(vg find -n 10 -c 1 x.vg | vg kmers -n -k 11 - | sort -n -k 2 -k 3 | tail -1 | cut -f 2) \
    12 \
    "kmers correctly generated from all nodes"

is $(vg kmers -g -k 11 -t 1 x.vg | wc -l) 2168 "GCSA2 output produces the expected number of lines"

#is $(vg kmers -g -k 11 -t 1 x.vg | cut -f 1 | sort | uniq | wc -l) $(vg kmers -k 11 x.vg | cut -f 1 | sort | uniq | wc -l) "GCSA2 produces output for all kmers"

is $(vg kmers -g -k 11 -t 1 x.vg | grep AATAAGGCTTG | md5sum | cut -f 1 -d\ ) "72e6700f7a6906d2d34e3bf12de78e9f" "GCSA2 output works when next position is multiple"

is $(vg kmers -g -k 11 -t 1 x.vg | grep CATATTAGCCA | md5sum | cut -f 1 -d\ ) "2bbb55f269882959418f9f55e651cd2a" "GCSA2 output works when previous characters are multiple"

is $(vg construct -r small/x.fa -v small/x.vcf.gz| vg kmers -g -k 11 -t 1 - | grep AAGAATACAA | md5sum | cut -f 1 -d\ ) "b56ea597d9f876f99d24e25fe0c710c1" "GCSA2 output correctly represents repeated kmers at the same position"

rm x.vg
rm -rf x.vg.index

is $(vg kmers -n -k 11 -e 5 -d jumble/j.vg | wc -l) \
    12354 \
    "edge-max correctly bounds the number of kmers in a complex graph"

is $(vg construct -r small/x.fa -v small/x.vcf.gz| vg kmers -g -k 11 -t 1 -H 1000 -T 1001 - | grep '1000\|1001' | wc -l) 37 "head and tail nodes can be specified in GCSA2 output"

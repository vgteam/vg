#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 15

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

is $(vg kmers -g -k 11 -t 1 x.vg | wc -l) 4342 "GCSA2 output produces the expected number of lines"

is $(vg kmers -gB -k 11 -t 1 x.vg | wc -c) 111528 "GCSA2 binary output produces the expected number bytes"

#is $(vg kmers -g -k 11 -t 1 x.vg | cut -f 1 | sort | uniq | wc -l) $(vg kmers -k 11 x.vg | cut -f 1 | sort | uniq | wc -l) "GCSA2 produces output for all kmers"

is "$(vg kmers -g -k 11 -t 1 x.vg | grep AATAAGGCTTG | cut -f 4,5)" "$(printf 'A,G\t7:0,8:0')" "GCSA2 output works when next position is multiple"

is "$(vg kmers -g -k 11 -t 1 x.vg | grep CATATTAGCCA | cut -f 3)" "A,G" "GCSA2 output works when previous characters are multiple"

is "$(vg kmers -g -k 11 -t 1 x.vg | grep AAGAATACAA | cut -f1,3,4 | tr '\n\t' '  ')" "AAAGAATACAA G A,G AAGAATACAAA A G AAGAATACAAG A A " "GCSA2 output correctly represents repeated kmers at the same position"

rm x.vg
rm -rf x.vg.index

is "$(vg kmers -k 15 -e 1 jumble/j.vg -t 1 | cut -f1 | sort | tr '\n' ' ')" "CGGCCTGGCGCACAA CGGCCTGGCTCACAA TGGCCTGGCGCACAA TGGCCTGGCTCACAA " "edge-max can limit to paths with exactly one choice from any node's point of view"

is $(vg kmers -n -k 11 -e 2 -d jumble/j.vg | wc -l) \
    7096 \
    "edge-max correctly bounds the number of kmers in a complex graph"


is $(vg construct -r small/x.fa -v small/x.vcf.gz| vg kmers -g -k 11 -t 1 -H 1000 -T 1001 - | grep '1000\|1001' | wc -l) 76 "start/stop node IDs can be specified in GCSA2 output"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view - |head -10 | vg view -v - | vg mod -o - | vg kmers -k 16 - >/dev/null
is $? 0 "attempting to generate kmers longer than the longest path in a graph correctly yields no kmers"

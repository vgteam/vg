#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 10

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg kmers -k 11 - | cut -f 1 | sort | uniq | wc -l) \
    4250 \
    "correct numbers of kmers in the graph"
    
is $(vg kmers -k 15 reversing/reversing_edge.vg | grep "CAAATAAGTGTAATC" | wc -l) 1 "to_end edges are handled correctly"

is $(vg kmers -k 15 reversing/reversing_edge.vg | grep "AAATAAGTGTAATCA" | wc -l) 1 "from_start edges are handled correctly"

vg construct -v small/x.vcf.gz -r small/x.fa >x.vg

is $(vg kmers -g -k 11 -t 1 x.vg | wc -l) 4356 "GCSA2 output produces the expected number of lines"

is $(vg kmers -gB -k 11 -t 1 x.vg | wc -c) 111528 "GCSA2 binary output produces the expected number bytes"

#is $(vg kmers -g -k 11 -t 1 x.vg | cut -f 1 | sort | uniq | wc -l) $(vg kmers -k 11 x.vg | cut -f 1 | sort | uniq | wc -l) "GCSA2 produces output for all kmers"

is "$(vg kmers -g -k 11 -t 1 x.vg | grep AATAAGGCTTG | cut -f 4,5)" "$(printf 'A,G\t7:0,8:0')" "GCSA2 output works when next position is multiple"

is "$(vg kmers -g -k 11 -t 1 x.vg | grep CATATTAGCCA | cut -f 3)" "G,A" "GCSA2 output works when previous characters are multiple"

rm x.vg
rm -rf x.vg.index

vg mod -p -l 11 -e 8 -S jumble/j.vg > j.vg
is $(vg kmers -g -k 11 j.vg | wc -l) \
   26148 \
   "edge-max correctly bounds the number of kmers in a complex graph"
rm -f j.vg


is $(vg construct -r small/x.fa -v small/x.vcf.gz| vg kmers -g -k 11 -t 1 -H 1000 -T 1001 - | grep '1000\|1001' | wc -l) 76 "start/stop node IDs can be specified in GCSA2 output"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view - | head -10 | vg view -vF - | vg mod -o - | vg kmers -k 16 - >/dev/null
is $? 0 "attempting to generate kmers longer than the longest path in a graph correctly yields no kmers"

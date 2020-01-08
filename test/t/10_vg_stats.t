#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 11

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
#is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

nodes=$(vg stats -z z.vg | head -1 | cut -f 2)
real_nodes=$(vg view -j z.vg | jq -c '.node[]' | wc -l)
is $nodes $real_nodes "vg stats reports the expected number of nodes"

edges=$(vg stats -z z.vg | tail -1 | cut -f 2)
real_edges=$(vg view -j z.vg | jq -c '.edge[]' | wc -l)
is $edges $real_edges "vg stats reports the expected number of edges"

graph_length=$(vg stats -l z.vg | tail -1 | cut -f 2)
real_length=$(vg view -j z.vg | jq -r '.node[].sequence' | tr -d '\n' | wc -c)
is $graph_length $real_length "vg stats reports the expected graph length"

subgraph_count=$(vg stats -s z.vg | wc -l)
is $subgraph_count 1 "vg stats reports the correct number of subgraphs"

subgraph_length=$(vg stats -s z.vg | head -1 | cut -f 2)
is $subgraph_length $graph_length  "vg stats reports the correct subgraph length"

rm -f z.vg

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg
is $(vg stats -n 13 -d t.vg | cut -f 2) 38 "distance to head is correct"
is $(vg stats -n 13 -t t.vg | cut -f 2) 11 "distance to tail is correct"
rm -f t.vg

vg construct -m 1000 -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg map -x x.xg -g x.gcsa -T small/x-s1337-n100.reads >x.gam
is "$(vg stats -a x.gam x.vg | md5sum | cut -f 1 -d\ )" "$(md5sum correct/10_vg_stats/15.txt | cut -f 1 -d\ )" "aligned read stats are computed correctly"

is "$(vg stats -z x.vg)" "$(vg stats -z x.xg)" "basic stats agree between graph formats"

is "$(vg stats -a x.gam | grep 'Total alignments')" "Total alignments: 100" "stats can be computed for GAM files without graphs"
rm -f x.vg x.xg x.gcsa x.gam

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view -g - > tiny_names.gfa
printf "P\tref.1\t1+,3+,5+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
printf "P\talt1.1\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
vg view -Fv tiny_names.gfa > tiny_names.vg 
is $(vg stats -O tiny_names.vg | wc -l) 113 "a path overlap description of a test graph has the expected length"
rm -f tiny_names.gfa tiny_names.vg

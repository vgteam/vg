#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 9

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

is $(vg view -Fv msgas/q_redundant.gfa | vg stats -S - | md5sum | cut -f 1 -d\ ) 01fadb6a004ddb87e5fc5d056b565218 "perfect to and from siblings are determined"

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg
is $(vg stats -n 13 -d t.vg | cut -f 2) 38 "distance to head is correct"
is $(vg stats -n 13 -t t.vg | cut -f 2) 11 "distance to tail is correct"

vg construct -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg sim -s 1337 -n 100 -x x.xg >x.reads
vg map -x x.xg -g x.gcsa -r x.reads >x.gam
is "$(vg stats -a x.gam x.vg | md5sum | cut -f 1 -d\ )" "$(md5sum correct/10_vg_stats/15.txt | cut -f 1 -d\ )" "aligned read stats are computed correctly"
rm -f x.vg x.xg x.gcsa x.gam x.reads


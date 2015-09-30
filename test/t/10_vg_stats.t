#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 6

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
#is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

nodes=$(vg stats -z z.vg | head -1 | cut -f 2)
is $nodes 84553 "vg stats reports the expected number of nodes"

edges=$(vg stats -z z.vg | tail -1 | cut -f 2)
is $edges 115357 "vg stats reports the expected number of edges"

graph_length=$(vg stats -l z.vg | tail -1 | cut -f 2)
is $graph_length 1029257 "vg stats reports the expected graph length"

subgraph_count=$(vg stats -s z.vg | wc -l)
is $subgraph_count 1 "vg stats reports the correct number of subgraphs"

subgraph_length=$(vg stats -s z.vg | head -1 | cut -f 2)
is $subgraph_length $graph_length  "vg stats reports the correct subgraph length"

rm -f z.vg

is $(vg view -v msgas/q_redundant.gfa | vg stats -S - | md5sum | cut -f 1 -d\ ) 01fadb6a004ddb87e5fc5d056b565218 "perfect to and from siblings are determined"

#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 6

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg stats -z - | grep nodes | cut -f 2) 210 "construction produces the right number of nodes"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg stats -z - | grep edges | cut -f 2) 291 "construction produces the right number of edges"

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

size=$(ls -s z.vg | cut -f 1 -d\ )
is $size 2740 "the 1mb graph has the expected file size"

nodes=$(vg stats -z z.vg | head -1 | cut -f 2)
is $nodes 84646 "the 1mb graph has the expected number of nodes"

edges=$(vg stats -z z.vg | tail -1 | cut -f 2)
is $edges 115515 "the 1mb graph has the expected number of edges"

rm -f z.vg

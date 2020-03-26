#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 4

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg

num_nodes=$(vg view -g x.vg | grep ^S | wc -l)
num_edges=$(vg view -g x.vg | grep ^L | wc -l)

#echo $num_nodes

is $(vg concat x.vg x.vg | vg view -g - | grep ^S | wc -l) $(echo "$num_nodes * 2" | bc) "concat doubles the number of nodes"
is $(vg concat x.vg x.vg | vg view -g - | grep ^L | wc -l) $(echo "$num_edges * 2 + 1" | bc) "concat doubles the number of edges + 1"

rm -f x.vg

vg view -Jv ./reversing/reversing_path.json  > reversing.vg

num_nodes=$(vg view -g reversing.vg | grep ^S | wc -l)
num_edges=$(vg view -g reversing.vg | grep ^L | wc -l)

is $(vg concat reversing.vg reversing.vg -p | vg view -g - | grep ^S | wc -l) $(echo "$num_nodes * 2" | bc) "concat -p doubles the number of nodes on reversing graph"
# without -p, the heads/tails are backwards so you get something uglier
is $(vg concat reversing.vg reversing.vg -p | vg view -g - | grep ^L | wc -l) $(echo "$num_edges * 2 + 1" | bc) "concat -p doubles the number of edges + 1 on reversing graph"

rm -f reversing.vg


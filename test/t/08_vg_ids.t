#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 8

num_nodes=$(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -c - | vg view -g - | grep ^S | wc -l)

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -i 1000 - | vg ids -c - | vg view -g - | grep ^S | cut -f 2 | head -1) 1 "minimum node as expected (topological sort and id compaction correct)"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -i 1000 - | vg ids -c - | vg view -g - | grep ^S | cut -f 2 | tail -1) $num_nodes "maximum node is as expected (topological sort and id compaction correct)"

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg ids -j x.vg y.vg z.vg

last=$(vg view -g y.vg | grep ^S | tail -1 | cut -f 2)
first=$(vg view -g z.vg | grep ^S | head -1 | cut -f 2)

is $first $(echo "$last + 1" | bc) "correctly generated joint id space for several graphs"

rm x.vg y.vg z.vg

vg ids -s cyclic/self_loops.vg > sorted.vg
is $? 0 "can sort and re-number a graph with self loops"

vg ids -s cyclic/all.vg > sorted.vg
is $? 0 "can sort and renumber a complex cyclic graph"

is $(vg ids -s ids/unordered.vg | vg view -j - | jq -c '.edge[] | select(.from > .to)' | wc -l) 0 "sorting removes back-edges in a DAG"

rm sorted.vg

is $(vg ids -s ids/unordered.vg | vg view -j - | jq -c '.node[1] == {"id":2,"sequence":"T"}') "true" "sorting assigns node IDs in topological order"

vg ids -s graphs/snp1kg-brca2-unsorted.vg | vg validate -
is $? 0 "can handle graphs with out-of-order mappings"

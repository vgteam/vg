#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10

num_nodes=$(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -c - | vg view -g - | grep ^S | wc -l)

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -i 1000 - | vg ids -c - | vg view -g - | grep ^S | cut -f 2 | sort -n | head -1) 1 "minimum node as expected (id compaction correct)"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -i 1000 - | vg ids -c - | vg view -g - | grep ^S | cut -f 2 | sort -n | tail -1) $num_nodes "maximum node is as expected (id compaction correct)"

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg ids -j x.vg y.vg z.vg

last=$(vg view -g y.vg | grep ^S | cut -f 2 | sort -n | tail -1)
first=$(vg view -g z.vg | grep ^S | cut -f 2 | sort -n | head -1)

is $first $(echo "$last + 1" | bc) "correctly generated joint id space for several graphs"

rm x.vg y.vg z.vg

vg ids -s cyclic/self_loops.vg > sorted.vg
is $? 0 "can sort and re-number a graph with self loops"

vg ids -s cyclic/all.vg > sorted.vg
is $? 0 "can sort and renumber a complex cyclic graph"

is $(vg ids -s ids/unordered.vg | vg view -j - | jq -r -c '.edge[] | select((.from | tonumber) > (.to | tonumber))' | wc -l) 0 "sorting removes back-edges in a DAG"

rm sorted.vg

is $(vg ids -s ids/unordered.vg | vg view -j - | jq -r -c '.node[1] | (.sequence == "T" and (.id | tostring) == "2")') "true" "sorting assigns node IDs in topological order"

is $(vg ids -s ids/unordered.vg | vg stats -r - | awk '{print $2}') $(vg stats -r ids/unordered.vg | awk '{print $2}') "sorting does not affect id range"
is $(vg convert ids/unordered.vg -v | vg ids -s - | vg stats -r - | awk '{print $2}') $(vg stats -r ids/unordered.vg | awk '{print $2}') "sorting does not affect id range of vg"
is $(vg convert ids/unordered.vg -a | vg ids -s - | vg stats -r - | awk '{print $2}') $(vg stats -r ids/unordered.vg | awk '{print $2}') "sorting does not affect id range of hg"
# this test relies on id sorting being implemented in pg
#is $(vg convert ids/unordered.vg -p | vg ids -s - | vg stats -r - | awk '{print $2}') $(vg stats -r ids/unordered.vg | awk '{print $2}') "sorting does not affect id range of pg"

# this test now breaks under the current VG.paths semantics, which require our paths to record the exact match lengths of the nodes
#vg ids -s graphs/snp1kg-brca2-unsorted.vg | vg validate -
#is $? 0 "can handle graphs with out-of-order mappings"

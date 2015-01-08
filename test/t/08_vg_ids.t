#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 3

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

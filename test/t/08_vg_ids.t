#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 2

num_nodes=$(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -c - | vg view -g - | grep ^S | wc -l)

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -i 1000 - | vg ids -c - | vg view -g - | grep ^S | cut -f 2 | head -1) 1 "minimum node as expected (topological sort and id compaction correct)"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg ids -i 1000 - | vg ids -c - | vg view -g - | grep ^S | cut -f 2 | tail -1) $num_nodes "maximum node is as expected (topological sort and id compaction correct)"


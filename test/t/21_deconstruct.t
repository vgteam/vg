#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

## Test VG .to_superbubbles()
is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg deconstruct -s - | wc -l) 4 "vg deconstruct produces the expected number of superbubbles in a simple graph."

## Make sure deconstruct successfully finds the right nodes in a superbubble of a simple graph.
#is $(vg construct -r tiny/tiny.fa -v tiny.vcf.gz | vg deconstruct -s - > x.txt && diff x.txt expected_x.txt | wc -l) 0 "vg deconstruct produces the expected superbubble format for a simple graph."

## Test if deconstruct works on a graph that must be DAGified.
#is $(vg construct -r COMPLEXGRAPH -v ANOTHERONE -m 50 | vg deconstruct -s - | md5sum) SUPERBUBBLEFILEHASH "deconstruct finds the expected superbubbles in a DAGified graph and uses the IDs from the original graph space."

## Test if we get the expected VCF from all paths in graph
#is $(vg construct -r tiny/tiny.fa -v tiny.vcf.gz | vg deconstruct - | md5sum) aaee5555 "vg deconstruct produces the expected vcf output on a simple graph."

## Test deconstruction on a specific path
# is $(vg construct -r multipath/mp.fa -v multipath/mp.vcf.gz | vg deconstruct -p - | md5sum AFJAKJFIJ "deconstruct can deconstruct a specific path."

## Test deconstruction on a set of paths specified in a file
# is $(vg construct -r multipath/mp.fa -v multipath/mp.vcf.gz | vg deconstruct -P - | md5sum AFJAKJFIJ "deconstruct can deconstruct multiple paths specified in a file."

## Test if the depth filtering for deconstruct works with a small GAM file
# is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg deconstruct -d 10 -a tiny/alignment.gam - | md5sum) AFKJADK "deconstruct can depth filter using aGAM."

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 8

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -t 1 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 549b183e90036f767acdd10e7d5ba125 "MSGA produces the expected graph for GRCh38 HLA-V"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -t 4 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 549b183e90036f767acdd10e7d5ba125 "graph for GRCh38 HLA-V is unaffected by the number of alignment threads"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 10000 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 300 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ )  "varying alignment bandwidth does not affect output graph"

is $(vg msga -f msgas/s.fa -k 16 -b s1 -B 20 | vg view - | md5sum | cut -f 1 -d\ ) d8a7332e4cca570e0d8325606e6a1e25 "banded alignment can detect and include large deletions in the graph"

is $(vg msga -f msgas/s.fa -k 16 -b s2 -B 20 | vg view - | md5sum | cut -f 1 -d\ ) d8a7332e4cca570e0d8325606e6a1e25 "assembly around indels works regardless of the base sequence that is used"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
is $(vg msga -g t.vg -s CAAATTTTCTGGAGTTCTAT | vg stats -s - | wc -l) 1 "soft clips at node boundaries (start) are included correctly"
is $(vg msga -g t.vg -s TTCTATAATATG | vg stats -s - | wc -l) 1 "soft clips at node boundaries (end) are included correctly"
rm t.vg

is $(vg msga -f msgas/s.fa -k 16 -b s1 -B 20 | vg view -j - | jq -M -c --sort-keys '{"node": .node, "edge": .edge}') $(vg msga -f msgas/s.fa -f msgas/s-rev.fa -k 16 -b s1 -B 20 | vg view -j - | jq -M -c --sort-keys '{"node": .node, "edge": .edge}') "adding in existing sequences in reverse doesn't change graph"

# TODO: you should get the same graph no matter the orientation of the input sequences, but you don't. This is probably due to the orientation-dependence of the kmer search.


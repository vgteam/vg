#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 15

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 | vg mod -U 10 - | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 16e56f0090b310d2b1479d49cf790324 "MSGA produces the expected graph for GRCh38 HLA-V"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 4 | vg mod -U 10 - | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 16e56f0090b310d2b1479d49cf790324 "graph for GRCh38 HLA-V is unaffected by the number of alignment threads"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 --no-mem-threader| vg mod -U 10 - | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 5456c2fb7c9a9c80f6a4cfd4dcc8d4ec "MSGA produces the expected graph for GRCh38 HLA-V with the MEM threader off"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 10000 -N | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 300 -N | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ )  "varying alignment bandwidth does not affect output graph"

is $(vg msga -f msgas/s.fa -B 20 | vg mod -U 10 - | vg mod -c - | vg view - | md5sum | cut -f 1 -d\ ) ceadc1a150fb2cafee88288c1bf45cc6 "msga alignment can detect and include large deletions in the graph"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
is $(vg msga -g t.vg -s CAAATTTTCTGGAGTTCTAT -N | vg stats -s - | wc -l) 1 "soft clips at node boundaries (start) are included correctly"
is $(vg msga -g t.vg -s TTCTATAATATG -N | vg stats -s - | wc -l) 1 "soft clips at node boundaries (end) are included correctly"
rm t.vg

is $(vg msga -f msgas/s.fa -b s1 -B 20 -N | vg view -j - | jq -M -c --sort-keys '{"node": .node, "edge": .edge}') $(vg msga -f msgas/s.fa -f msgas/s-rev.fa -b s1 -B 20 -N | vg view -j - | jq -M -c --sort-keys '{"node": .node, "edge": .edge}') "adding in existing sequences in reverse doesn't change graph"

is $((for seq in $(vg msga -f msgas/w.fa -b x -K 16 | vg paths -x - | vg view -a - | jq .sequence | sed s/\"//g ); do grep $seq msgas/w.fa ; done) | wc -l) 2 "the paths of the graph encode the original sequences used to build it"

vg msga -f msgas/w.fa -b x -K 16 -B 20 | vg validate -
is $? 0 "even when banding the paths of the graph encode the original sequences used to build it"

vg msga -f GRCh38_alts/FASTA/HLA/K-3138.fa -B 256 -K 11 -X 1 -E 4 -Q 22 | vg validate -
is $? 0 "HLA K-3138 correctly includes all input paths"

vg msga -f msgas/cycle.fa -b s1 -B 16 -t 1 | vg validate -
is $? 0 "a difficult cyclic path can be included to produce a valid graph"

is $(vg msga -f msgas/cycle.fa -b s1 -B 16 -t 1 | vg mod -D - | vg mod -U 10 - | vg view - | grep ^S | wc -l) 2 "a cyclic path can be normalized"

vg msga -f msgas/l.fa -b a1 -B 8 -m 8 | vg validate -
is $? 0 "edges in cycles with two nodes are correctly included"

vg msga -f GRCh38_alts/FASTA/HLA/B-3106.fa -B 64 -K 11 -X 1 | vg validate -
is $? 0 "HLA B-3106 is assembled into a valid graph"

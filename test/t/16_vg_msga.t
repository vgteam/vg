#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 16

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 | vg mod -n - | vg mod -n - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) e3a28b1662bf849060abc050cc7a4b66 "MSGA produces the expected graph for GRCh38 HLA-V"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 4 | vg mod -n - | vg mod -n - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) e3a28b1662bf849060abc050cc7a4b66 "graph for GRCh38 HLA-V is unaffected by the number of alignment threads"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 10000 -N | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 300 -N | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ )  "varying alignment bandwidth does not affect output graph"

is $(vg msga -f msgas/s.fa -b s1 -B 20 -N | vg view - | md5sum | cut -f 1 -d\ ) 456ec5be516a364ea8238241d80ebd1f "banded alignment can detect and include large deletions in the graph"

is $(vg msga -f msgas/s.fa -b s2 -B 20 -N | vg view - | md5sum | cut -f 1 -d\ ) 456ec5be516a364ea8238241d80ebd1f "assembly around indels works regardless of the base sequence that is used"

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

vg msga -f msgas/cycle.fa -b s1 -B 20 -t 1 | vg validate -
is $? 0 "a difficult cyclic path can be included to produce a valid graph"

is $(vg msga -f msgas/cycle.fa -b s1 -B 20 -t 1 -N | vg view - | wc -l) 35 "a cyclic path can be normalized"

is $(vg msga -f msgas/cycle.fa -b s1 -B 20 -t 1 -P 0.9 | vg mod -n - | vg mod -n - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) \
   $(vg msga -f msgas/cycle.fa -b s1 -B 16 -t 1 -P 0.9 | vg mod -n - | vg mod -n - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) \
   "bandwidth does not affect a cyclic graph, meaning normalization is correct"

vg msga -f msgas/l.fa -b a1 -B 8 -m 8 | vg validate -
is $? 0 "edges in cycles with two nodes are correctly included"

vg msga -f GRCh38_alts/FASTA/HLA/B-3106.fa -B 64 -K 11 -X 1 | vg validate -
is $? 0 "HLA B-3106 is assembled into a valid graph"

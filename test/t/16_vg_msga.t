#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 12

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 4 -k 16 | vg mod -U 10 - | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 | vg mod -U 10 - | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) "graph for GRCh38 HLA-V is unaffected by the number of alignment threads"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 | vg mod -U 10 - | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 16e56f0090b310d2b1479d49cf790324 "MSGA produces the expected graph for GRCh38 HLA-V"

#is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -w 512 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -w 1024 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ )  "varying alignment bandwidth does not affect output graph"
## TODO this test is bad??
is $(vg msga -f msgas/s.fa -w 32 | vg mod -U 10 - | vg mod -c - | vg view - | md5sum | cut -f 1 -d\ ) 603a3c92233e27fc018aa7d3b399b5ac "msga alignment can detect and include large deletions in the graph"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
is $(vg msga -g t.vg -s CAAATTTTCTGGAGTTCTAT -N | vg stats -s - | wc -l) 1 "soft clips at node boundaries (start) are included correctly"
is $(vg msga -g t.vg -s TTCTATAATATG -N | vg stats -s - | wc -l) 1 "soft clips at node boundaries (end) are included correctly"
rm t.vg

vg msga -f msgas/s.fa -b s1 -w 20 | vg mod -U 10 - | vg mod -c - >s.vg
vg msga -g s.vg -f msgas/s-rev.fa -w 20 | vg mod -U 10 - | vg mod -c - >s+rev.vg
is $(vg view s.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg view s+rev.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) "adding in existing sequences in reverse doesn't change graph"
rm -f s.vg s+rev.vg

is $((for seq in $(vg msga -f msgas/w.fa -b x -K 16 | vg paths -x - | vg view -a - | jq .sequence | sed s/\"//g ); do grep $seq msgas/w.fa ; done) | wc -l) 2 "the paths of the graph encode the original sequences used to build it"

vg msga -f msgas/w.fa -b x -K 16 -w 20 | vg validate -
is $? 0 "even when banding the paths of the graph encode the original sequences used to build it"

vg msga -f GRCh38_alts/FASTA/HLA/K-3138.fa -w 256 -W 64 -E 4 | vg validate -
is $? 0 "HLA K-3138 correctly includes all input paths"

vg msga -f msgas/cycle.fa -b s1 -w 32 -e 4 -P 0.95 -t 1 | vg validate -
is $? 0 "a difficult cyclic path can be included to produce a valid graph"

vg msga -f msgas/l.fa -b a1 -w 16 | vg validate -
is $? 0 "edges in cycles with two nodes are correctly included"

vg msga -f GRCh38_alts/FASTA/HLA/B-3106.fa -w 256 -E 4 -B 4 -W 64 -P 0.9 | vg validate -
is $? 0 "HLA B-3106 is assembled into a valid graph"

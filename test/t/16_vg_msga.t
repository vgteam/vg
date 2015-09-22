#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 5

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 | md5sum | cut -f 1 -d\ ) a1bd86103690ffa0d8dada5af1566a49 "MSGA produces the expected graph for GRCh38 HLA-V"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -t 4 | md5sum | cut -f 1 -d\ ) a1bd86103690ffa0d8dada5af1566a49 "graph for GRCh38 HLA-V is unaffected by the number of alignment threads"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 10000 | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 300 | md5sum | cut -f 1 -d\ ) "varying alignment bandwidth does not affect output graph"

is $(vg msga -f msgas/s.fa -k 16 -b s1 -B 20 | md5sum | cut -f 1 -d\ ) 210d2773a036ad3426c42470b1e1e048 "banded alignment can detect and include large deletions in the graph"

is $(vg msga -f msgas/s.fa -k 16 -b s2 -B 20 | md5sum | cut -f 1 -d\ ) b20a61094c66a634a2914688b670a2b0 "assembly around indels works regardless of the base sequence that is used"

#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 2

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 | md5sum | awk '{print $1}') c0494f1307d9ab2606d75fa3e05b0909 "MSGA produces the expected graph for GRCh38 HLA-V"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 10000 | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 300 | md5sum | cut -f 1 -d\ ) "varying alignment bandwidth does not affect output graph"

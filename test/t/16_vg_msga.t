#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 1

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa | md5sum | awk '{print $1}') e0ddd2aa4c0586f1584b5ed0cf5d0cd4 "MSGA produces the expected graph for GRCh38 HLA-V"

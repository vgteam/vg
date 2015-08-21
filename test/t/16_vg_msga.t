#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 1

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa | md5sum | awk '{print $1}') 4aced985f916e96cd9279ff9b1d90d97 "MSGA produces the expected graph for GRCh38 HLA-V"

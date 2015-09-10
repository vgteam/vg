#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 1

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 | md5sum | awk '{print $1}') 1d28f1d2003ded0b3983dd989fdd5fe3 "MSGA produces the expected graph for GRCh38 HLA-V"

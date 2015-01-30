#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 1

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg sim -l 100 -n 100 - | wc -l) 100 \
    "vg sim creates the correct number of reads"


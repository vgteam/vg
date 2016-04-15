#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 3

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg

is $(vg sim -l 100 -n 100 -x x.xg | wc -l) 100 \
    "vg sim creates the correct number of reads"

is $(vg sim -l 100 -n 100 -a -x x.xg | vg view -a - | wc -l) 100 \
   "alignments may be generated rather than read sequences"

is $(vg sim -s 33232 -l 100 -n 100 -J -x x.xg | grep is_reverse | wc -l) 50 "alignments are produced on both strands"

rm -f x.vg x.xg

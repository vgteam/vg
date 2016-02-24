#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 3

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg sim -l 100 -n 100 - | wc -l) 100 \
    "vg sim creates the correct number of reads"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg sim -l 100 -n 100 -a - | vg view -a - | wc -l) 100 \
   "alignments may be generated rather than read sequences"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg sim -s 33232 -l 100 -n 100 -a - | vg view -a - | grep is_reverse | wc -l) 46 "alignments are produced on both strands"

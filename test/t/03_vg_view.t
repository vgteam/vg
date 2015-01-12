#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 2

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -d - | wc -l) 504 "view produces the expected number of lines of dot output"
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - | wc -l) 502 "view produces the expected number of lines of GFA output"


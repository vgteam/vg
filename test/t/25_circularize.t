#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg circularize -p x - | vg view -j - | grep is_circular | wc -l) 1 "a path may be circularized"

#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg circularize -p x - > circular.vg

is $(vg view -j circular.vg | jq -c '.path[] | select(.is_circular)' | wc -l) 1 "a path may be circularized"

vg index -x circular.xg circular.vg
vg convert circular.xg -v  > extracted.vg

is $(vg view -j extracted.vg | jq -c '.path[] | select(.is_circular)' | wc -l) 1 "a circular path survives a round trip to/from xg"

rm -f circular.vg circular.xg extracted.vg


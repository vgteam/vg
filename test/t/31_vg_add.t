#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 5

vg construct -r add/ref.fa > ref.vg
vg add -v add/benedict.vcf ref.vg > benedict.vg
is "$?" "0" "vg add can create a graph"

vg add -v add/rename.vcf -n chrR=ref ref.vg > benedict2.vg
is "$?" "0" "vg add can create a graph with contig renames"

diff benedict.vg benedict2.vg
is "$?" "0" "vg add produces the same graph from VCFs with different contig names"

vg construct -r small/x.fa > x-ref.vg
vg add -v small/x.vcf.gz x-ref.vg > x.vg
is "$?" "0" "vg add can create a slightly larger graph"

is "$(vg view -c x.vg | jq -c '.path[].mapping[] | select(.rank | not)' | wc -l)" "0" "ranks are calculated for emitted paths"

rm -rf ref.vg benedict.vg benedict.vg x-ref.vg x.vg

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 8

vg construct -r add/ref.fa > ref.vg
vg add -v add/benedict.vcf ref.vg > benedict.vg
is "$?" "0" "vg add can create a graph"

vg add -v add/rename.vcf -n chrR=ref ref.vg > benedict2.vg
is "$?" "0" "vg add can create a graph with contig renames"

diff benedict.vg benedict2.vg
is "$?" "0" "vg add produces the same graph from VCFs with different contig names"

vg add -v add/benedict.vcf ref.vg > no-n.vg
vg construct -r add/refN.fa > refN.vg
vg add -v add/benedict.vcf refN.vg > with-n.vg

is "$(vg view -j with-n.vg | jq '.node[].id' | wc -l)" "$(vg view -j no-n.vg | jq '.node[].id' | wc -l)" "having reference Ns does not affect the graph topology"

vg construct -r small/x.fa > x-ref.vg
vg add -v small/x.vcf.gz x-ref.vg > x.vg
is "$?" "0" "vg add can create a slightly larger graph"

is "$(vg view -c x.vg | jq -c '.path[].mapping[] | select(.rank | not)' | wc -l)" "0" "ranks are calculated for emitted paths"

is "$(vg view -Jv add/backward.json | vg add -v add/benedict.vcf - | vg stats -N -)" "5" "graphs with backward nodes can be added to"

is "$(vg view -Jv add/backward_and_forward.json | vg add -v add/benedict.vcf - | vg mod --unchop - | vg stats -N -)" "5" "graphs with backward and forward nodes can be added to"

rm -rf ref.vg benedict.vg benedict.vg x-ref.vg x.vg refN.vg no-n.vg with-n.vg


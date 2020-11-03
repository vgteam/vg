#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 12

vg construct -r add/ref.fa > ref.vg
vg add -v add/benedict.vcf ref.vg > benedict.vg
is "$?" "0" "vg add can create a graph"

vg add -v add/rename.vcf -n chrR=ref ref.vg > benedict2.vg
is "$?" "0" "vg add can create a graph with contig renames"

diff benedict.vg benedict2.vg
is "$?" "0" "vg add produces the same graph from VCFs with different contig names"

vg convert -p ref.vg >ref.pg
vg add -v add/benedict.vcf ref.pg > benedict3.vg
is "$?" "0" "vg add can create a graph from a PackedGraph"

diff benedict.vg benedict3.vg
is "$?" "0" "vg add produces the same graph from the same input in different formats"

vg add -v add/separated.vcf ref.vg > no-n.vg
vg construct -r add/refN.fa > refN.vg
vg add -v add/separated.vcf refN.vg > with-n.vg

is "$(vg view -j with-n.vg | jq '.node[].id' | wc -l)" "$(vg view -j no-n.vg | jq '.node[].id' | wc -l)" "having reference Ns does not affect the graph topology"

vg construct -r add/ngap.fa > ngap.vg
vg add -v add/ngap-offset.vcf ngap.vg > ngap-add.vg
(( EXPECTED_BASES = "$(cat add/ngap.fa | grep -v '>' | tr -d '\n' | wc -c)" + "$(cat add/ngap-offset.vcf | grep -v '#' | wc -l)" ))

is "$(vg stats -l ngap-add.vg | cut -f2)" "${EXPECTED_BASES}" "adding variants adds only the alt bases near large N gaps" 

vg construct -r small/x.fa > x-ref.vg
vg add -v small/x.vcf.gz x-ref.vg > x.vg
is "$?" "0" "vg add can create a slightly larger graph"

is "$(vg view -c x.vg | jq -c '.path[].mapping[] | select(.rank | not)' | wc -l)" "0" "ranks are calculated for emitted paths"

is "$(vg view -Jv add/backward.json | vg add -v add/benedict.vcf - | vg mod --unchop - | vg stats -N -)" "5" "graphs with backward nodes can be added to"

is "$(vg view -Jv add/multi.json | vg add -v add/multi.vcf - | vg mod --unchop - | vg stats -N -)" "5" "graphs with multiple overlapping paths nodes can be added to"

is "$(vg view -Jv add/backward_and_forward.json | vg add -v add/benedict.vcf - | vg mod --unchop - | vg stats -N -)" "5" "graphs with backward and forward nodes can be added to"

rm -rf ref.vg ref.pg benedict.vg benedict2.vg benedict3.vg x-ref.vg x.vg refN.vg no-n.vg with-n.vg ngap.vg ngap-add.vg


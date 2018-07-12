#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg

vg mod -N t.vg >t.ref.vg
vg index -x t.xg t.vg
vg index -x t.ref.xg t.ref.vg

is "$(vg sim -s 7331 -n 10 -l 50 -x t.xg -a | vg annotate -n -x t.ref.xg -a - | awk '{ if ($5 < 50) print }' | wc -l)" "10" "we can detect when reads contain non-reference variation"

rm -f t.vg t.ref.vg t.xg t.ref.xg

vg view -Jv cyclic/circular_path.json > circular_path.vg
vg index -x circular_path.xg circular_path.vg
vg annotate -p -x circular_path.xg -b cyclic/circular_path_origin.bed > circular_path_origin.gam

is "$(vg view -aj circular_path_origin.gam | jq -c '.path.mapping[].position')" "$(printf '{"node_id":"1","offset":"5"}\n{"node_id":"1"}')" "annotations derived from BED can span the end/start joins of circular paths"

rm -f circular_path.vg circular_path.xg circular_path_origin.gam

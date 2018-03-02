#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

# split xy2 graph into an x graph and a y graph, then merge their gbwts
vg construct -m 32 -a -r small/xy.fa -v small/xy2.vcf.gz -R x -C >xy2_x.vg
vg construct -m 32 -a -r small/xy.fa -v small/xy2.vcf.gz -R y -C >xy2_y.vg
vg ids -j xy2_x.vg xy2_y.vg
vg index -x xy2_x.xg -v small/xy2.vcf.gz --gbwt-name xy2_x.gbwt xy2_x.vg
vg index -x xy2_y.xg -v small/xy2.vcf.gz --gbwt-name xy2_y.gbwt xy2_y.vg
is $(vg gbwt --merge xy2_x.gbwt xy2_y.gbwt --fast --output xy2.gbwt -p | grep Sequences | tail -1 | awk '{ print $2 }' | grep 8 | wc -l) 1  "merged gbwt has correct number of sequences"

rm -f xy2_x.vg xy2_x.xg xy2_x.gbwt xy2_y.vg xy2_y.xg xy2_y.gbwt

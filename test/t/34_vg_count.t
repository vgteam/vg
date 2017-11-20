#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

vg construct -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -s 420 -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
vg count -x flat.xg -o 2snp.gam.counts -g 2snp.gam
is $(vg count -x flat.xg -di 2snp.gam.counts  | cut -f 3 | grep -v ^0$ | wc -l) 2 "allele observation counting detects 2 SNPs"

vg pileup -p flat.vg 2snp.gam >2snp.gam.vgpu
is $(vg view -l 2snp.gam.vgpu|  jq '.node_pileups[].base_pileup[].num_bases' | awk '{ print NR-1, $0 }' | head | md5sum | cut -f 1 -d\ )\
   $(vg count -x flat.xg -di 2snp.gam.counts | awk '{ print $1, $2 }' | head | md5sum | cut -f 1 -d\ ) "pileup counts agree with graph coverage"

rm -f flat.vg 2snp.vg 2snp.xg 2snp.sim flat.gcsa flat.gcsa.lcp flat.xg 2snp.xg 2snp.gam 2snp.gam.counts 2snp.gam.vgpu

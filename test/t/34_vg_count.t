#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg construct -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -s 420 -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
vg count -x flat.xg -o 2snp.gam.counts -g 2snp.gam
is $(vg count -x flat.xg -di 2snp.gam.counts  | cut -f 3 | grep -v ^0$ | wc -l) 2 "allele observation counting detects 2 SNPs"

rm -f flat.vg 2snp.vg 2snp.xg 2snp.sim flat.gcsa flat.gcsa.lcp flat.xg 2snp.xg 2snp.gam 2snp.gam.counts

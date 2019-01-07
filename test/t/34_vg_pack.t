#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 6

vg construct -m 1000 -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam -e
is $(vg pack -x flat.xg -di 2snp.gam.cx -e | tail -n+2 | cut -f 5 | grep -v ^0$ | wc -l) 2 "allele observation packing detects 2 SNPs"

vg augment -a pileup flat.vg 2snp.gam -P 2snp.gam.vgpu >/dev/null
is $(vg view -l 2snp.gam.vgpu|  jq '.node_pileups[].base_pileup[] | (.num_bases // 0)' | awk '{ print NR-1, $0 }' | head | md5sum | cut -f 1 -d\ )\
   $(vg pack -x flat.xg -di 2snp.gam.cx -e | awk '{ print $3, $4 }' | tail -n+2 | head | md5sum | cut -f 1 -d\ ) "pileup packs agree with graph coverage"
   
vg pack -x flat.xg -i 2snp.gam.cx -i 2snp.gam.cx -i 2snp.gam.cx -o 2snp.gam.cx.3x

is $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx.3x | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")"/3 | bc) $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")" | bc) "graph coverages are merged from multiple .cx indexes"

is $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx.3x | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")"/3 | bc) $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")" | bc) "edit records are merged from multiple .cx indexes"

x=$(vg pack -x flat.xg -di 2snp.gam.cx | wc -c )
vg pack -x flat.xg -o 2snp.gam.cx -b 10 -g 2snp.gam
y=$(vg pack -x flat.xg -di 2snp.gam.cx | wc -c )
is $x $y "binned edit accumulation does not affect the result"

vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx.3x -i 2snp.gam.cx -i 2snp.gam.cx -i 2snp.gam.cx
x=$(vg pack -x flat.xg -di 2snp.gam.cx.3x | wc -c)
cat 2snp.gam 2snp.gam 2snp.gam | vg pack -x flat.xg -o 2snp.gam.cx -g -
y=$(vg pack -x flat.xg -di 2snp.gam.cx.3x | wc -c)

is $x $y "pack index merging produces the expected result"

rm -f flat.vg 2snp.vg 2snp.xg 2snp.sim flat.gcsa flat.gcsa.lcp flat.xg 2snp.xg 2snp.gam 2snp.gam.cx 2snp.gam.cx.3x 2snp.gam.vgpu

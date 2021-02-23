#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 6

vg construct -m 10 -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg sim -l 30 -x 2snp.vg -n 30 -a -s 1 >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam
# total read bases (30 * 30) / total graph bases 50 = 18
is $(vg depth flat.vg -g 2snp.gam | awk '{print $1}') 18 "vg depth gets correct depth from gam"
is $(vg depth flat.xg -k 2snp.gam.cx -b 100000 | awk '{print int($4)}') 18 "vg depth gets correct depth from pack"
is $(vg depth flat.xg -k 2snp.gam.cx -b 10 | wc -l) 5 "vg depth gets correct number of bins"
vg convert flat.vg -G 2snp.gam | gzip > 2snp.gaf.gz
is $(vg depth flat.vg -a 2snp.gaf.gz | awk '{print $1}') 18 "vg depth gets correct depth from gaf"
vg augment flat.vg 2snp.gam -i > flat-aug.vg
is $(vg depth flat-aug.vg | awk '{print $1}' | uniq | wc -l) $(vg paths -Lv flat-aug.vg | wc -l) "vg depth of paths reports all paths"
is $(vg depth flat-aug.vg -P x | awk '{print $1}' | uniq | wc -l) 1 "vg depth of paths reports just path with selected prefix"
rm -f flat.vg flat.gcsa flat.xg 2snp.vg 2snp.sim 2snp.gam 2snp.gam.cx 2snp.gaf.gz flat-aug.vg

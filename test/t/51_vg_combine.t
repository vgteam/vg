#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10

vg construct -r small/x.fa -v small/x.vcf.gz  >  x.vg
vg construct -r small/x.fa -v small/x.vcf.gz | vg paths -d -v - > y.vg
vg construct -r small/x.fa -v small/x.vcf.gz | vg paths -d -v - > z.vg

vg combine x.vg y.vg z.vg > xyz.vg

is $(vg stats -z xyz.vg | grep nodes | awk '{print $2}') 645 "combined graph as correct total nodes"
is $(vg stats -z xyz.vg | grep edges | awk '{print $2}') 888 "combined graph as correct total edges"
vg paths -Ev x.vg > x.paths
vg paths -Ev xyz.vg > xyz.paths
diff x.paths xyz.paths
is $? 0  "combined graph has same input path"

rm -f x.vg y.vg z.vg xyz.vg x.paths xyz.paths

vg construct -r small/x.fa -v small/x.vcf.gz | vg convert -p - >  x.vg
vg construct -r small/x.fa -v small/x.vcf.gz | vg paths -d -v - | vg convert -a - > y.vg
vg construct -r small/x.fa -v small/x.vcf.gz | vg paths -d -v - | vg convert -p - > z.vg

vg combine x.vg y.vg z.vg > xyz.vg

is $(vg stats -z xyz.vg | grep nodes | awk '{print $2}') 645 "combined handle graph as correct total nodes"
is $(vg stats -z xyz.vg | grep edges | awk '{print $2}') 888 "combined handle graph as correct total edges"
vg paths -Ev x.vg > x.paths
vg paths -Ev xyz.vg > xyz.paths
diff x.paths xyz.paths
is $? 0  "combined handle graph has same input path"

rm -f x.vg y.vg z.vg xyz.vg x.paths xyz.paths

vg construct -r small/x.fa -v small/x.vcf.gz | vg convert -p - >  x.vg
vg construct -r small/x.fa -v small/x.vcf.gz | vg convert -a - > y.vg
vg construct -r small/x.fa -v small/x.vcf.gz | vg convert -v - > z.vg

vg combine -p x.vg y.vg z.vg > xyz.vg

is $(vg stats -z xyz.vg | grep nodes | awk '{print $2}') 645 "path combined graph as correct total nodes"
is $(vg stats -z xyz.vg | grep edges | awk '{print $2}') 890 "path combined path graph as correct total edges"
is $(vg paths -Ev xyz.vg | wc -l) 1 "path combined graph as correct number of paths"
is $(vg paths -Ev xyz.vg | awk '{print $2}') 3003 "combined path has correct size"

rm -f x.vg y.vg z.vg xyz.vg



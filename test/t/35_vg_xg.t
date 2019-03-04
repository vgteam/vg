#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg
vg xg -i x.xg -X y.vg

is $? 0 "converter successfully takes in xg file and outputs vg file"

vg mod -E x.vg | vg view - | grep -v P | sort > x.gfa
vg mod -E y.vg | vg view - | grep -v P | sort > y.gfa
diff x.gfa y.gfa

is $? 0 "files are the same"

rm -f x.xg x.vg y.vg x.gfa y.gfa

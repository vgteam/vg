#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 3

is "$(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg paths --list -v -)" "x" "path listing works from a vg"

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg

is "$(vg paths --list -x x.xg)" "x" "path listing works from an xg"
is "$(vg paths --threads --list -x x.xg -g x.gbwt | wc -l)" "2" "thread listing works from a gbwt"

rm -f x.xg x.gbwt x.vg

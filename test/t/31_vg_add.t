#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg construct -r add/ref.fa > ref.vg
vg add -v add/benedict.vcf ref.vg > benedict.vg
is "$?" "0" "vg add can create a graph"

rm -rf ref.vg benedict.vg

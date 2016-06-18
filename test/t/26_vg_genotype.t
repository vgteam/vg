#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 2

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.vg.xg tiny.vg
vg sim -s 1337 -n 100 -x tiny.vg.xg -l 30 > reads.txt
vg map -r reads.txt -k 8 -V tiny.vg > tiny.gam
vg index -d tiny.gam.index -N tiny.gam

vg genotype tiny.vg tiny.gam.index > /dev/null
is "$?" "0" "vg genotype runs successfully"

vg genotype tiny.vg tiny.gam.index -v > /dev/null
is "$?" "0" "vg genotype runs successfully when emitting vcf"

rm -Rf tiny.vg tiny.vg.xg tiny.gam.index tiny.gam reads.txt




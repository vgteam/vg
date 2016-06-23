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

vg construct -r tiny/tiny.fa >flat.vg
vg index -x flat.xg -g flat.gcsa -k 8 flat.vg
(echo '>flat1'; vg sim -n 1 -l 50 -e 0.05 -s 69 -x flat.xg ) >flat1.fa
vg sim -n 1 -l 50 -e 0.05 -s 69 -x flat.xg -a >flat1.fa.gam
vg construct -r flat1.fa >flat1.vg
vg index -x flat1.xg flat1.vg
vg sim -n 30 -l 50 -e 0.005 -s 7372 -x flat1.xg -a >flat1.sim
(echo '>flat2'; vg sim -n 1 -l 50 -e 0.05 -s 77 -x flat.xg ) >flat2.fa
vg sim -n 1 -l 50 -e 0.05 -s 77 -x flat.xg -a >flat2.fa.gam
vg construct -r flat2.fa >flat2.vg
vg index -x flat2.xg flat2.vg
vg sim -n 30 -l 50 -e 0.005 -s 8675309 -x flat2.xg -a >flat2.sim
vg map -x flat.xg -g flat.gcsa -G <(cat flat1.sim flat2.sim) >flat.gam
vg index -d flat.gam.index -N flat.gam
vg genotype flat.vg flat.gam.index -v
vg mod -i flat.gam flat.vg | vg view -dp - | dot2browser

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 6

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

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.vg.xg tiny.vg
# Simulate 0 reads
vg sim -s 1337 -n 0 -x tiny.vg.xg -l 30 > reads.txt
vg map -r reads.txt -k 8 -V tiny.vg > tiny.gam
vg index -d tiny.gam.index -N tiny.gam

is "$(vg genotype tiny.vg tiny.gam.index -Sp --ref notARealPath 2>&1 | grep 'Found 0 superbubbles' | wc -l)" "1" "vg genotype finds no superbubbles for an empty subset"

is "$(vg genotype tiny.vg tiny.gam.index -vSp 2>&1 | grep 'Found 9 superbubbles' | wc -l)" "1" "vg genotype finds few superbubbles for a subset of just the reference"

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
vg genotype flat.vg flat.gam.index -C >flat.loci
cat tiny/tiny.fa flat1.fa flat2.fa >flats.fa
vg msga -f flats.fa -b x | vg mod -D - | vg mod -n - | vg mod -c - >flat_msga.vg
vg mod -Q flat.loci flat.vg | vg mod -D - | vg mod -n - | vg mod -c - >flat_mod.vg

is $(vg view flat_mod.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) \
   $(vg view flat_msga.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) \
   "called genotypes are correct in a small simulated example"

vg view -q flat.loci | vg view -qJz - | vg view -q - >/dev/null
is "$?" "0" "genotype format can be converted to and from JSON"

rm -rf flat.vg flat.xg flat.gcsa.lcp flat.gcsa flat1.fa flat1.fa.gam flat1.vg flat1.xg flat1.sim flat2.fa flat2.fa.gam flat2.vg flat2.xg flat2.sim flat.gam flat.gam.index flats.fa flats.fa.fai flat.loci flat_msga.vg flat_mod.vg

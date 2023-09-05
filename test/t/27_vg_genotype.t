#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 5

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.vg.xg -g tiny.vg.gcsa -k 16 tiny.vg
vg sim -n 100 -x tiny.vg.xg -l 30 > reads.txt
vg map -T reads.txt -g tiny.vg.gcsa -x tiny.vg.xg > tiny.gam

vg genotype tiny.vg tiny.gam > /dev/null
is "$?" "0" "vg genotype runs successfully"

vg genotype tiny.vg tiny.gam -v > /dev/null
is "$?" "0" "vg genotype runs successfully when emitting vcf"

rm -Rf tiny.vg tiny.vg.xg tiny.gam reads.txt

vg construct -r tiny/tiny.fa >flat.vg
vg index -x flat.xg -g flat.gcsa -k 8 flat.vg
cat tiny/flat-s69-n1-l50-e0.05.fa >flat1.fa
cat tiny/flat-s69-n1-l50-e0.05.gam >flat1.fa.gam
vg construct -r flat1.fa >flat1.vg
vg index -x flat1.xg flat1.vg
cat tiny/flat1-s7372-n30-l50-e0.005.gam >flat1.sim
cat tiny/flat-s77-n1-l50-e0.05.fa >flat2.fa
cat tiny/flat-s77-n1-l50-e0.05.gam >flat2.fa.gam
vg construct -r flat2.fa >flat2.vg
vg index -x flat2.xg flat2.vg
cat tiny/flat2-s8675309-n30-l50-e0.005.gam >flat2.sim
vg map -x flat.xg -g flat.gcsa -G <(cat flat1.sim flat2.sim) >flat.gam
vg genotype flat.vg flat.gam -t 1 >flat.loci
cat tiny/tiny.fa flat1.fa flat2.fa >flats.fa
vg msga -f flats.fa -b x | vg paths -d -v - | vg mod -n - | vg mod -c - >flat_msga.vg
vg augment flat.vg -L flat.loci | vg paths -d -v - | vg mod -n - | vg mod -c - >flat_mod.vg

is $(vg view flat_mod.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) \
   $(vg view flat_msga.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) \
   "called genotypes are correct in a small simulated example"

vg view -q flat.loci | vg view -qJz - | vg view -q - >/dev/null
is "$?" "0" "genotype format can be converted to and from JSON"

rm -rf flat.vg flat.xg flat.gcsa.lcp flat.gcsa flat1.fa flat1.fa.gam flat1.vg flat1.xg flat1.sim flat2.fa flat2.fa.gam flat2.vg flat2.xg flat2.sim flat.gam flats.fa flats.fa.fai flat.loci flat_msga.vg flat_mod.vg

vg construct -v call/bigins.vcf.gz -r tiny/tiny.fa > bigins.vg
vg index -x bigins.vg.xg -g bigins.vg.gcsa -k 16 bigins.vg
cat call/bigins-s1337-n100-l12.reads > reads.txt
vg map -T reads.txt -g bigins.vg.gcsa -x bigins.vg.xg > bigins.gam
is "$(vg genotype bigins.vg bigins.gam -t 1 -v | grep GACGTTACAATGAGCCCTACAGACATATC | wc -l)" "1" "genotype finds big insert" 

rm -rf bigins.vg bigins.vg.xg bigins.vg.gcsa bigins.vg.gcsa.lcp bigins.gam reads.txt 


#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 5

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.vg.xg tiny.vg
vg sim -a -s 1337 -n 100 -x tiny.vg.xg -l 30 > reads.gam
vg map -G reads.gam -k 8 -V tiny.vg > tiny.gam
vg index -d tiny.gam.index -N tiny.gam
vg genotype tiny.vg tiny.gam.index >tiny.loci
is $(vg locify -g tiny.gam.index -x tiny.vg.xg -l tiny.loci -f -n -o loci.sorted | vg view -a - | wc -l) 100 "locify produces output for each input alignment"
is $(cat loci.sorted | wc -l) 6 "the sorted list of loci is the right length"
is $(head -1 loci.sorted) "1+0_1+0" "the first locus is as expected"
is $(head -4 loci.sorted | tail -1) "9+18_12+0" "a middle locus is as expected"
is $(tail -1 loci.sorted) "15+9_15+9" "the last locus is as expected"

rm -rf tiny.vg tiny.vg.xg tiny.vg reads.gam tiny.gam tiny.gam.index tiny.loci loci.sorted

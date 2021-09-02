#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 7

vg construct -r small/x.fa -v small/x.vcf.gz >s.vg
vg index -x s.xg -g s.gcsa s.vg
vg snarls -T s.vg > s.snarls
vg index -s s.snarls -j s.dist s.vg
vg sim -n 1000 -l 100 -e 0.01 -i 0.005 -x s.xg -a >s.sim

is $(vg map -x s.xg -g s.gcsa -G s.sim --surject-to sam | vg inject -x s.xg - | vg gamcompare - s.sim | vg view -a - | wc -l) 1000 "gamcompare completes"

is $(vg gamcompare --range 10 s.sim s.sim | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l) 1000 "gamcompare says the truth is correctly mapped"

# Map a couple adjacent reads with multi-positioning
vg map -x s.xg  -g s.gcsa -s "AATCTCTCTGAACTTCAGTTTAATTATC" > read1.gam
vg annotate -a read1.gam -p -x s.xg > read1.single.gam
vg annotate -a read1.gam -m -x s.xg > read1.multi.gam
vg map -x s.xg  -g s.gcsa -s "TCTAATATGGAGATGATACTACTGACAG" > read2.gam
vg annotate -a read2.gam -p -x s.xg > read2.single.gam
vg annotate -a read2.gam -m -x s.xg > read2.multi.gam

is "$(vg gamcompare -r 30 read1.single.gam read2.single.gam | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l)" "1" "Reads annotated with leftmost position only are close enough at a long distance"
is "$(vg gamcompare -r 30 -d s.dist read1.gam read2.gam | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l)" "1" "Reads are close enough at a long distance with a distance index"
is "$(vg gamcompare -r 10 read1.single.gam read2.single.gam | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l)" "0" "Reads annotated with leftmost position only are too far apart at a short distance"
is "$(vg gamcompare -r 10 -d s.dist read1.gam read2.gam | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l)" "0" "Reads are too far apart at a short distance with a distance index"
is "$(vg gamcompare -r 10 read1.multi.gam read2.multi.gam | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l)" "1" "Reads annotated with multiple positions only are close enough at a short distance"

rm -f read1.gam read1.single.gam read1.multi.gam read2.gam read2.single.gam read2.multi.gam 

rm -f s.vg s.xg s.gcsa s.gcsa.lcp s.sim s.snarls s.dist

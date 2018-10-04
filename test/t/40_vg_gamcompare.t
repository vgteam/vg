#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 2

vg construct -r small/x.fa -v small/x.vcf.gz >s.vg
vg index -x s.xg -g s.gcsa s.vg
vg sim -n 1000 -l 100 -e 0.01 -i 0.005 -x s.xg -a >s.sim

is $(vg map -x s.xg -g s.gcsa -G s.sim --surject-to sam | vg inject -x s.xg - | vg gamcompare - s.sim | vg view -a - | wc -l) 1000 "gamcompare completes"

is $(vg gamcompare --range 10 s.sim s.sim | vg view -aj - | jq -c 'select(.correctly_mapped)' | wc -l) 1000 "gamcompare says the truth is correctly mapped"

rm -f s.vg s.xg s.gcsa s.gcsa.lcp s.sim

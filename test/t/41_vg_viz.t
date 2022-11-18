#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 5

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg
vg index -x t.xg -g t.gcsa t.vg
vg map -G <(vg sim -l 30 -n 100 -x t.xg -a) -d t | vg pack -x t.xg -o t.cx -g -

vg viz -x t.xg -o t.svg -i t.cx -n alignments
is $(echo $(wc -c t.svg | cut -f 1 -d\ )' > 0' | bc) 1 "vg viz runs"

rm -f t.png
vg viz -x t.xg -o t.png
is "${?}" "0" "vg viz succeeds when drawing a PNG"
stat t.png 2>/dev/null >/dev/null
is "${?}" "0" "vg viz actually creates an output PNG"

rm -f cactus.png
vg viz -x graphs/cactus-BRCA2.gfa -o cactus.png
is "${?}" "0" "vg viz succeeds when drawing a larger PNG from GFA"
stat cactus.png 2>/dev/null >/dev/null
is "${?}" "0" "vg viz actually creates a larger PNG from GFA"

rm -f t.vg t.xg t.gcsa t.gcsa.lcp t.cx t.svg t.png cactus.png

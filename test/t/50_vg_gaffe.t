#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg
vg snarls --include-trivial x.vg > x.snarls
vg index -s x.snarls -j x.dist x.vg
vg minimizer -k 29 -w 11 -g x.gbwt -i x.min x.xg

vg gaffe -x x.xg -H x.gbwt -m x.min -d x.dist -f reads/small.middle.ref.fq > mapped.gam
is "${?}" "0" "a read can be mapped with all indexes specified without crashing"

rm -f x.vg x.xg x.gbwt x.snarls x.min x.dist mapped.gam

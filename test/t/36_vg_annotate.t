#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg

vg mod -N t.vg >t.ref.vg
vg index -x t.xg t.vg
vg index -x t.ref.xg t.ref.vg

is $(vg sim -s 7331 -n 10 -l 50 -x t.xg -a | vg annotate -n -x t.ref.xg -a - | awk '{ if ($5 < 50) print }' | wc -l) 10 "we can detect when reads contain non-reference variation"

rm -f t.vg t.ref.vg t.xg t.ref.xg

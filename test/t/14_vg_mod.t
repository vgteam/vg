#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 2

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^P | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^S | wc -l) \
    "vg mod yields a graph with only a particular path"

is $(vg mod -o graphs/orphans.vg | vg view - | wc -l) 8 "orphan edge removal works"

vg construct -r tiny/tiny.fa >t.vg
vg index -s -k 11 t.vg
vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg view -a -
vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view -

rm t.vg
rm -rf t.vg.index

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 4

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg

is $(vg snarls -u t.vg | head -1 | cut -f 3) 1,2,3,4,5,6, "a ultrabubble's internal nodes are correctly reported"
rm -f t.vg

is $(cat graphs/missed_bubble.gfa | vg view -Fv - | vg snarls -u - | grep 79,80,81,299, | wc -l) 1 "ultrabubbles are detected even when the graph initially has reversing edges"

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg

vg view -J -v snarls/snarls.json > snarls.vg
is $(vg snarls snarls.vg -r st.pb | vg view -R - | wc -l) 3 "vg snarls made right number of protobuf Snarls"
is $(vg view -E st.pb | wc -l) 6 "vg snarls made right number of protobuf SnarlTraversals"

rm -f snarls.vg st.pb 


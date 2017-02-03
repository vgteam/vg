#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 8

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg

is $(vg snarls -b t.vg | head -1 | cut -f 3) 1,2,3,4,5,6, "a superbubble's internal nodes are correctly reported"
is $(vg snarls -u t.vg | head -1 | cut -f 3) 1,2,3,4,5,6, "a ultrabubble's internal nodes are correctly reported"
rm -f t.vg

is $(cat graphs/missed_bubble.gfa | vg view -Fv - | vg snarls -b - | grep 79,80,81,299, | wc -l) 1 "superbubbles are detected even when the graph initially has reversing edges"
is $(cat graphs/missed_bubble.gfa | vg view -Fv - | vg snarls -u - | grep 79,80,81,299, | wc -l) 1 "ultrabubbles are detected even when the graph initially has reversing edges"

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg

vg snarls z.vg -b > sb.txt
vg snarls z.vg -u > cb.txt
is $(diff sb.txt cb.txt | wc -l) 0 "superbubbles and cactus bubbles identical for 1mb1kgp"

is $(vg snarls z.vg -r st.pb | vg view -R - | wc -l) 27303 "vg snarls made right number of protobuf Snarls"
is $(vg view -E st.pb | wc -l) 57570 "vg snarls made right number of protobug SnarlTraversals"

rm -f z.vg sb.txt cb.txt st.pb 

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg mod -X 1 - > tiny.vg

vg snarls tiny.vg -b > sb.txt
vg snarls tiny.vg -u > cb.txt
is $(diff sb.txt cb.txt | wc -l) 0 "superbubbles and cactus bubbles identical for atomized tiny"

rm sb.txt cb.txt tiny.vg

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -j - | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "alignment traverses the correct path"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -j - | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 96 "alignment score is as expected"

is $(vg align -js $(cat mapsoftclip/70211809-70211845.seq) mapsoftclip/70211809-70211845.vg | jq -c '.path .mapping[0] .position .node_id') 70211814 "alignment does not contain excessive soft clips"

is $(vg align -js $(cat mapsoftclip/113968116:113968146.seq ) mapsoftclip/113968116:113968146.vg | jq -c ".score") 274 "alignment does not overflow when using 8x16bit vectors"

is $(vg align -js $(cat mapsoftclip/280136066-280136088.seq) mapsoftclip/280136066-280136088.vg | jq -c '.path .mapping[0] .position .node_id') 280136076 "Ns do not cause excessive soft clipping"

is $(vg align -js GGCTATGTCTGAACTAGGAGGGTAGAAAGAATATTCATTTTGGTTGCCACAAACCATCGAAACAAAGATGCAGGTCATTGATGTAAAACTACAGTTAGTTCCTACTGACTCCTTTTCAGCTTCTCTTCATTGCTATGAGCCAGCGTCTCCT graphs/59867692-59867698.vg | jq '.path.mapping[0].position.node_id') 59867694 "nodes are only referenced if they have mappings"

vg construct -r tiny/tiny.fa >t.vg
seq=CAAATAAGGCTTGGAAATGTTCTGGAGTTCTATTATATTCCAACTCTCTT
vg align -s $seq t.vg | vg mod -i - t.vg >t2.vg
is $(vg align -s $seq -Q query t2.vg | vg mod -i - -P t2.vg | vg view - | grep query | wc -l) 4 "align can use query names and outputs GAM"
rm t.vg t2.vg


is $(vg align -s TATATATATACCCCCCCCC -j cyclic/all.vg | jq ".path.mapping[].position.node_id" | tr '\n' ',' | grep "5,6" | wc -l) 1  "alignment to cyclic graphs works"

vg align -s ACGT -j cyclic/reverse_self.vg >/dev/null
is $? 0  "graphs where duplicated nodes need flipping can be used for alignment"

vg align -s AGTCCTTGAAAGAGGGCAAAATAAACTGTTAGTAGAGCCAGGTCTGAAAACAACACTTTCTTGC inverting/m.vg >/dev/null
is $? 0 "node flipping doesn't destroy the alignment"

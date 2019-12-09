#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 20

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz > x.vg

is $(vg align x.vg -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --full-l-bonus 0 -j | tr ',' '\n' | grep node_id | grep "72\|73\|76\|77" | wc -l) 4 "alignment traverses the correct path"

is $(vg align x.vg -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --full-l-bonus 0 -j | jq -r '.score') 48 "alignment score is as expected"

is $(vg align x.vg -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --full-l-bonus 5 -j | jq -r '.score') 58 "full length bonus works and is included by default"

is $(vg align  x.vg -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTA --full-l-bonus 5 --pinned --pin-left -j - | jq -r '.score') 37 "bonuses are included on only one end for pinned alignments"

is $(vg align  x.vg --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 --full-l-bonus 0 -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -j - | jq -r '.score') 96 "scoring parameters are respected"

is $(vg align  x.vg --score-matrix 2_2.mat --gap-open 3 --gap-extend 1 --full-l-bonus 0 -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -j - | jq -r '.score') 96 "score-matrix file should give same results as --match 2 --mismatch 2: scoring parameters are respected"

rm -f x.vg

is $(vg align -js $(cat mapsoftclip/70211809-70211845.seq) --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 --full-l-bonus 0 mapsoftclip/70211809-70211845.vg | jq -r -c '.path .mapping[0] .position .node_id') 70211814 "alignment does not contain excessive soft clips under lenient scoring"

is $(vg align --score-matrix 2_2.mat -js $(cat mapsoftclip/70211809-70211845.seq) --gap-open 3 --gap-extend 1 --full-l-bonus 0 mapsoftclip/70211809-70211845.vg | jq -r -c '.path .mapping[0] .position .node_id') 70211814 "score-matrix file should give same results as --match 2 --mismatch 2: alignment does not contain excessive soft clips under lenient scoring"

is $(vg align -js $(cat mapsoftclip/113968116:113968146.seq ) --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 --full-l-bonus 0 mapsoftclip/113968116:113968146.vg | jq -r -c ".score") 274 "alignment score does not overflow at 255 when using 8x16bit vectors"

is $(vg align --score-matrix 2_2.mat -js $(cat mapsoftclip/113968116:113968146.seq ) --gap-open 3 --gap-extend 1 --full-l-bonus 0 mapsoftclip/113968116:113968146.vg | jq -r -c ".score") 274 "score-matrix file should give same results as --match 2 --mismatch 2: alignment score does not overflow at 255 when using 8x16bit vectors"

is $(vg align -js $(cat mapsoftclip/280136066-280136088.seq) mapsoftclip/280136066-280136088.vg | jq -r -c '.path .mapping[0] .position .node_id') 280136076 "Ns do not cause excessive soft clipping"

is $(vg align -js GGCTATGTCTGAACTAGGAGGGTAGAAAGAATATTCATTTTGGTTGCCACAAACCATCGAAACAAAGATGCAGGTCATTGATGTAAAACTACAGTTAGTTCCTACTGACTCCTTTTCAGCTTCTCTTCATTGCTATGAGCCAGCGTCTCCT graphs/59867692-59867698.vg | jq -r '.path.mapping[0].position.node_id') 59867694 "nodes are only referenced if they have mappings"

vg construct -m 1000 -r tiny/tiny.fa >t.vg
seq=CAAATAAGGCTTGGAAATGTTCTGGAGTTCTATTATATTCCAACTCTCTT
vg align -s $seq t.vg | vg augment t.vg - -i -S >t2.vg
is $(vg align -s $seq -Q query t2.vg | vg augment t2.vg - -i -B -S | vg view - | grep "query" | cut -f 3 | grep -o "[0-9]\+" | wc -l) 4 "align can use query names and outputs GAM"
rm t.vg t2.vg


is $(vg align -s TATATATATACCCCCCCCC -j cyclic/all.vg | jq -r ".path.mapping[].position.node_id" | tr '\n' ',' | grep "5,6" | wc -l) 1  "alignment to cyclic graphs works"

vg align -s ACGT -j cyclic/reverse_self.vg >/dev/null
is $? 0  "graphs where duplicated nodes need flipping can be used for alignment"

vg align -s AGTCCTTGAAAGAGGGCAAAATAAACTGTTAGTAGAGCCAGGTCTGAAAACAACACTTTCTTGC inverting/m.vg >/dev/null
is $? 0 "node flipping doesn't destroy the alignment"

vg align -s ATTTTTAACTCCATGTTTGAGAAACATTTAATAATGTAATGTGTTTGTGGCACAGCAGGAGTAC graphs/difficult-inv.vg >/dev/null
is $? 0 "alignment correctly handles an inversion"

vg align -s AAACATACATTTTC graphs/exploding.vg >/dev/null
is $? 0 "the exploding graph doesn't blow up"

is $(vg align -s GTAATGGTAATGGATATGTTGGGCTTTTTTCTTT -j -p graphs/f.vg | jq -r '.path | length') 1 "pinning doesn't cause problems for gssw"

is $(vg align -s GTAATGGTAATGGATATGTTGGGCTTTTTTCTTT -j -p -L graphs/f.vg | jq -r '.path | length') 1 "left pinning is correctly handled"

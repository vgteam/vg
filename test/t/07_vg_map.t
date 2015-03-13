#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 4

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s -k 11 x.vg

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg -J | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "global alignment traverses the correct path"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg -J | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 96 "alignment score is as expected"

vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -d x.vg.index >/dev/null
is $? 0 "vg map takes -d as input without a variant graph"

#is $(vg map -s TCAGATTCTCATCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG x.vg | grep to_length | grep '"sequence": "T"' | wc -l) 1 "vg map can align across a SNP"

vg map -s TCAGATTCTCATCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG x.vg

is $(vg map -r <(vg sim -n 69 -n 1000 -l 100 x.vg) x.vg | vg view -a - | jq -c '.score == 200 // [.score, .sequence]' | grep -v true | wc -l) 0 "alignment is perfect on a small graph"

rm x.vg
rm -rf x.vg.index

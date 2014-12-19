#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 2

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG - | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "alignment traverses the correct path"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG - | tr ',' '\n' | grep score | awk '{ print $2 }') 96 "alignment score is as expected"

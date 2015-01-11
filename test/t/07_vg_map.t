#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 2

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s -k 11 x.vg

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "global alignment traverses the correct path"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 96 "alignment score is as expected"

rm x.vg
rm -rf x.vg.index

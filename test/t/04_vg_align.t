#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 5

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG - | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "alignment traverses the correct path"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG - | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 96 "alignment score is as expected"

is $(vg align -s $(cat mapsoftclip/70211809-70211845.seq) mapsoftclip/70211809-70211845.vg | jq -c '.path .mapping[0] .position .node_id') 70211814 "alignment does not contain excessive soft clips"

is $(vg align -s $(cat mapsoftclip/113968116:113968146.seq ) mapsoftclip/113968116:113968146.vg | jq -c ".score") 274 "alignment does not overflow when using 8x16bit vectors"

is $(vg align -s $(cat mapsoftclip/280136066-280136088.seq) mapsoftclip/280136066-280136088.vg | jq -c '.path .mapping[0] .position .node_id') 280136076 "Ns do not cause excessive soft clipping"

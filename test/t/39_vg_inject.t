#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 8

vg construct -r small/x.fa > j.vg
vg index -x j.xg j.vg
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index -k 11 -g x.gcsa -x x.xg x.vg

#samtools view -f 0 small/x.sam | awk '$6 ~ /100M/' > p.bam
#samtools view -f 0 small/x.sam | awk '$6 !~ /100M/' > pm.bam
samtools view -f 16 small/x.sam > i.bam
#samtools view -f 16 small/x.sam | awk '$6 ~ /100M/' > i.bam
#samtools view -f 16 small/x.sam | awk '$6 !~ /100M/' > im.bam

is $(vg inject -p x -x x.xg small/x.sam | vg view -a - | wc -l) \
    1000 "reads are generated"

is $(vg inject -p x -x x.xg small/x.sam | vg surject -p x -x x.xg -t 1 - | vg view -a - | wc -l) \
    1000 "vg inject works for all reads included in the original bam"

is "$(samtools view -S small/x.sam | cut -f 1 | sort)" "$(vg inject -p x -x x.xg small/x.sam | vg surject -p x -x x.xg -t 1 -s - | sort -k 4n | cut -f 1)" \
    "vg inject retains read names"

is "$(samtools view -S small/x.sam | cut -f 4)" "$(vg inject -p x -x x.xg small/x.sam | vg surject -p x -x x.xg -t 1 -s - | cut -f 4 | sort)" \
    "vg inject works perfectly for the position of alignment"

is "$(samtools view -S small/x.sam | cut -f 5)" "$(vg inject -p x -x x.xg small/x.sam | vg surject -p x -x x.xg -t 1 -s - | sort -k 4n | cut -f 5)" \
    "vg inject works perfectly for the mappign quality of alignment"

is "$(samtools view -S small/x.sam | cut -f 6)" "$(vg inject -p x -x x.xg small/x.sam | vg surject -p x -x x.xg -t 1 -s - | sort -k 4n | cut -f 6)" \
    "vg inject works perfectly for the cigar of alignment"

is "$(samtools view -S small/x.sam | cut -f 6)" "$(vg inject -x x.xg small/x.sam | vg surject -p x -x x.xg -t 1 -s - | sort -k 4n | cut -f 6)" \
    "vg inject doesn't need to be told which path to use"

is "$(vg inject -p -x x.xg i.bam | vg view -a - | jq .path.mapping[0].position.is_reverse | grep -v true | wc -l)" \
    0 "vg inject works perfectly for the reads flagged as is_reverse"

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 14

vg construct -r small/x.fa > j.vg
vg index -x j.xg j.vg
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index -k 11 -g x.gcsa -x x.xg x.vg

#samtools view -f 0 small/x.bam | awk '$6 ~ /100M/' > p.bam
#samtools view -f 0 small/x.bam | awk '$6 !~ /100M/' > pm.bam
#samtools view -O bam -f 16 small/x.bam > small/i.bam
#samtools view -f 16 small/x.bam | awk '$6 ~ /100M/' > i.bam
#samtools view -f 16 small/x.bam | awk '$6 !~ /100M/' > im.bam

is $(vg inject -x x.xg small/x.bam | vg view -a - | wc -l) \
    1000 "reads are generated"

is $(vg inject -x x.xg small/x.bam | vg surject -x x.xg -t 1 - | vg view -a - | wc -l) \
    1000 "vg inject works for all reads included in the original bam"

is "$(samtools view small/x.bam | cut -f 1 | sort)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 1 | sort)" \
    "vg inject retains read names"

is "$(samtools view small/x.bam | cut -f 4)" "$(vg inject -x j.xg small/x.bam | vg surject -p x -x j.xg -t 1 -s - | samtools view -S - | cut -f 4 | sort -n)" \
    "vg inject works perfectly for the position of alignment in the simple graph"

is "$(samtools view small/x.bam | cut -f 4)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 4 | sort -n)" \
    "vg inject works perfectly for the position of alignment"

is "$(samtools view small/x.bam | cut -f 5 | grep 60 | wc -l)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 5 | grep 60 | wc -l)" \
    "vg inject works perfectly for the mapping quality of alignment"

is "$(samtools view small/x.bam | sort -k 1 | cut -f 6)" "$(vg inject -x j.xg small/x.bam | vg surject -p x -x j.xg -t 1 -s - | samtools view -S - | sort -k 1 | cut -f 6)" \
    "vg inject works perfectly for the cigar of alignment in the simple graph"

is "$(samtools view small/x.bam | sort -k 1 | cut -f 6)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | sort -k 1 | cut -f 6)" \
    "vg inject works perfectly for the cigar of alignment"

is "$(vg inject -x x.xg small/i.bam | vg view -a - | wc -l)" 470 "vg inject preserves all reads"

is "$(vg inject -x x.xg small/i.bam | vg view -a - | jq .path.mapping[0].position.is_reverse | grep -v true | wc -l)" \
    0 "vg inject works perfectly for the reads flagged as is_reverse"

cat <(samtools view -H small/x.bam) <(printf "name\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGT\tHHHHHHHHHHHH\n") > unmapped.sam
is "$(vg inject -x x.xg unmapped.sam | vg view -aj - | grep "path" | wc -l)" 0 "vg inject does not make an alignment for an umapped read"
is "$(echo $?)" 0 "vg inject does not crash on an unmapped read"

is $(vg inject -x x.xg small/x.bam -o GAF | wc -l) \
    1000 "vg inject supports GAF output"

samtools cat small/x.bam small/i.bam > all.bam
is "$(vg inject -x x.xg all.bam | vg validate -A -c -a - x.xg 2>&1 | grep invalid | wc -l)" 0 "vg inject produces valid alignments of the read and reference"

rm j.vg j.xg x.vg x.gcsa x.gcsa.lcp x.xg unmapped.sam all.bam

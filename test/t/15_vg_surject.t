#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 5

vg construct -r small/x.fa >j.vg
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s -k 11 x.vg

is $(vg map -r <(vg sim -s 1337 -n 100 j.vg) x.vg | vg surject -p x -d x.vg.index - | vg view -a - | jq .score | grep 200 | wc -l) \
    100 "vg surject works perfectly for perfect reads derived from the reference"

is $(vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg surject -p x -d x.vg.index - | vg view -a - | wc -l) \
    100 "vg surject works for every read simulated from a dense graph"

is $(vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg surject -p x -d x.vg.index -s - | grep -v ^@ | wc -l) \
    100 "vg surject produces valid SAM output"

is $(vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg surject -p x -d x.vg.index -b - | samtools view - | wc -l) \
    100 "vg surject produces valid BAM output"

is $(vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg surject -p x -d x.vg.index -c - | samtools view - | wc -l) \
    100 "vg surject produces valid CRAM output"

rm -rf j.vg x.vg x.vg.index

vg index -s -k 27 -e 7 graphs/fail.vg

read=TTCCTGTGTTTATTAGCCATGCCTAGAGTGGGATGCGCCATTGGTCATCTTCTGGCCCCTGTTGTCGGCATGTAACTTAATACCACAACCAGGCATAGGTGAAAGATTGGAGGAAAGATGAGTGACAGCATCAACTTCTCTCACAACCTAG
#is $(vg map -s $read graphs/fail.vg | vg surject -i graphs/GRCh37.path_names -d graphs/fail.vg.index -s - | grep $read | wc -l) 1 "surjection of a cigar with an insertion and a deletion"

rm -rf graphs/fail.vg.index

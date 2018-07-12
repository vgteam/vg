#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

## Check whether depth filtering passes all alignments when the minimum depth is zero.

# Check downsampling

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg 
vg sim -l 100 -n 100 -p 50 -x x.xg -s 1 -a > sim.gam
vg surject -x x.xg -s -i sim.gam > sim.sam

vg filter -d 123.5 -t 10 sim.gam > filtered.gam
samtools view -s 123.5 sim.sam > filtered.sam

is "$(vg view -aj filtered.gam | jq -rc '.name' | sed 's/_[12]//g' | sort | md5sum)" "$(cat filtered.sam | grep -v "^@" | cut -f1 | sort | md5sum)" "samtools and vg filter agree on how to select downsampled paired reads"

vg sim -l 100 -n 100 -x x.xg -s 1 -a > sim.gam
vg surject -x x.xg -s sim.gam > sim.sam

vg filter -d 456.2 -t 10 sim.gam > filtered.gam
samtools view -s 456.2 sim.sam > filtered.sam

is "$(vg view -aj filtered.gam | jq -rc '.name' | sort | md5sum)" "$(cat filtered.sam | grep -v "^@" | cut -f1 | sort | md5sum)" "samtools and vg filter agree on how to select downsampled unpaired reads"   

rm -f x.vg x.xg sim.gam sim.sam filtered.gam filtered.sam


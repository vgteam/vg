#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg  x.vg
vg sim -x x.xg -l 100 -n 5000 -e 0.01 -i 0.001 -a > x.gam

# sanity check: does passing no options preserve input
is $(vg filter x.gam | vg view -a - | jq . | grep mapping | wc -l) 5000 "vg filter with no options preserves input."

# basic chunking tests
printf "x\t2\t8\nx\t8\t20\ny\t0\t1\nx\t150\t500\nx\t0\t100000000\n" > chunks.bed
vg filter -x x.xg -R chunks.bed -B filter_chunk x.gam

# right number of chunks
is $(ls -l filter_chunk-*.gam | wc -l) 5 "vg filter makes right number of chunks."

# is chunk 0 (2-3) comprised of nodes 1,2,4? 
is $(vg view -a filter_chunk-0.gam | jq -c '.path.mapping[].position' | jq 'select ((.node_id == "1") or (.node_id == "2") or (.node_id == "4"))' | grep node | sed s/,// | sort | uniq | wc -l) 3 "vg filter left chunk has all left nodes"

# check that chunk 4 is off to the right a bit
is $(vg view -a filter_chunk-3.gam | jq -c '.path.mapping[].position' | jq 'select (((.node_id | tonumber) < 4))' | wc -l) 0 "vg filter right chunk has no left nodes"

# check that chunk 5 is everything
is $(vg view -a filter_chunk-4.gam | jq . | grep mapping | wc -l) 5000 "vg filter big chunk has everything"

# Downsampling works
SAMPLED_COUNT=$(vg filter x.gam --downsample 0.5 | vg view -a - | jq . | grep mapping | wc -l)
OUT_OF_RANGE=0
if [[ "${SAMPLED_COUNT}" -lt 2000 || "${SAMPLED_COUNT}" -gt 3000 ]]; then
    # Make sure it's in a reasonable range for targeting 50%.
    # We won't get 50% always because it is random sampling.
    # Sometimes it might even be outside this range!
    # But my binomial calculator said the probability of that was NaN so I bet it won't happen
    OUT_OF_RANGE=1
fi

is "${OUT_OF_RANGE}" "0" "vg filter downsamples correctly"

rm -f x.gam filter_chunk*.gam chunks.bed

cp small/x-s1-l100-n100-p50.gam paired.gam
cp small/x-s1-l100-n100.gam single.gam

# Check downsampling against samtools 1.0+
# If the installed samtools isn't new enough, fall back on precalculated hashes
PAIRED_HASH=a31e81e05f86224b5aec73f6f8e2d9a9
SINGLE_HASH=854be75583698d41023bf073d26833d6
if samtools 2>&1 | grep "Version" | grep -v ": 0" >/dev/null; then
    # We can calculate hashes ourselves
    vg surject -x x.xg -s -i paired.gam > paired.sam
    samtools view -s 123.5 -S paired.sam > filtered.sam
    
    PAIRED_HASH="$(cat filtered.sam | grep -v "^@" | cut -f1 | sort | md5sum | cut -f1 -d' ')"
    
    vg surject -x x.xg -s single.gam > single.sam
    samtools view -s .2 -S single.sam > filtered.sam
    
    SINGLE_HASH="$(cat filtered.sam | grep -v "^@" | cut -f1 | sort | md5sum | cut -f1 -d' ')"
fi

vg filter -d 123.5 -t 10 paired.gam > filtered.gam

is "$(vg view -aj filtered.gam | jq -rc '.name' | sed 's/_[12]//g' | sort | md5sum | cut -f1 -d' ')" "${PAIRED_HASH}"  "samtools 1.0+ and vg filter agree on how to select downsampled paired reads"

# "0." and "." syntax need to mean the same thing.
vg filter -d 0.2 -t 10 single.gam > filtered.gam

is "$(vg view -aj filtered.gam | jq -rc '.name' | sort | md5sum | cut -f1 -d' ')" "${SINGLE_HASH}" "samtools 1.0+ and vg filter agree on how to select downsampled unpaired reads"

vg annotate -p -x x.xg -a paired.gam > paired.annotated.gam
is "$(vg filter -X "[a-f]" paired.annotated.gam | vg view -aj - | wc -l)" "200" "reads with refpos annotations not matching an exclusion regex are let through"
is "$(vg filter -X "[w-z]" paired.annotated.gam | vg view -aj - | wc -l)" "0" "reads with refpos annotations matching an exclusion regex are removed"

rm -f x.vg x.xg paired.gam paired.sam paired.annotated.gam single.gam single.sam filtered.gam filtered.sam
                                                               

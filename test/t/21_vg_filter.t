#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 13

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg  x.vg
vg sim -x x.xg -l 100 -n 5000 -e 0.01 -i 0.001 -a > x.gam

# sanity check: does passing no options preserve input
is $(vg filter x.gam | vg view -a - | jq . | grep mapping | wc -l) 5000 "vg filter with no options preserves input."

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

is "$(vg filter -U x.gam | vg view -aj - | wc -l)" "0" "negating a non-filter results in no reads"

is "$(cat <(vg filter -U -d 123.5 -t 10 paired.gam | vg view -aj -) <(vg filter -d 123.5 -t 10 paired.gam | vg view -aj -)  | wc -l)" "$(vg view -aj paired.gam | wc -l)" "a filter and its complement should form the entire file"

vg index -x g.xg -g g.gcsa -k 16 graphs/refonly-lrc_kir.vg
vg mpmap -x g.xg -g g.gcsa -f reads/grch38_lrc_kir_paired.fq -n dna -B -i -I 10 -D 50 -F GAM -t 1 > g.proper.gam
vg mpmap -x g.xg -g g.gcsa -f reads/grch38_lrc_kir_paired.fq -n dna -B -i -I 10 -D 1 -F GAM -t 1 > g.improper.gam

is "$(vg filter -p g.proper.gam | vg view -aj - | wc -l)" 2 "properly paired read passes proper filter"
is "$(vg filter -p g.improper.gam | vg view -aj - | wc -l)" 0 "improperly paired read fails proper filter"

rm g.xg g.gcsa g.gcsa.lcp g.proper.gam g.improper.gam


is "$(echo '{"sequence": "GATTACA", "name": "read1", "annotation": {"features": ["test"]}, "fragment_next": {"name": "read2"}}{"sequence": "CATTAG", "name": "read2", "fragment_prev":{"name": "read1"}}' | vg view -JGa - | vg filter -F "test" -i - | vg view -aj - | wc -l)" "0" "read pairs can be tropped by feature"
is "$(echo '{"sequence": "GATTACA", "name": "read1", "annotation": {"features": ["test"]}, "fragment_next": {"name": "read2"}}{"sequence": "CATTAG", "name": "read2", "fragment_prev":{"name": "read1"}}' | vg view -JGa - | vg filter -F "test" -I - | vg view -aj - | wc -l)" "2" "read pairs can be kept if only one read fails"

is "$(echo '{"sequence": "GATTACA", "name": "read1"}' | vg view -JGa - | vg filter -P - | vg view -aj - | wc -l)" "0" "unmapped reads can be filtered out"


rm -f x.gam filter_chunk*.gam chunks.bed
rm -f x.vg x.xg paired.gam paired.sam paired.annotated.gam single.gam single.sam filtered.gam filtered.sam
                                                               

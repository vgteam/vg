#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 5

vg view -J -v pileup/tiny.json > tiny.vg

# Compare output of pileup on tiny.vg and pileup/alignment.json
# with pileup/truth.json, which has been manually vetted.
# Will also test some vg view functionality. 
vg view -J -a -G pileup/alignment.json > alignment.gam
vg augment tiny.vg alignment.gam -P tiny.gpu > /dev/null
vg view tiny.gpu -l -j | jq . > tiny.gpu.json
is $(jq --argfile a tiny.gpu.json --argfile b pileup/truth.json -n '($a == $b)') true "vg augment -P produces the expected output for test case on tiny graph."
rm -f alignment.gam tiny.gpu tiny.gpu.json

# Make sure well-supported edits are augmented in
vg view -J -a -G pileup/edits.json > edits.gam
vg augment tiny.vg edits.gam -A edits-embedded.gam > augmented.vg

# We want 3 edits with no sequence per read, and we have 12 reads in this file.
is "$(vg view -aj edits-embedded.gam | jq -c '.path.mapping[].edit[].sequence' | grep null | wc -l)" "36" "augmentation embeds reads fully for well-supported SNPs"
is "$(vg stats -N augmented.vg)" "18" "adding a well-supported SNP by augmentation adds 3 more nodes"

rm -f edits.gam edits-embedded.gam augmented.vg

# Make sure every edit is augmented in
vg view -J -a -G pileup/edit.json > edit.gam
vg augment tiny.vg edit.gam -A edit-embedded.gam > augmented.vg

# This file only has one read.
is "$(vg view -aj edit-embedded.gam | jq -c '.path.mapping[].edit[].sequence' | grep null | wc -l)" "3" "augmentation embeds reads fully for probable errors"
is "$(vg stats -N augmented.vg)" "18" "adding a probable error by augmentation adds 3 more nodes"

rm -f edit.gam edit-embedded.gam augmented.vg

rm -f tiny.vg



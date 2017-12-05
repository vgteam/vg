#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 3

vg view -J -v pileup/tiny.json > tiny.vg

# Compare output of pileup on tiny.vg and pileup/alignment.json
# with pileup/truth.json, which has been manually vetted.
# Will also test some vg view functionality. 
vg view -J -a -G pileup/alignment.json > alignment.gam
vg augment tiny.vg alignment.gam -P tiny.gpu > /dev/null
vg view tiny.gpu -l -j | jq . > tiny.gpu.json
is $(jq --argfile a tiny.gpu.json --argfile b pileup/truth.json -n '($a == $b)') true "vg pileup produces the expected output for test case on tiny graph."
rm -f alignment.gam tiny.gpu tiny.gpu.json

# Make sure every edit is augmented in
vg view -J -a -G pileup/edit.json > edit.gam
vg augment tiny.vg edit.gam -A edit-embedded.gam > augmented.vg

is "$(vg view -aj edit-embedded.gam | jq -c '.path.mapping[].edit[].sequence' | grep null | wc -l)" "3" "augmentation embeds reads fully"
is "$(vg stats -N augmented.vg)" "18" "adding a SNP by augmentation adds 3 more nodes"

rm -f edit.gam edit-embedded.gam augmented.vg

rm -f tiny.vg



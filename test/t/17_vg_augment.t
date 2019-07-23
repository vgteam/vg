#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 13

vg view -J -v pileup/tiny.json > tiny.vg

# Compare output of pileup on tiny.vg and pileup/alignment.json
# with pileup/truth.json, which has been manually vetted.
# Will also test some vg view functionality. 
vg view -J -a -G pileup/alignment.json > alignment.gam
vg augment -a pileup tiny.vg alignment.gam -P tiny.gpu > /dev/null
vg view tiny.gpu -l -j | jq . > tiny.gpu.json
is $(jq --argfile a tiny.gpu.json --argfile b pileup/truth.json -n '($a == $b)') true "vg augment -P produces the expected output for test case on tiny graph."
rm -f alignment.gam tiny.gpu tiny.gpu.json

# Make sure well-supported edits are augmented in
vg view -J -a -G pileup/edits.json > edits.gam
vg augment -a direct tiny.vg edits.gam -A edits-embedded.gam > augmented.vg

# We want 3 edits with no sequence per read, and we have 12 reads in this file.
is "$(vg view -aj edits-embedded.gam | jq -c '.path.mapping[].edit[].sequence' | grep null | wc -l)" "36" "direct augmentation embeds reads fully for well-supported SNPs"
is "$(vg stats -N augmented.vg)" "18" "adding a well-supported SNP by direct augmentation adds 3 more nodes"

rm -f edits.gam edits-embedded.gam augmented.vg

# Make sure every edit is augmented in
vg view -J -a -G pileup/edit.json > edit.gam
vg augment -a direct tiny.vg edit.gam -A edit-embedded.gam > augmented.vg

# This file only has one read.
is "$(vg view -aj edit-embedded.gam | jq -c '.path.mapping[].edit[].sequence' | grep null | wc -l)" "3" "direct augmentation embeds reads fully for probable errors"
is "$(vg stats -N augmented.vg)" "18" "adding a probable error by direct augmentation adds 3 more nodes"

rm -f edit.gam edit-embedded.gam augmented.vg

rm -f tiny.vg

# These are tests that used to be in 14_vg_mod.t

vg construct -m 1000 -r tiny/tiny.fa >t.vg
vg index -k 11 -g t.idx.gcsa -x t.idx.xg t.vg

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | wc -l) 1 "path inclusion does not modify the graph when alignment is a perfect match"

is $(vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | wc -l) 5 "path inclusion with a complex variant introduces the right number of nodes"

# checks that we get a node with the id 4, which is the ref-matching dual to the deletion
is $(vg map -s CAAAAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | grep 4 | grep T | wc -l) 1 "path inclusion works for deletions"

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGCAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | wc -l) 4 "SNPs can be included in the graph"

rm t.vg
rm -rf t.idx.xg t.idx.gcsa

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
vg align -s GGGGGGGAAATTTTCTGGAGTTCTATTATATTCCAAAAAAAAAA t.vg >t.gam
is $(vg augment -i t.vg t.gam | vg view - | grep ^S | grep $(vg augment -i t.vg t.gam | vg stats  -H - | awk '{ print $3}') | cut -f 3) GGGGG "a soft clip at read start becomes a new head of the graph"
is $(vg augment -i t.vg t.gam | vg view - | grep ^S | grep $(vg augment -i t.vg t.gam | vg stats  -T - | awk '{ print $3}') | cut -f 3) AAAAAAAA "a soft clip at read end becomes a new tail of the graph"
rm -rf t.vg t.gam

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg map -x x.xg -g x.gcsa -G small/x-s1337-n100-e0.01-i0.005.gam -t 1 >x.gam
vg augment -Z x.trans -i x.vg x.gam >x.mod.vg
is $(vg view -Z x.trans | wc -l) 1290 "the expected graph translation is exported when the graph is edited"
rm -rf x.vg x.xg x.gcsa x.reads x.gam x.mod.vg x.trans

vg construct -m 1000 -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
is $(vg augment flat.vg 2snp.gam -i | vg mod -D - | vg mod -n - | vg view - | grep ^S | wc -l) 7 "editing the graph with many SNP-containing alignments does not introduce duplicate identical nodes"
rm -f flat.vg flat.gcsa flat.xg 2snp.vg 2snp.xg 2snp.sim 2snp.gam

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 21

vg view -J -v pileup/tiny.json > tiny.vg

# Make sure well-supported edits are augmented in
vg view -J -a -G pileup/edits.json > edits.gam
vg augment -a direct tiny.vg edits.gam -A edits-embedded.gam > augmented.vg

# We want 3 edits with no sequence per read, and we have 12 reads in this file.
is "$(vg view -aj edits-embedded.gam | jq -c '.path.mapping[].edit[].sequence' | grep null | wc -l)" "36" "direct augmentation embeds reads fully for well-supported SNPs"
is "$(vg stats -N augmented.vg)" "18" "adding a well-supported SNP by direct augmentation adds 3 more nodes"

# Run again but with packed logic.  output should be identical with min threshold of 1
vg augment -a direct tiny.vg edits.gam -A edits-embedded.gam -m 1 > augmented.m1.vg
is "$(vg stats -N augmented.m1.vg)" "18" "adding a well-supported SNP by direct augmentation adds 3 more nodes with -m 1"

rm -f edits.gam edits-embedded.gam augmented.vg augmented.m1.vg

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

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTAATATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i -m 2 | vg view - | grep ^S | wc -l) 1 "path inclusion does not modify the graph when alignment has a SNP but doesnt meet the coverage threshold"

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTAATATATTCCAACTCTCTG -V read -d t.idx | vg augment t.vg - -i -m 2 -A read_aug.gam | vg view - | grep ^P | awk '{print $4}' | uniq) "50M" "path inclusion does not modify the included path when alignment has a SNP but doesnt meet the coverage threshold"

is $(vg view -a read_aug.gam | jq. | grep edit | wc) 1 "output GAM has single edit when SNP was filtered out due to coverage"

is $(vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | wc -l) 5 "path inclusion with a complex variant introduces the right number of nodes"

# checks that we get a node with the id 4, which is the ref-matching dual to the deletion
is $(vg map -s CAAAAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | grep 4 | grep T | wc -l) 1 "path inclusion works for deletions"

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGCAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg augment t.vg - -i | vg view - | grep ^S | wc -l) 4 "SNPs can be included in the graph"

rm t.vg
rm -rf t.idx.xg t.idx.gcsa read_aug.gam

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
vg align -s GGGGGGGAAATTTTCTGGAGTTCTATTATATTCCAAAAAAAAAA t.vg >t.gam
is $(vg augment -i t.vg t.gam | vg view - | grep ^S | grep $(vg augment -i t.vg t.gam | vg stats  -H - | awk '{ print $3}') | cut -f 3) GGGGG "a soft clip at read start becomes a new head of the graph"
is $(vg augment -i t.vg t.gam | vg view - | grep ^S | grep $(vg augment -i t.vg t.gam | vg stats  -T - | awk '{ print $3}') | cut -f 3) AAAAAAAA "a soft clip at read end becomes a new tail of the graph"
vg align -s AAATTTTCTGGAGTTCTAT t.vg >> t.gam
vg find -x t.vg -n 9 -c 1 > n9.vg
vg augment n9.vg t.gam -s -A n9_aug.gam > /dev/null
is $(vg view -a n9_aug.gam | wc -l) "1" "augment -s works as desired"
rm -rf t.vg t.gam n9.vg n9_aug.gam

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg map -x x.xg -g x.gcsa -G small/x-s1337-n100-e0.01-i0.005.gam -t 1 >x.gam
vg augment -Z x.trans -i x.vg x.gam >x.mod.vg
is $(vg view -Z x.trans | wc -l) 1288 "the expected graph translation is exported when the graph is edited"
rm -rf x.vg x.xg x.gcsa x.reads x.gam x.mod.vg x.trans

vg construct -m 1000 -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
is $(vg augment flat.vg 2snp.gam -i | vg mod -D - | vg mod -n - | vg view - | grep ^S | wc -l) 7 "editing the graph with many SNP-containing alignments does not introduce duplicate identical nodes"

vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTATCTGGAGTTCTATTATATCCCAACTCTCTG/' | vg view -Fv - >2err.vg
vg sim -l 30 -x 2err.vg -n 10 -a >2err.sim
vg map -g flat.gcsa -x flat.xg -G 2err.sim -k 8 >2err.gam
cat 2snp.gam 2err.gam > 4edits.gam
vg augment flat.vg 2snp.gam | vg view - | grep S | awk '{print $3}' | sort >  2snp_default.nodes
vg augment flat.vg 2snp.gam -m 1 | vg view - | grep S | awk '{print $3}' | sort >  2snp_m1.nodes
diff 2snp_default.nodes 2snp_m1.nodes
is "$?" 0 "augmenting 2 snps with -m 1 produces the same nodes as default"
vg augment flat.vg 4edits.gam -m 11 | vg view - | grep S | awk '{print $3}' | sort > 4edits_m11.nodes
diff 2snp_default.nodes 4edits_m11.nodes
is "$?" 0 "augmenting 2 snps and 2 errors with -m 11 produces the same nodes as with just the snps"

vg augment flat.vg 2snp.gam | vg view - | grep S | awk '{print $3}' | sort > vg_augment.nodes
vg convert flat.vg -p > flat.pg
vg augment flat.pg 2snp.gam | vg convert -v - | vg view - | grep S | awk '{print $3}' | sort > packed_graph_augment.nodes
diff vg_augment.nodes packed_graph_augment.nodes
is "$?" 0 "augmenting a packed graph produces same results as a vg graph"
vg convert flat.vg -a > flat.hg
vg augment flat.hg 2snp.gam | vg convert -v - | vg view - | grep S | awk '{print $3}' | sort > hash_graph_augment.nodes
diff vg_augment.nodes hash_graph_augment.nodes
is "$?" 0 "augmenting a hash graph produces same results as a vg graph"

rm -f flat.vg flat.gcsa flat.xg flat.pg flat.hg 2snp.vg 2snp.xg 2snp.sim 2snp.gam vg_augment.nodes packed_graph_augment.nodes hash_graph_augment.nodes
rm -f 2err.sim 2err.gam 4edits.gam 2snp_default.nodes 2snp_m1.nodes 4edits_m11.nodes

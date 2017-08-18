#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 41

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep "^P" | cut -f 3 | grep -o "[0-9]\+" |  wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep "^S" | wc -l) \
    "vg mod yields a graph with only a particular path"

is $(vg mod -o graphs/orphans.vg | vg view - | wc -l) 8 "orphan edge removal works"

vg construct -r tiny/tiny.fa >t.vg
vg index -k 11 -g t.idx.gcsa -x t.idx.xg t.vg

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 1 "path inclusion does not modify the graph when alignment is a perfect match"

is $(vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 5 "path inclusion with a complex variant introduces the right number of nodes"

# checks that we get a node with the id 4, which is the ref-matching dual to the deletion
is $(vg map -s CAAAAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | grep 4 | grep T | wc -l) 1 "path inclusion works for deletions"

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGCAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 4 "SNPs can be included in the graph"

rm t.vg
rm -rf t.idx.xg t.idx.gcsa

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 - | vg view -g - | sort | md5sum | awk '{ print $1 }') ce5d5ffaa71fea6a25cb4a5b836ccb89 "graph complexity reduction works as expected"

is $( vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 -t 16 - | vg mod -S -l 200 - | vg view - | grep ^S | wc -l) 186 "short subgraph pruning works"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
vg align -s GGGGGGGAAATTTTCTGGAGTTCTATTATATTCCAAAAAAAAAA t.vg >t.gam
is $(vg mod -i t.gam t.vg | vg view - | grep ^S | grep $(vg mod -i t.gam t.vg | vg stats  -H - | awk '{ print $3}') | cut -f 3) GGGGG "a soft clip at read start becomes a new head of the graph"
is $(vg mod -i t.gam t.vg | vg view - | grep ^S | grep $(vg mod -i t.gam t.vg | vg stats  -T - | awk '{ print $3}') | cut -f 3) AAAAAAAA "a soft clip at read end becomes a new tail of the graph"
rm -rf t.vg t.gam

is $(vg mod -n msgas/q_redundant.vg | vg view - | grep ^S | wc -l) 4 "normalization produces the correct number of nodes"

vg mod -n msgas/q_redundant.vg | vg validate -
is $? 0 "normalization produces a valid graph"

vg mod -u msgas/q_redundant.vg | vg validate -
is $? 0 "unchop produces a valid graph"

is $(vg mod -n msgas/q_redundant.vg | vg stats -l - | cut -f 2) 154 "normalization removes redundant sequence in the graph"

is $(vg view -Fv graphs/normalize_me.gfa | vg mod -n - | vg view - | md5sum | cut -f 1 -d\ ) 1dc872cf9cfa9bf110064b37afa25c7b "normalization doesn't introduce cycles and does remove redundancy in bubbles"

# shows that after mod we have == numbers of path annotations and nodes
# in this one-path graph
is $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep "^P" | cut -f 3 | grep -o "[0-9]\+" |wc -l) \
   $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep "^S" | wc -l) \
   "vg mod removes non-path nodes and edge"

set -o pipefail
vg view -Jv reversing/reversing_path.json | vg mod -X 3 - | vg validate -
is "$?" "0" "chopping a graph works correctly with reverse mappings"
set +o pipefail

is $(vg msga -w 20 -f msgas/s.fa | vg mod -X 5 -| vg mod -u - | vg validate - && vg msga -w 20 -f msgas/s.fa | vg mod -X 5 - | vg mod -u - | vg paths -x - | vg view -a - | jq '.sequence' | sort | md5sum | cut -f 1 -d\ ) 2f785068c91dbe84177c1fd679b6f133 "unchop correctly handles paths"

is $(vg view -Jv msgas/inv-mess.json | vg mod -u - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -u - | md5sum | cut -f 1 -d\ ) 99caa2e7716596c7597535a6f0bc9c6e "unchop correctly handles a graph with an inversion"

is "$(vg view -Jv reversing/double_reversing.json | vg mod -u - | vg stats -z - | grep "nodes" | cut -f2)" "1" "unchop handles doubly-reversing edges"

is $(vg view -Jv msgas/inv-mess.json | vg mod -n - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -n - | md5sum | cut -f 1 -d\ ) 02e484f8bc83bf4f37ecda0ad684b235 "normalization works on a graph with an inversion"

vg msga -w 20 -f msgas/s.fa > s.vg
vg msga -g s.vg -s TCAGATTCTCATCCCTCCTCAAGGGCTTCTGTAGCTTTGATGTGGAGTAGTTCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -w 16 -N > s2.vg
vg mod -u s2.vg >/dev/null
is $? 0 "mod successfully unchops a difficult graph"
rm -f s.vg s2.vg

is $(vg msga -w 20 -f msgas/s.fa  | vg mod -r s1 - | vg view - | grep ^P | cut -f 3 | sort | uniq | wc -l) 1 "a single path may be retained"

is $(vg msga -w 20 -f msgas/s.fa | vg mod -r s1 - | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) $(vg msga -w 20 -f msgas/s.fa  | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) "path filtering does not modify the graph"

vg msga -f msgas/l.fa -b a1 -w 16 | vg mod -X 8 - | vg validate -
is $? 0 "chopping self-cycling nodes retains the cycle"

vg mod -U 3 graphs/atgclinv2.vg | vg validate -
is $? 0 "unrolling works and produces a valid graph"

vg mod -U 3 graphs/atgclinv2.vg | vg mod -f 50 - | vg validate -
is $? 0 "unfolding works and produces a valid graph"

is $(comm -12 <(vg mod -d 7 graphs/atgclinv2.vg | vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) <(cat graphs/atgclinv2.vg| vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) | wc -l) $(cat graphs/atgclinv2.vg| vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq | wc -l) "dagify-unroll produces a graph with the same kmers as the original graph"

is $(comm -12 <(vg mod -f 7 graphs/atgclinv2.vg | vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) <(cat graphs/atgclinv2.vg| vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) | wc -l) $(cat graphs/atgclinv2.vg| vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq | wc -l) "unfold produces a graph with the same kmers as the original graph"

vg mod -f 6 graphs/atgclinv2.vg | vg mod -d 5 - | vg validate -
is $? 0 "dagify handles a graph with two strongly connected components"

is $(vg mod -f 6 graphs/atgclinv2.vg | vg mod -d 5 - | vg stats -c - | wc -l) $(vg mod -f 6 graphs/atgclinv2.vg | vg mod -d 5 - | vg stats -N -) "unfold followed by dagify produces a graph with no cycles"

is $(comm -12 <(vg mod -f 7 graphs/atgclinv2.vg | vg mod -d 7 - | vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) <(cat graphs/atgclinv2.vg| vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq ) | wc -l) $(cat graphs/atgclinv2.vg| vg kmers -g -k 7 - | cut -f 1 | grep -v '\#' | grep -v '\$' | sort | uniq | wc -l) "unfold followed by dagify produces a graph with the same kmers as the original graph"

vg mod -w 100 graphs/ununrollable.vg | vg validate -
is $? 0 "dagify unrolls the un-unrollable graph"

vg mod -s graphs/not-simple.vg | vg validate -
is $? 0 "sibling simplification does not disrupt paths"

vg msga -f msgas/cycle.fa -b s1 -w 32 -t 1 | vg mod -N - | vg mod -D -| vg mod -U 10 - | vg mod -c - >c.vg
is $(cat c.vg| vg mod -X 30 - | vg mod -w 100 - | vg stats -N -) 31 "dagify correctly calculates the minimum distance through the unrolled component"
is $(cat c.vg | vg mod -X 10 - | vg mod -w 50 -L 400 - | vg stats -l - | cut -f 2) 422 "dagify only takes one step past our component length limit"
rm -f c.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg sim -s 1337 -n 100 -e 0.01 -i 0.005 -x x.xg -a >x.sim
vg map -x x.xg -g x.gcsa -G x.sim -t 1 >x.gam
vg mod -Z x.trans -i x.gam x.vg >x.mod.vg
is $(vg view -Z x.trans | wc -l) 1280 "the expected graph translation is exported when the graph is edited"
rm -rf x.vg x.xg x.gcsa x.reads x.gam x.mod.vg x.trans

vg construct -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -s 420 -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim >2snp.gam
is $(vg mod -i 2snp.gam flat.vg | vg mod -D - | vg mod -n - | vg view - | grep ^S | wc -l) 7 "editing the graph with many SNP-containing alignments does not introduce duplicate identical nodes"
rm -f flat.vg 2snp.vg 2snp.xg 2snp.sim 2snp.gam

# Note the math (and subsetting) only works out on a flat alleles graph
vg construct -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
vg mod -v small/x.vcf.gz x.vg >x.sample.vg
hom_sites=$(gunzip -c small/x.vcf.gz | grep -v "^#" | grep "1|1" | wc -l)
aug_nodes=$(vg stats x.vg -N)
sample_nodes=$(vg stats x.sample.vg -N)
is ${sample_nodes} $((aug_nodes - hom_sites)) "subsetting a flat-alleles graph to a sample graph works"
rm -f x.vg x.sample.vg 

is "$(vg view -Fv overlaps/two_snvs_assembly1.gfa | vg mod --bluntify - | vg stats -l - | cut -f2)" "315" "bluntifying overlaps results in a graph with duplicated overlapping bases removed"

is "$(vg view -Fv overlaps/two_snvs_assembly4.gfa | vg mod --bluntify - | vg stats -l - | cut -f2)" "335" "bluntifying overlaps works in a more complex graph"

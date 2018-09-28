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

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 - | vg stats -E - ) 285 "graph complexity reduction works as expected"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 -t 16 - | vg mod -S -l 200 - | vg stats -l - | cut -f 2) 983 "short subgraph pruning works"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
vg align -s GGGGGGGAAATTTTCTGGAGTTCTATTATATTCCAAAAAAAAAA t.vg >t.gam
is $(vg mod -i t.gam t.vg | vg view - | grep ^S | grep $(vg mod -i t.gam t.vg | vg stats  -H - | awk '{ print $3}') | cut -f 3) GGGGG "a soft clip at read start becomes a new head of the graph"
is $(vg mod -i t.gam t.vg | vg view - | grep ^S | grep $(vg mod -i t.gam t.vg | vg stats  -T - | awk '{ print $3}') | cut -f 3) AAAAAAAA "a soft clip at read end becomes a new tail of the graph"
rm -rf t.vg t.gam

is $(vg mod -n msgas/q_redundant.vg | vg view - | grep ^S | wc -l) 4 "normalization produces the correct number of nodes"

is $(vg view -vF graphs/redundant-snp.gfa | vg mod -n - | vg view - | grep ^S | wc -l) 4 "normalization removes redundant SNP alleles"

vg mod -n msgas/q_redundant.vg | vg validate -
is $? 0 "normalization produces a valid graph"

vg mod -u msgas/q_redundant.vg | vg validate -
is $? 0 "unchop produces a valid graph"

is $(vg mod -n msgas/q_redundant.vg | vg stats -l - | cut -f 2) 154 "normalization removes redundant sequence in the graph"

is $(vg view -Fv graphs/normalize_me.gfa | vg mod -n - | vg view - | sort | md5sum | cut -f 1 -d\ ) 4212cbf054655a4608f3b0fee4b4a59a "normalization doesn't introduce cycles and does remove redundancy in bubbles"

# shows that after mod we have == numbers of path annotations and nodes
# in this one-path graph
is $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep "^P" | cut -f 3 | grep -o "[0-9]\+" |wc -l) \
   $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep "^S" | wc -l) \
   "vg mod removes non-path nodes and edge"

set -o pipefail
vg view -Jv reversing/reversing_path.json | vg mod -X 3 - | vg validate -
is "$?" "0" "chopping a graph works correctly with reverse mappings"
set +o pipefail

is $(vg msga -w 20 -f msgas/s.fa | vg mod -X 5 -| vg mod -u - | vg validate - && vg msga -w 20 -f msgas/s.fa | vg mod -X 5 - | vg mod -u - | vg paths -v - -X | vg view -a - | jq '.sequence' | sort | md5sum | cut -f 1 -d\ ) 2f785068c91dbe84177c1fd679b6f133 "unchop correctly handles paths"

is "$(vg view -Jv msgas/inv-mess.json | vg mod -u - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -u - | vg view - | sort | diff - msgas/inv-mess-unchopped.gfa)" "" "unchop correctly handles a graph with an inversion"

is "$(vg view -Jv reversing/double_reversing.json | vg mod -u - | vg stats -z - | grep "nodes" | cut -f2)" "1" "unchop handles doubly-reversing edges"

is "$(vg view -Jv msgas/inv-mess.json | vg mod -n - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -n - | vg view - | sort | diff - msgas/inv-mess-normalized.gfa)" "" "normalization works on a graph with an inversion"

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

vg msga -g <(vg msga -f msgas/cycle.fa -b s1 -w 32 -t 1 | vg mod -D - | vg mod -U 10 -) -f msgas/cycle.fa -t 1 | vg mod -N - | vg mod -D - | vg mod -U 10 - >c.vg
is $(cat c.vg | vg mod -w 100 - | vg stats -N -) 4 "dagify correctly calculates the minimum distance through the unrolled component"
is $(cat c.vg | vg mod -w 100 - | vg stats -l - | cut -f 2) 200 "dagify produces a graph of the correct size"
rm -f c.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg map -x x.xg -g x.gcsa -G small/x-s1337-n100-e0.01-i0.005.gam -t 1 >x.gam
vg mod -Z x.trans -i x.gam x.vg >x.mod.vg
is $(vg view -Z x.trans | wc -l) 1288 "the expected graph translation is exported when the graph is edited"
rm -rf x.vg x.xg x.gcsa x.reads x.gam x.mod.vg x.trans

vg construct -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
is $(vg mod -i 2snp.gam flat.vg | vg mod -D - | vg mod -n - | vg view - | grep ^S | wc -l) 7 "editing the graph with many SNP-containing alignments does not introduce duplicate identical nodes"
rm -f flat.vg flat.gcsa flat.xg 2snp.vg 2snp.xg 2snp.sim 2snp.gam

# Note the math (and subsetting) only works out on a flat alleles graph
vg construct -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
vg mod -v small/x.vcf.gz x.vg >x.sample.vg
is $(vg stats x.sample.vg -l | cut -f 2) 1041 "subsetting a flat-alleles graph to a sample graph works"
rm -f x.vg x.sample.vg 

is $(vg mod -M 5 jumble/j.vg|  vg stats -s - | wc -l) 7 "removal of high-degree nodes results in the expected number of subgraphs"

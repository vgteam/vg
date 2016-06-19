#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 37

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^P | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^S | wc -l) \
    "vg mod yields a graph with only a particular path"

is $(vg mod -o graphs/orphans.vg | vg view - | wc -l) 8 "orphan edge removal works"

vg construct -r tiny/tiny.fa >t.vg
vg index -s -k 11 -d t.idx t.vg

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 1 "path inclusion does not modify the graph when alignment is a perfect match"

is $(vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 5 "path inclusion with a complex variant introduces the right number of nodes"

# checks that we get a node with the id 4, which is the ref-matching dual to the deletion
is $(vg map -s CAAAAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | grep 4 | grep T | wc -l) 1 "path inclusion works for deletions"

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGCAGTTCTATTATATTCCAACTCTCTG -d t.idx | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 4 "SNPs can be included in the graph"

rm t.vg
rm -rf t.idx

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 - | vg view -g - | sort | md5sum | awk '{ print $1 }') 7cd68c575e236202fab8522724d7e9df "graph complexity reduction works as expected"

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
is $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep ^P |wc -l) \
   $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep ^S |wc -l) \
   "vg mod removes non-path nodes and edge"

is "$(vg view -Jv reversing/reversing_path.json | vg mod -X 3 - | vg validate - && echo 'Graph is valid')" "Graph is valid" "chopping a graph works correctly with reverse mappings"

is $(vg msga -B 20 -f msgas/s.fa | vg mod -X 5 -| vg mod -u - | vg validate - && vg msga -B 20 -f msgas/s.fa | vg mod -X 5 - | vg mod -u - | vg paths -x - | vg view -a - | jq '.sequence' | sort | md5sum | cut -f 1 -d\ ) 2f785068c91dbe84177c1fd679b6f133 "unchop correctly handles paths"

is $(vg view -Jv msgas/inv-mess.json | vg mod -u - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -u - | md5sum | cut -f 1 -d\ ) 0e7a50bb7367d9f84fbc9bd78378d70f "unchop correctly handles a graph with an inversion"

is $(vg view -Jv msgas/inv-mess.json | vg mod -n - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -n - | md5sum | cut -f 1 -d\ ) 84138fe8bd1015fb6b80278c2ed8f7c6 "normalization works on a graph with an inversion"

vg msga -g s.vg -s TCAGATTCTCATCCCTCCTCAAGGGCTTCTGTAGCTTTGATGTGGAGTAGTTCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -B 16 -Nz | vg mod -u - >/dev/null
is $? 0 "mod successfully unchops a difficult graph"

is $(vg msga -B 20 -f msgas/s.fa  | vg mod -r s1 - | vg view - | grep ^P | cut -f 3 | sort | uniq | wc -l) 1 "a single path may be retained"

is $(vg msga -B 20 -f msgas/s.fa | vg mod -r s1 - | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) $(vg msga -B 20 -f msgas/s.fa  | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) "path filtering does not modify the graph"

vg msga -f msgas/l.fa -b a1 -B 8 | vg mod -X 8 - | vg validate -
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

vg msga -f msgas/cycle.fa -b s1 -B 20 -t 1 | vg mod -D - | vg mod -n - | vg mod -n - >c.vg
is $(cat c.vg| vg mod -X 30 - | vg mod -w 100 - | vg stats -N -) 36 "dagify correctly calculates the minimum distance through the unrolled component"
is $(cat c.vg | vg mod -X 10 - | vg mod -w 50 -L 400 - | vg stats -l - | cut -f 2) 400 "dagify only takes one step past our component length limit"
rm -f c.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 11 x.vg
vg sim -s 1337 -n 100 -x x.xg >x.reads
vg map -x x.xg -g x.gcsa -r x.reads -L 10 >x.gam
vg mod -Z x.trans -i x.gam x.vg >x.mod.vg
is $(vg view -Z x.trans | jq -c --sort-keys . | sort | md5sum | cut -f 1 -d\ ) be1b8d52a3928d91efcec802205e8971 "a valid graph translation is exported when the graph is edited"
rm -rf x.vg x.xg x.gcsa x.reads x.gam x.mod.vg x.trans

vg construct -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -s 420 -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg map -V flat.vg -k 8 -G 2snp.sim >2snp.gam
is $(vg mod -i 2snp.gam flat.vg | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f -1 -d\ ) d47169ce8fb4251904d1d2238a44555c "editing the graph with many SNP-containing alignments does not introduce duplicate identical nodes"
rm -f flat.vg 2snp.vg 2snp.xg 2snp.sim 2snp.gam

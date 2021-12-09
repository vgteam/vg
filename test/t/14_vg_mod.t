#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order

plan tests 32

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep "^P" | cut -f 3 | grep -o "[0-9]\+" |  wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep "^S" | wc -l) \
    "vg mod yields a graph with only a particular path"

is $(vg mod graphs/orphans.vg | vg view - | wc -l) 8 "orphan edge removal is automatic"

is $(vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 - | vg stats -E - ) 285 "graph complexity reduction works as expected"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 -t 16 - | vg mod -S -l 200 - | vg stats -l - | cut -f 2) 983 "short subgraph pruning works"

is $(vg mod -U 10 msgas/q_redundant.vg | vg view - | grep ^S | wc -l) 4 "normalization produces the correct number of nodes"

is $(vg view -vF graphs/redundant-snp.gfa | vg mod -n - | vg view - | grep ^S | wc -l) 4 "normalization removes redundant SNP alleles"

vg mod -n msgas/q_redundant.vg | vg validate -
is $? 0 "normalization produces a valid graph"

vg mod -U 10 msgas/q_redundant.vg | vg validate -
is $? 0 "looped normalization produces a valid graph"

vg mod -u msgas/q_redundant.vg | vg validate -
is $? 0 "unchop produces a valid graph"

is $(vg mod -U 10 msgas/q_redundant.vg | vg stats -l - | cut -f 2) 154 "normalization removes redundant sequence in the graph"

is $(vg view -Fv graphs/normalize_me.gfa | vg mod -U 10 - | vg view - | sort | md5sum | cut -f 1 -d\ ) $(md5sum graphs/normalize_me.norm.gfa | cut -f 1 -d\ ) "normalization doesn't introduce cycles and does remove redundancy in bubbles"

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

is "$(vg view -Jv msgas/inv-mess.json | vg mod -U 10 - | vg validate - && vg view -Jv msgas/inv-mess.json | vg mod -U 10 - | vg view - | sort | diff - msgas/inv-mess-normalized.gfa)" "" "normalization works on a graph with an inversion"

vg msga -w 20 -f msgas/s.fa > s.vg
vg msga -g s.vg -s TCAGATTCTCATCCCTCCTCAAGGGCTTCTGTAGCTTTGATGTGGAGTAGTTCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -w 16 -N > s2.vg
vg mod -u s2.vg >/dev/null
is $? 0 "mod successfully unchops a difficult graph"
rm -f s.vg s2.vg

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

vg msga -g <(vg msga -f msgas/cycle.fa -b s1 -w 19 -O 18 -k 4 -t 1 | vg paths -d -v - | vg mod -U 10 -) -f msgas/cycle.fa -t 1 | vg mod -N - | vg paths -d -v - | vg mod -U 10 - >c.vg
is $(cat c.vg | vg mod -w 100 - | vg stats -N -) 4 "dagify correctly calculates the minimum distance through the unrolled component"
is $(cat c.vg | vg mod -w 100 - | vg stats -l - | cut -f 2) 200 "dagify produces a graph of the correct size"
rm -f c.vg

# Note the math (and subsetting) only works out on a flat alleles graph
vg construct -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
vg mod -v small/x.vcf.gz x.vg >x.sample.vg
is "$(vg stats x.sample.vg -l | cut -f 2)" "1041" "subsetting a flat-alleles graph to a sample graph works"
rm -f x.vg x.sample.vg 

is $(vg mod -M 5 jumble/j.vg|  vg stats -s - | wc -l) 7 "removal of high-degree nodes results in the expected number of subgraphs"

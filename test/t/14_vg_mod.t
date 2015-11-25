#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 20

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^P | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^S | wc -l) \
    "vg mod yields a graph with only a particular path"

is $(vg mod -o graphs/orphans.vg | vg view - | wc -l) 8 "orphan edge removal works"

vg construct -r tiny/tiny.fa >t.vg
vg index -s -k 11 t.vg

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 1 "path inclusion does not modify the graph when alignment is a perfect match"

is $(vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 5 "path inclusion with a complex variant introduces the right number of nodes"

# checks that we get a node with the id 4, which is the ref-matching dual to the deletion
is $(vg map -s CAAAAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | grep 4 | grep T | wc -l) 1 "path inclusion works for deletions"

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGCAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 4 "SNPs can be included in the graph"

rm t.vg
rm -rf t.vg.index

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 10 -e 3 - | vg view -g - | sort | md5sum | awk '{ print $1 }') a395b6f07549724b68f236dee5e0e084 "graph complexity reduction works as expected"

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

is $(vg view -v graphs/normalize_me.gfa | vg mod -n - | vg view - | md5sum | cut -f 1 -d\ ) a38ac699119df12630ce747e3bf6b5fe "normalization doesn't introduce cycles and does remove redundancy in bubbles"

# shows that after mod we have == numbers of path annotations and nodes
# in this one-path graph
is $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep ^P |wc -l) \
   $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep ^S |wc -l) \
   "vg mod removes non-path nodes and edge"

is "$(vg view -Jv reversing/reversing_path.json | vg mod -X 3 - | vg validate - && echo 'Graph is valid')" "Graph is valid" "chopping a graph works correctly with reverse mappings"

is $(vg view -Jv msgas/inv-mess.json | vg mod -u - | md5sum | cut -f 1 -d\ ) 9684fb6d14ffe5cb8e21cc11abbf04d0 "unchop correctly handles a graph with an inversion"

is $(vg view -Jv msgas/inv-mess.json | vg mod -n - | md5sum | cut -f 1 -d\ ) 29a99770b1d3a7cdbc0ceb73159d6c1f "normalization works on a graph with an inversion"

vg msga -g s.vg -s TCAGATTCTCATCCCTCCTCAAGGGCTTCT$(revcomp AACTACTCCACATCAAAGCTAC)CCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -k 16 -B 16 -Nz | vg mod -u - >/dev/null
is $? 0 "mod successfully unchops a difficult graph"

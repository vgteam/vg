#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 17

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

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 8 -e 3 - | vg view -g - | sort | md5sum | awk '{ print $1 }') 7c426b504a33ba7985372b978650c3dc "graph complexity reduction works as expected"

is $( vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 8 -e 3 -t 16 - | vg mod -S -l 200 - | vg view - | grep ^S | wc -l) 152 "short subgraph pruning works"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
vg align -s GGGGGGGAAATTTTCTGGAGTTCTATTATATTCCAAAAAAAAAA t.vg >t.gam
is $(vg mod -i t.gam t.vg | vg view - | grep ^S | grep $(vg mod -i t.gam t.vg | vg stats  -H - | awk '{ print $3}') | cut -f 3) GGGGG "a soft clip at read start becomes a new head of the graph"
is $(vg mod -i t.gam t.vg | vg view - | grep ^S | grep $(vg mod -i t.gam t.vg | vg stats  -T - | awk '{ print $3}') | cut -f 3) AAAAAAAA "a soft clip at read end becomes a new tail of the graph"
rm -rf t.vg t.gam

is $(vg mod -n msgas/q_redundant.vg | vg view - | grep ^S | wc -l) 4 "normalization produces the correct number of nodes"

is $(vg mod -n msgas/q_redundant.vg | vg stats -l - | cut -f 2) 154 "normalization removes redundant sequence in the graph"

is $(vg view -v graphs/normalize_me.gfa | vg mod -n - | vg view - | md5sum | cut -f 1 -d\ ) a38ac699119df12630ce747e3bf6b5fe "normalization doesn't introduce cycles and does remove redundancy in bubbles"

# shows that after mod we have == numbers of path annotations and nodes
# in this one-path graph
is $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep ^P |wc -l) \
   $(vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg mod -N - | vg view - | grep ^S |wc -l) \
   "vg mod removes non-path nodes and edge"

is $(vg view -v msgas/inv-mess.gfa | vg mod -u - | md5sum | cut -f 1 -d\ ) 75e29bc73b25beaf0450a6aeacef4f67 "unchop correctly handles a graph with an inversion"

is $(vg view -v msgas/inv-mess.gfa | vg mod -n - | md5sum | cut -f 1 -d\ ) 5b1c2a8467f5f1b1b1c8da10f62d48ab "normalization works on a graph with an inversion"

vg msga -g s.vg -s TCAGATTCTCATCCCTCCTCAAGGGCTTCT$(revcomp AACTACTCCACATCAAAGCTAC)CCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -k 16 -B 16 -Nz | vg mod -u - >/dev/null
is $? 0 "mod successfully unchops a difficult graph"

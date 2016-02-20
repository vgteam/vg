#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 16

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -s x.vg
is $? 0 "indexing nodes and edges of graph"

# note that we use "negatives" here even if it isn't so by default
vg index -n -k 11 x.vg
is $? 0 "indexing 11mers"

node_matches=$(vg find -k TAAGGTTTGAA -c 0 x.vg | vg view -g - | grep "^S" | cut -f 2 | grep -E '1$|2$|9$|5$|6$|8$' | wc -l)
is $node_matches 6 "all expected nodes found via kmer find"

edge_matches=$(vg find -k TAAGGTTTGAA -c 0 x.vg | vg view -g - | grep "^L" | cut -f 2 | grep -E '1$|2$|8$|5$|6$' | wc -l)
is $edge_matches 5 "all expected edges found via kmer find"

is $(vg find -n 2 -n 3 -c 1 x.vg | vg view -g - | wc -l) 15 "multiple nodes can be picked using vg find"

is $(vg find -S AGGGCTTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGT x.vg | vg view - | wc -l) 33 "vg find returns a correctly-sized graph when seeking a sequence"

is $(vg find -S AGGGCTTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGT -j 11 x.vg | vg view - | wc -l) 33 "vg find returns a correctly-sized graph when using jump-kmers"

is $(vg find -p x:0-100 x.vg | vg view -g - | wc -l) 42 "vg find returns a subgraph corresponding to particular reference coordinates"

is $(vg find -p x:0-100 x.vg | vg view -j - | jq ".node[].sequence" | tr -d '"\n' | wc -c) 100 "vg find returns a path of the correct length"

is $(vg find -p x:0-100 -c 1 x.vg | vg view -g - | wc -l) 70 "larger graph is returned when the reference path is queried with context"

is $(vg find -p x -c 10 x.vg | vg view -g - | wc -l) $(vg view -g x.vg | wc -l) "entire graph is returned when the reference path is queried with context"

is $(vg find -s 10 x.vg | wc -l) 1 "we can find edges on start"

is $(vg find -e 10 x.vg | wc -l) 1 "we can find edges on end"

rm -rf x.vg.index

vg index -x x.idx x.vg 2>/dev/null
is $(vg find -x x.idx -p x:200-300 -c 2 | vg view - | grep CTACTGACAGCAGA | cut -f 2) 72 "a path can be queried from the xg index"
is $(vg find -x x.idx -n 203 -c 1 | vg view - | grep CTACCCAGGCCATTTTAAGTTTCCTGT | wc -l) 1 "a node near another can be obtained using context from the xg index"

rm -f x.idx
rm -f x.vg


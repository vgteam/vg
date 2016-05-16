#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 24

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -s -d x.idx x.vg
is $? 0 "indexing nodes and edges of graph"

# note that we use "negatives" here even if it isn't so by default
vg index -n -k 11 -d x.idx x.vg
is $? 0 "indexing 11mers"

node_matches=$(vg find -k TAAGGTTTGAA -c 0 -d x.idx | vg view -g - | grep "^S" | cut -f 2 | grep -E '1$|2$|9$|5$|6$|8$' | wc -l)
is $node_matches 6 "all expected nodes found via kmer find"

edge_matches=$(vg find -k TAAGGTTTGAA -c 0 -d x.idx | vg view -g - | grep "^L" | cut -f 2 | grep -E '1$|2$|8$|5$|6$' | wc -l)
is $edge_matches 5 "all expected edges found via kmer find"

is $(vg find -n 2 -n 3 -c 1 -d x.idx | vg view -g - | wc -l) 15 "multiple nodes can be picked using vg find"

is $(vg find -S AGGGCTTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGT -d x.idx | vg view - | wc -l) 33 "vg find returns a correctly-sized graph when seeking a sequence"

is $(vg find -S AGGGCTTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGT -j 11 -d x.idx | vg view - | wc -l) 33 "vg find returns a correctly-sized graph when using jump-kmers"

is $(vg find -p x:0-100 -d x.idx | vg view -g - | wc -l) 42 "vg find returns a subgraph corresponding to particular reference coordinates"

is $(vg find -p x:0-100 -d x.idx | vg view -j - | jq ".node[].sequence" | tr -d '"\n' | wc -c) 100 "vg find returns a path of the correct length"

is $(vg find -p x:0-100 -c 1 -d x.idx | vg view -g - | wc -l) 70 "larger graph is returned when the reference path is queried with context"

is $(vg find -p x -c 10 -d x.idx | vg view -g - | wc -l) $(vg view -g x.vg | wc -l) "entire graph is returned when the reference path is queried with context"

is $(vg find -s 10 -d x.idx | wc -l) 1 "we can find edges on start"

is $(vg find -e 10 -d x.idx | wc -l) 1 "we can find edges on end"

rm -rf x.idx

vg index -x x.idx x.vg 2>/dev/null
is $(vg find -x x.idx -p x:200-300 -c 2 | vg view - | grep CTACTGACAGCAGA | cut -f 2) 72 "a path can be queried from the xg index"
is $(vg find -x x.idx -n 203 -c 1 | vg view - | grep CTACCCAGGCCATTTTAAGTTTCCTGT | wc -l) 1 "a node near another can be obtained using context from the xg index"

vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(( for seq in $(vg sim -l 50 -n 100 -x x.xg); do vg find -M $seq -g x.gcsa; done ) | jq length | grep ^1$ | wc -l) 100 "each perfect read contains one maximal exact match"

vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(vg find -n 2 -n 3 -c 1 -L -x x.idx | vg view -g - | wc -l) 15 "vg find -L finds same number of nodes (with -c 1)"
rm -f x.idx x.xg x.gcsa x.gcsa.lcp x.vg

vg index -x m.xg inverting/m.vg
is $(vg find -n 174 -c 200 -L -x m.xg | vg view -g - | grep S | wc -l) 7 "vg find -L only follows alternating paths"
is $(vg find -n 2308 -c 10 -L -x m.xg | vg view -g - | grep S | wc -l) 10 "vg find -L tracks length"
is $(vg find -n 2315 -n 183 -n 176 -c 1 -L -x m.xg | vg view -g - | grep S | wc -l) 7 "vg find -L works with more than one input node"
rm m.xg

vg construct -rmem/h.fa >h.vg
vg index -g h.gcsa -k 16 h.vg
is $(vg find -M ACCGTTAGAGTCAG -g h.gcsa) '[["ACC",["1:-32"]],["CCGTTAG",["1:5"]],["GTTAGAGT",["1:19"]],["TAGAGTCAG",["1:40"]]]' "we find the 4 canonical SMEMs from @lh3's bwa mem poster"
rm -f h.gcsa h.gcsa.lcp h.vg

vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 >giab.vg
vg index -x giab.xg -g giab.gcsa -k 11 giab.vg
is $(vg find -M ATTCATNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) a7bce59dd511e6fb003720b8d5a788a0 "we can find the right MEMs for a sequence with Ns"
is $(vg find -M ATTCATNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) $(vg find -M ATTCATNNNNNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) "we find the same MEMs sequences with different lengths of Ns"
rm -f giab.vg giab.xg giab.gcsa


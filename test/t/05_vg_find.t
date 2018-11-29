#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 25

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -x x.xg x.vg 2>/dev/null
is $(vg find -x x.xg -p x:200-300 -c 2 | vg view - | grep CTACTGACAGCAGA | cut -f 2) 72 "a path can be queried from the xg index"
is $(vg find -x x.xg -n 203 -c 1 | vg view - | grep CTACCCAGGCCATTTTAAGTTTCCTGT | wc -l) 1 "a node near another can be obtained using context from the xg index"

vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(( for seq in $(vg sim -l 50 -n 100 -x x.xg); do vg find -M $seq -g x.gcsa; done ) | jq length | grep ^1$ | wc -l) 100 "each perfect read contains one maximal exact match"

vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(vg find -n 1 -n 3 -D -x x.xg ) 8 "vg find -D finds approximate distance between 2 adjacent nodes"
is $(vg find -n 1 -n 2 -D -x x.xg ) 0 "vg find -D finds approximate distance between node and adjacent snp"
is $(vg find -n 17 -n 20 -D -x x.xg ) 7 "vg find -D jumps deletion"
# The correct distance is 6 still from the SNP, but renumbering ref vs alt makes
# the heuristic add 1 getting to the reference path
# TODO: improve heuristic
is $(vg find -n 16 -n 20 -D -x x.xg ) 19 "vg find -D jumps deletion from snp"

is $(vg find -n 2 -n 3 -c 1 -L -x x.xg | vg view -g - | grep "^S" | wc -l) 5 "vg find -L finds same number of nodes (with -c 1)"

is $(vg find -r 6:2 -L -x x.xg | vg view -g - | grep S | wc -l) 3 "vg find -L works with -r "

rm -f x.idx x.xg x.gcsa x.gcsa.lcp x.vg

vg index -x m.xg inverting/m.vg
is $(vg find -n 174 -c 200 -L -x m.xg | vg view -g - | grep S | wc -l) 7 "vg find -L only follows alternating paths"
is $(vg find -n 2308 -c 10 -L -x m.xg | vg view -g - | grep S | wc -l) 10 "vg find -L tracks length"
is $(vg find -n 2315 -n 183 -n 176 -c 1 -L -x m.xg | vg view -g - | grep S | wc -l) 7 "vg find -L works with more than one input node"
rm m.xg

vg construct -m 1000 -rmem/h.fa >h.vg
vg index -g h.gcsa -k 16 h.vg
is $(vg find -M ACCGTTAGAGTCAG -g h.gcsa) '[["ACC",["1:-32"]],["CCGTTAG",["1:5"]],["GTTAGAGT",["1:19"]],["TAGAGTCAG",["1:40"]]]' "we find the 4 canonical SMEMs from @lh3's bwa mem poster"
rm -f h.gcsa h.gcsa.lcp h.vg

vg construct -m 1000 -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 >giab.vg
vg index -x giab.xg -g giab.gcsa -k 11 giab.vg
is $(vg find -M ATTCATNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) $(md5sum correct/05_vg_find/28.txt | cut -f -1 -d\ ) "we can find the right MEMs for a sequence with Ns"
is $(vg find -M ATTCATNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) $(vg find -M ATTCATNNNNNNNNAGTTAA -g giab.gcsa | md5sum | cut -f -1 -d\ ) "we find the same MEMs sequences with different lengths of Ns"
rm -f giab.vg giab.xg giab.gcsa

#vg construct -r mem/r.fa > r.vg
#vg index -x r.xg -g r.gcsa -k 16 r.vg
#is $(vg find -M GAGGCAGTGAAGAGATCGTGGGAGGGAC -R 10 -Z 10 -g r.gcsa -x r.xg) '[["GAGGCAGTGAAGAGATCGTGGGAG",["1:20"]],["GCAGTGAAGAGATCGTGGGAGGGAC",["1:43"]],["CAGTGAAGAGATCGTGGGAGG",["1:0"]]]' "we can find sub-MEMs and only return hits that are not within containing MEMs"
#is $(vg find -M GAGGCAGTGAAGAGATCGTGGGAGGGAC -R 10 -Z 10 -f -g r.gcsa -x r.xg) '[["GAGGCAGTGAAGAGATCGTGGGAG",["1:20"]],["GCAGTGAAGAGATCGTGGGAGGGAC",["1:43"]],["CAGTGAAGAGATCGTGGGAGG",["1:0"]]]' "we get the same (sufficiently long) MEMs with the fast sub-MEM option"
#rm -rf r.vg r.xg r.gcsa* mem/r.fa.fai

#vg view -J mem/s.json -v > s.vg
#vg index -x s.xg -g s.gcsa -k 16 s.vg
#is $(vg find -M ACGTGCCGTTAGCCAGTGGGTTAG -R 10 -Z 10 -x s.xg -g s.gcsa) '[["ACGTGCCGTTAGCCAGTGGGTTAG",["3:11"]],["AGCCAGTGGGTTA",["1:0","2:0"]]]' "we can find sub-MEMs in a hard case of de-duplication"
#is $(vg find -M ACGTGCCGTTAGCCAGTGGGTTAG -R 10 -Z 10 -f -x s.xg -g s.gcsa) '[["ACGTGCCGTTAGCCAGTGGGTTAG",["3:11"]],["AGCCAGTGGGTTA",["1:0","2:0"]]]' "we can find the same (sufficiently long) sub-MEMs the hard de-duplication case with the fast sub-MEM option"
#rm -rf s.vg s.xg s.gcsa*

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 11 x.vg
vg map -x x.xg -g x.gcsa -T small/x-s1337-n100.reads >x.gam
vg index -d x.db -N x.gam
is $(vg find -o 127 -d x.db | vg view -a - | wc -l) 6 "the index can return the set of alignments mapping to a particular node"
is $(vg find -A <(vg find -N <(seq 37 52 ) -x x.xg ) -d x.db | vg view -a - | wc -l) 15 "a subgraph query may be used to obtain a particular subset of alignments"

vg gamsort -i x.sorted.gam.gai x.gam > x.sorted.gam
is $(vg find -o 127 --sorted-gam x.sorted.gam | vg view -a - | wc -l) 6 "the GAM index can return the set of alignments mapping to a node"
is $(vg find -A <(vg find -N <(seq 37 52 ) -x x.xg ) --sorted-gam x.sorted.gam | vg view -a - | wc -l) 15 "a subgraph query may be used to obtain a particular subset of alignments from a sorted GAM"

rm -rf x.db x.gam x.sorted.gam x.sorted.gam.gai

is $(vg find -G small/x-s1337-n1.gam -x x.xg | vg view - | grep ATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTT | wc -l) 1 "the index can be queried using GAM alignments"

rm -rf x.vg x.xg x.gcsa

vg construct -m 1000 -r tiny/tiny.fa -v tiny/tiny.vcf.gz >tiny.vg
vg index -x tiny.xg tiny.vg 
is $(vg find -x tiny.xg -n 12 -n 13 -n 14 -n 15 | vg view - | grep ^L | wc -l) 4 "find gets connected edges between queried nodes by default"
echo 12 13 >get.nodes
echo 14 >>get.nodes
echo 15 >>get.nodes
is $(vg find -x tiny.xg -N get.nodes | vg view - | grep ^S | wc -l) 4 "find gets nodes provided in a node file list"
rm -rf tiny.xg tiny.vg get.nodes

echo '{"node": [{"id": 1, "sequence": "A"}, {"id": 2, "sequence": "A"}], "edge": [{"from": 1, "to": 2}], "path": [{"name": "ref", "mapping": [{"position": {"node_id": 1}}]}]}' | vg view -Jv - >test.vg
vg index -x test.xg test.vg
is "$(vg find -x test.xg -p ref:0 -c 10 | vg view -j - | jq '.edge | length')" "1" "extracting by path adds no extra edges"

rm -f test.vg test.xg

vg construct -m 1000 -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a 2> /dev/null | vg view -j - | sed s/_alt/alt/g | vg view -Jv - >w.vg
vg index -x w.xg w.vg
is $(( cat w.vg | vg mod -D - ; vg find -x w.xg -Q alt ) | vg paths -L -v - | wc -l) 38 "pattern based path extraction works"
rm -f w.xg w.vg

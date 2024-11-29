#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 29

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -x x.xg x.vg 2>/dev/null
is $(vg find -x x.xg -p x:200-300 -c 2 | vg view - | grep CTACTGACAGCAGA | cut -f 2) 72 "a path can be queried from the xg index"
# todo: I complain because context expansion no longer puts ranks in paths, so edge check does't see that path is discontinuous
is $(vg find -x x.xg -n 203 -c 1 | vg view - | grep CTACCCAGGCCATTTTAAGTTTCCTGT | wc -l) 1 "a node near another can be obtained using context from the xg index"

vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(( for seq in $(vg sim -l 50 -n 100 -x x.xg); do vg find -M $seq -g x.gcsa; done ) | jq length | grep ^1$ | wc -l) 100 "each perfect read contains one maximal exact match"

vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(vg find -n 1 -n 3 -D -x x.xg ) 8 "vg find -D finds approximate distance between 2 adjacent node starts"
is $(vg find -n 1 -n 2 -D -x x.xg ) 8 "vg find -D finds approximate distance between node start and adjacent snp"
is $(vg find -n 17 -n 20 -D -x x.xg ) 7 "vg find -D jumps deletion"
is $(vg find -n 16 -n 20 -D -x x.xg ) 7 "vg find -D jumps deletion from other allele in snp"

is $(vg find -n 2 -n 3 -c 1 -L -x x.xg | vg view -g - | grep "^S" | wc -l) 5 "vg find -L finds same number of nodes (with -c 1)"

is $(vg find -r 6:5 -L -x x.xg | vg view -g - | grep S | wc -l) 3 "vg find -L works with -r.  it scans from start position of first node in range "

rm -f x.idx x.xg x.gcsa x.gcsa.lcp x.vg

vg index -x m.xg inverting/m.vg
is $(vg find -n 2308 -c 10 -L -x m.xg | vg view -g - | grep S | wc -l) 5 "vg find -L tracks length from start position of input node"
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
rm -f giab.vg giab.xg giab.gcsa{,.lcp}

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

vg gamsort -i x.sorted.gam.gai x.gam > x.sorted.gam
is $(vg find -o 127 --sorted-gam x.sorted.gam | vg view -a - | wc -l) 6 "the GAM index can return the set of alignments mapping to a node"
is $(vg find -A <(vg find -N <(seq 37 52 ) -x x.xg ) --sorted-gam x.sorted.gam | vg view -a - | wc -l) 15 "a subgraph query may be used to obtain a particular subset of alignments from a sorted GAM"

rm -rf x.gam x.sorted.gam x.sorted.gam.gai

is $(vg find -G small/x-s1337-n1.gam -x x.xg | vg view - | grep ATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTT | wc -l) 1 "the index can be queried using GAM alignments"

rm -rf x.vg x.xg x.gcsa{,.lcp}

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
cat w.vg | vg paths -d -v - > part1.tmp
vg find -x w.xg -Q alt > part2.tmp
is "$(vg combine -c part1.tmp part2.tmp | vg paths -L -v - | wc -l)" "$(vg view -j w.vg | jq -c '.path[] | select(.name | contains("alt"))' | wc -l)" "pattern based path extraction works"
rm -f w.xg w.vg part1.tmp part2.tmp

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg
vg index -x t.xg t.vg
is $(vg find -x t.xg -E -p x:30-35 | vg view - | grep ^S | wc -l) 4 "path DAG range query works"

vg find -x t.xg -E -p x:30-35 -p x:10-20 -W t.
is $((vg view t.x:30:35.vg; vg view t.x:10:20.vg) | wc -l) 20 "we can extract a set of targets to separate files"

echo x 30 36 | tr ' ' '\t' >t.bed
echo x 10 21 | tr ' ' '\t' >>t.bed
vg find -x t.xg -E -R t.bed -W q.
is $((vg view q.x:10:20.vg; vg view q.x:30:35.vg) | md5sum | cut -f 1 -d\ ) $((vg view t.x:10:20.vg ; vg view t.x:30:35.vg)| md5sum | cut -f 1 -d\ ) "the same extraction can be made using BED input"

is $(vg find -x t.xg -E -p x:30-35 -p x:10-20 -K 5 | wc -l) 36 "we see the expected number of kmers in the given targets"

rm -f t.xg t.vg t.x:30:35.vg t.x:10:20.vg q.x:30:35.vg q.x:10:20.vg t.bed

vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -x x.xg x.vg
vg gbwt -v small/xy2.vcf.gz -o x.gbwt -x x.vg
is $(vg find -p x -x x.xg -K 16 -H x.gbwt | cut -f 5 | sort | uniq -c  | tail -n 1 | awk '{ print $1 }') 1510 "we find the expected number of kmers with haplotype frequency equal to 2"
rm -f x.vg x.xg x.gbwt

vg convert -gp tiny/tiny.gfa | vg find -x - -n 1 -c 2 | vg convert -f - | vg ids -s - | sort > found1.gfa
vg find -x tiny/tiny.gfa -n 1 -c 2| vg ids -s - | sort > found2.gfa
diff found1.gfa found2.gfa
is $? 0 "GFA i/o for find -n consistent with converting both ways"


# Find nodes that map to the provided ids
vg construct -m 32 -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg gbwt -v small/xy2.vcf.gz -o x.gbwt -x x.vg
vg prune -u -m x.mapping -g x.gbwt -e 1 x.vg > x.unfolded.vg

rm -f expected.gfa
printf "H\tVN:Z:1.1\n" >> expected.gfa
printf "S\t72\tTGGAGTTCTATTATATTCC\n" >> expected.gfa
printf "S\t76\tTGGAGTTCTATTATATTCC\n" >> expected.gfa
printf "S\t97\tTCT\n" >> expected.gfa
printf "S\t98\tTCT\n" >> expected.gfa
printf "S\t82\tTGGAGTTCTATTATATTCC\n" >> expected.gfa
vg find -n 5 -n 23 --mapping x.mapping -x x.unfolded.vg 2> /dev/null | vg view -g - > found.gfa
cmp found.gfa expected.gfa
is $? 0 "find nodes that map to the provided node ids"

rm -f x.vg x.gbwt x.mapping x.unfolded.vg
rm -f expected.gfa found.gfa

#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 52

vg mod -U 10 msgas/hla_v.vg | vg mod -c - > hla_v.vg

#flatten it to one path
printf "gi|568815551:1054737-1055734\t0\t1000\n" > region.bed
vg clip hla_v.vg -b region.bed > clip_flat.vg
vg validate clip_flat.vg
is "$?" 0 "clipped graph is valid"
for i in `vg view clip_flat.vg | grep '^S' | awk '{print $1}'` ; do vg view clip_flat.vg -n $i | grep ^P | grep "gi|568815551:1054737-1055734"; done | wc -l > step_count
vg view clip_flat.vg | grep ^S | wc -l > node_count
diff step_count node_count
is "$?" 0 "every step in clipped graph belongs to reference path"
is "$(vg paths -Ev hla_v.vg -Q "gi|568815551:1054737-1055734" | awk '{ print $2 }')" "$(vg stats -l clip_flat.vg | awk '{ print $2 }')" "clipped graph has same length as ref path"

rm -f region.bed step_count node_count

#clip out one snarl
printf "gi|157734152:29563108-29564082\t90\t92\n" > region.bed
vg clip hla_v.vg -b region.bed > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(vg view clip.vg | grep ^S | wc -l)" "49" "Just one node filtered"

rm -f region.bed clip.vg

# clip out one edge
printf "gi|568815564:1054403-1055400\t150\t153\n" > region.bed
vg clip hla_v.vg -b region.bed > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(vg view clip.vg | grep ^L | wc -l)" "65" "Just one edge filtered"

rm -f region.bed clip.vg

# clip out low coverage node
vg clip hla_v.vg -d 4 -P "gi|568815551:1054737-1055734" > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(vg view clip.vg | grep ^S | wc -l)" "49" "Just one node filtered"

rm -f clip.vg

# clip out out-of-bounds low coverage node
printf "gi|568815551:1054737-1055734\t5\t25\n" > region.bed
vg clip hla_v.vg -b region.bed -d 4 > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
vg view hla_v.vg | sort > hla.gfa
vg view clip.vg | sort > clip.gfa
diff hla.gfa clip.gfa
is "$?" 0 "clipping bad region changes nothing"

rm -f clip.vg hla.gfa clip.gfa

# clip out in-bounds low coverage node
printf "gi|568815551:1054737-1055734\t600\t650\n" > region.bed
vg clip hla_v.vg -b region.bed -d 4 > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is "$(vg view clip.vg | grep ^S | wc -l)" "49" "Just one node filtered"

rm -f region.bed clip.vg

rm -f hla_v.vg

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view - | sort  >  tiny.gfa
cp tiny.gfa tiny-stubs.gfa
printf "S\t0\tA\n" >> tiny-stubs.gfa
printf "L\t0\t+\t1\t+\t0M\n" >> tiny-stubs.gfa
printf "S\t100\tA\n" >> tiny-stubs.gfa
printf "L\t0\t+\t100\t+\t0M\n" >> tiny-stubs.gfa
printf "S\t200\tA\n" >> tiny-stubs.gfa
printf "L\t5\t+\t200\t+\t0M\n" >> tiny-stubs.gfa
printf "S\t300\tA\n" >> tiny-stubs.gfa
printf "L\t200\t+\t300\t+\t0M\n" >> tiny-stubs.gfa
vg clip tiny.gfa -s -P x | sort > tiny-nostubs.gfa
diff tiny.gfa tiny-nostubs.gfa
is "$?" 0 "stub clipping removed all stubs"
diff <(vg clip tiny/tiny.gfa -s -P x | vg paths -v - -E | sort) <(vg clip tiny/tiny.gfaz -s -P x | vg paths -v - -E | sort)
is "$?" 0 "clip preserves path sequences for matching GFA and GFAZ inputs"

printf "x\t5\t25\n" > region.bed
is "$(vg clip tiny-stubs.gfa -s -b region.bed | vg stats -N -)" "17" "region clipping filtered out only 2 / 4 stub nodes"

printf "L\t100\t+\t2\t-\t0M\n" >> tiny-stubs.gfa
printf "L\t15\t+\t13\t-\t0M\n" >> tiny-stubs.gfa
is "$(vg clip tiny-stubs.gfa -sS -P x | vg stats -HT - | sort -nk 2 | awk '{print $2}' | head -1)" "1" "Correct head after path stubbification"
is "$(vg clip tiny-stubs.gfa -sS -P x | vg stats -HT - | sort -nk 2 | awk '{print $2}' | tail -1)" "15" "Correct tail after path stubbification"


rm -f tiny.gfa tiny-stubs.gfa region.bed tiny-nostubs.gfa

# a reference path that starts on the reverse strand: its outer side is the node's right side
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t9\tTTTT\n" > rev-stub.gfa
printf "L\t1\t-\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t9\t+\t1\t-\t0M\nP\tx\t1-,2+,3+\t*\n" >> rev-stub.gfa
is "$(vg clip rev-stub.gfa -S -P x | vg paths -L -x -)" "x" "path stubbification leaves a reverse-strand reference path whole"
is "$(vg clip rev-stub.gfa -S -P x | vg view - | grep -c $'^L\t9\t')" "0" "path stubbification cuts the stub off a reverse-strand reference start"

# a reference path that ends on the reverse strand: its outer side is the node's left side
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t9\tTTTT\n" > rev-end.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t-\t0M\nL\t3\t-\t9\t+\t0M\nP\tx\t1+,2+,3-\t*\n" >> rev-end.gfa
is "$(vg clip rev-end.gfa -S -P x | vg paths -L -x -)" "x" "path stubbification leaves a path ending on the reverse strand whole"
is "$(vg clip rev-end.gfa -S -P x | vg view - | grep -c $'^L\t3\t-\t9\t')" "0" "path stubbification cuts the stub off a reverse-strand reference end"

# reference x is 1 2 3 4.  10 hangs off 1 but also enters 3, 11 hangs off 1, 12-13 hang off 4, and
# 14 runs off 4 and back into 1
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\n" > ends.gfa
printf "S\t10\tA\nS\t11\tC\nS\t12\tG\nS\t13\tT\nS\t14\tA\n" >> ends.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t10\t+\t1\t+\t0M\nL\t10\t+\t3\t+\t0M\n" >> ends.gfa
printf "L\t11\t+\t1\t+\t0M\nL\t4\t+\t12\t+\t0M\nL\t12\t+\t13\t+\t0M\nL\t4\t+\t14\t+\t0M\nL\t14\t+\t1\t+\t0M\n" >> ends.gfa
printf "P\tx\t1+,2+,3+,4+\t*\nP\th1\t10+,1+,2+\t*\nP\th2\t10+,3+\t*\nP\th3\t11+,1+\t*\n" >> ends.gfa
printf "P\th4\t4+,12+,13+\t*\nP\th5\t3+,4+,14+,1+,2+\t*\n" >> ends.gfa
is "$(vg clip ends.gfa -S -P x | vg stats -s - | wc -l)" "4" "path stubbification cuts off everything hanging off the ends"
is "$(vg clip ends.gfa -Sk -P x | vg stats -s - | wc -l)" "1" "connected path stubbification keeps one component"
is "$(vg clip ends.gfa -Sk -P x | vg stats -E -)" "8" "connected path stubbification only cuts the edges that close cycles through the ends"
is "$(vg clip ends.gfa -Sk -P x | vg view - | grep -c $'^L\t10\t+\t1\t')" "0" "connected path stubbification cuts a re-entering stub off the end"

# 20 hangs off 1 with a self loop and no tip: it stays attached, without the loop
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t20\tTTA\n" > loop-stub.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t20\t+\t20\t+\t0M\nL\t20\t+\t1\t+\t0M\nP\tx\t1+,2+\t*\nP\th\t20+,20+,20+,1+,2+\t*\n" >> loop-stub.gfa
is "$(vg clip loop-stub.gfa -Sk -P x | vg stats -s - | wc -l)" "1" "connected path stubbification keeps a looped stub attached"
is "$(vg clip loop-stub.gfa -Sk -P x | vg stats -E -)" "2" "connected path stubbification cuts a tipless stub's loop to give it a tip"
vg clip loop-stub.gfa -sSk -P x > /dev/null 2>&1
is "$?" "1" "connected path stubbification cannot be combined with stub clipping"

# x and y both start at 1, so 20's one edge to them is on both their ends
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t20\tTTA\n" > shared-end.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t4\t+\t0M\nL\t20\t+\t20\t+\t0M\nL\t20\t+\t1\t+\t0M\n" >> shared-end.gfa
printf "P\tx\t1+,2+,3+\t*\nP\ty\t1+,4+\t*\nP\th\t20+,20+,1+,2+\t*\n" >> shared-end.gfa
is "$(vg clip shared-end.gfa -Sk -P x -P y | vg stats -s - | wc -l)" "1" "connected path stubbification keeps a looped stub on an end two paths share"

# the reference comes back through its last node, and y runs inside x
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGG\nS\t4\tTTTT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n" > ref-ends.gfa
printf "L\t3\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tx\t1+,2+,3+,3+,3+,4+\t*\nP\ty\t2+,3+\t*\n" >> ref-ends.gfa
is "$(vg clip ref-ends.gfa -S -P x -P y | vg paths -L -x - | sort | tr '\n' ' ')" "x y " "path stubbification never cuts an edge a reference path walks"

# 14 hangs off both ends of x, but most haplotypes run from 4 into it
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t14\tCATCAT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n" > two-ends.gfa
printf "L\t3\t+\t4\t+\t0M\nL\t4\t+\t14\t+\t0M\nL\t14\t+\t1\t+\t0M\nP\tx\t1+,2+,3+,4+\t*\nP\th1\t3+,4+,14+\t*\n" >> two-ends.gfa
printf "P\th2\t3+,4+,14+\t*\nP\th3\t3+,4+,14+,1+,2+\t*\n" >> two-ends.gfa
is "$(vg clip two-ends.gfa -Sk -P x | vg view - | grep -c $'^L\t4\t+\t14\t')" "1" "connected path stubbification keeps a piece on the end most haplotypes use"

# a ring with a bubble and no tips off 3: only the least used edge, 13 back to 10, needs to go
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t10\tCA\nS\t11\tG\nS\t12\tT\nS\t13\tAC\n" > ring.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t10\t+\t0M\nL\t10\t+\t11\t+\t0M\nL\t10\t+\t12\t+\t0M\n" >> ring.gfa
printf "L\t11\t+\t13\t+\t0M\nL\t12\t+\t13\t+\t0M\nL\t13\t+\t10\t+\t0M\nP\tx\t1+,2+,3+\t*\n" >> ring.gfa
printf "P\th1\t2+,3+,10+,11+,13+\t*\nP\th2\t2+,3+,10+,12+,13+\t*\nP\th3\t3+,10+,11+,13+,10+,12+,13+\t*\n" >> ring.gfa
is "$(vg clip ring.gfa -Sk -P x | vg view - | grep ^L | cut -f 2-5 | LC_ALL=C sort | tr '\t\n' ': ')" "1:+:2:+ 10:+:11:+ 10:+:12:+ 11:+:13:+ 12:+:13:+ 2:+:3:+ 3:+:10:+ " "connected path stubbification makes a tip with the cheapest cut that keeps the piece attached"

# 10-11 hangs off both ends, and once its edge from 4 goes, 10 is a tip: 11's loop can stay
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t10\tCA\nS\t11\tGT\n" > bridge.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t10\t+\t0M\nL\t10\t+\t11\t+\t0M\n" >> bridge.gfa
printf "L\t11\t+\t11\t+\t0M\nL\t11\t+\t1\t+\t0M\nP\tx\t1+,2+,3+,4+\t*\nP\th1\t3+,4+,10+,11+,11+,11+,1+,2+\t*\nP\th2\t11+,11+,11+,1+,2+\t*\n" >> bridge.gfa
is "$(vg clip bridge.gfa -Sk -P x | vg paths -L -x - | grep -c '^h2$')" "1" "connected path stubbification only makes tips for pieces that have none left"

# off 4, 8 leads to a tip at 10 but also to 9, a dead end that only loops: 9 needs a tip too
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t8\tCA\nS\t9\tGT\nS\t10\tAC\n" > dead-end.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t8\t+\t0M\nL\t8\t+\t10\t+\t0M\nL\t8\t+\t9\t+\t0M\n" >> dead-end.gfa
printf "L\t9\t+\t9\t+\t0M\nP\tx\t1+,2+,3+,4+\t*\nP\th1\t1+,2+,3+,4+,8+,10+\t*\nP\th2\t3+,4+,8+,9+,9+,9+\t*\n" >> dead-end.gfa
is "$(vg clip dead-end.gfa -Sk -P x | vg view - | grep ^L | cut -f 2-5 | LC_ALL=C sort | tr '\t\n' ': ')" "1:+:2:+ 2:+:3:+ 3:+:4:+ 4:+:8:+ 8:+:10:+ 8:+:9:+ " "connected path stubbification gives a looped dead end a tip"

# off 4, 18 hangs off the right side of 17 by both its sides: a dead end at 17's right side that only loops
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t16\tCA\nS\t17\tGT\nS\t18\tAC\n" > side-loop.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t17\t+\t0M\nL\t16\t+\t17\t+\t0M\nL\t17\t+\t18\t+\t0M\n" >> side-loop.gfa
printf "L\t17\t+\t18\t-\t0M\nP\tx\t1+,2+,3+,4+\t*\nP\th1\t1+,2+,3+,4+,17+,18+\t*\nP\th2\t16+,17+,18-,17-,4-,3-\t*\n" >> side-loop.gfa
is "$(vg clip side-loop.gfa -Sk -P x | vg stats -E -) $(vg clip side-loop.gfa -Sk -P x | vg stats -s - | wc -l)" "6 1" "connected path stubbification gives a dead end on one side of a node a tip"

# loops chained off the end of x through shared nodes (10 11 12, 11 13 14, 13 15 16): cutting the last one is enough
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\n" > necklace.gfa
for i in 10 11 12 13 14 15 16; do printf "S\t$i\tACGT\n"; done >> necklace.gfa
for e in "1 2" "2 3" "3 4" "4 10" "10 11" "11 12" "12 10" "11 13" "13 14" "14 11" "13 15" "15 16" "16 13"; do set -- $e; printf "L\t$1\t+\t$2\t+\t0M\n"; done >> necklace.gfa
printf "P\tx\t1+,2+,3+,4+\t*\nP\th\t4+,10+,11+,12+,10+,11+,13+,14+,11+,13+,15+,16+,13+,15+\t*\n" >> necklace.gfa
is "$(vg clip necklace.gfa -Sk -P x | vg view - | grep -c $'^L\t15\t+\t16\t') $(vg clip necklace.gfa -Sk -P x | vg stats -E -)" "0 12" "connected path stubbification cuts a chain of loops once, at its far end"

# two pieces on the end of x: 4 has a tip once 4-5 is cut, and 5 then only loops through the end
printf "H\tVN:Z:1.0\nS\t1\tACGTAC\nS\t2\tCCTTGA\nS\t3\tTTGACA\nS\t4\tGATTACA\nS\t5\tCATCAT\n" > two-pieces.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t-\t0M\nL\t3\t+\t5\t+\t0M\nL\t3\t+\t5\t-\t0M\nL\t4\t-\t5\t+\t0M\nP\tx\t1+,2+,3+\t*\n" >> two-pieces.gfa
is "$(vg clip two-pieces.gfa -Sk -P x | vg stats -E -) $(vg clip two-pieces.gfa -Sk -P x | vg clip - -Sk -P x | vg stats -E -)" "5 5" "connected path stubbification leaves nothing for a second run"

rm -f rev-stub.gfa rev-end.gfa ends.gfa loop-stub.gfa shared-end.gfa ref-ends.gfa two-ends.gfa ring.gfa bridge.gfa dead-end.gfa side-loop.gfa necklace.gfa two-pieces.gfa

# clipping A's private node 3 leaves two 8 bp pieces of A
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t5\tACGT\n" > frag.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\n" >> frag.gfa
printf "P\tREF#0#chr1\t1+,2+,4+,5+\t*\nP\tA#1#chr1\t1+,2+,3+,4+,5+\t*\nP\tB#1#chr1\t1+,2+,4+,5+\t*\n" >> frag.gfa
is "$(vg clip frag.gfa -d 2 -P REF | vg paths -L -x - | grep -c '^A#')" "2" "depth clipping keeps both pieces of a chopped path"
is "$(vg clip frag.gfa -d 2 -P REF -m 9 | vg paths -L -x - | grep -c '^A#')" "0" "min fragment length drops pieces that end at a clipped node"

rm -f frag.gfa

# deletion edges: REF is 1-10, of 10 bp nodes
for i in 1 2 3 4 5 6 7 8 9 10; do printf "S\t$i\tAAAAAAAAAA\n"; done > del-ref.gfa
for i in 1 2 3 4 5 6 7 8 9; do printf "L\t$i\t+\t$((i+1))\t+\t0M\n"; done >> del-ref.gfa
printf "P\tREF\t1+,2+,3+,4+,5+,6+,7+,8+,9+,10+\t*\n" >> del-ref.gfa

# 20 is an allele of 2 that H1 and H2 use, and BAD jumps from it to 9
cp del-ref.gfa del-allele.gfa
printf "S\t20\tCCCCCCCCCC\nL\t1\t+\t20\t+\t0M\nL\t20\t+\t3\t+\t0M\nL\t20\t+\t9\t+\t0M\n" >> del-allele.gfa
printf "P\tH1\t1+,20+,3+,4+,5+,6+,7+,8+,9+,10+\t*\nP\tH2\t1+,20+,3+,4+,5+,6+,7+,8+,9+,10+\t*\nP\tBAD\t1+,20+,9+,10+\t*\n" >> del-allele.gfa
is "$(vg clip del-allele.gfa -D 50 -P REF | vg view - | grep $'^L.*\t20\t' | cut -f 2-5 | LC_ALL=C sort | tr '\t\n' ': ')" "1:+:20:+ 20:+:3:+ " "deletion edge clipping cuts the jump, not the allele it leaves from"

# three haplotypes leave REF (15 4 9 8 13 6 7 5) at 4 and run through two bubbles to 5: one cut is enough
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do printf "S\t$i\tACGTTGCAAC\n"; done > del-bubble.gfa
for e in "15 4" "4 9" "9 8" "8 13" "13 6" "6 7" "7 5" "4 3" "4 2" "3 10" "3 1" "2 10" "2 1" "10 11" "10 12" "1 11" "1 12" "11 14" "12 14" "14 5"; do set -- $e; printf "L\t$1\t+\t$2\t+\t0M\n"; done >> del-bubble.gfa
printf "P\tREF\t15+,4+,9+,8+,13+,6+,7+,5+\t*\nP\tH1\t15+,4+,3+,10+,11+,14+,5+\t*\nP\tH2\t15+,4+,2+,1+,12+,14+,5+\t*\nP\tH3\t15+,4+,2+,1+,12+,14+,5+\t*\n" >> del-bubble.gfa
is "$(vg clip del-bubble.gfa -D 20 -c 10 -P REF | vg stats -E -)" "19" "deletion edge clipping makes one cut across a run of bubbles"

# 20-21 is an insertion H1-H10 share between 1 and 2, and X1 jumps through it from 9 to 10: cut X1's own edges
cp del-ref.gfa del-shared.gfa
printf "S\t20\tCCCCC\nS\t21\tGGGGG\nL\t1\t+\t20\t+\t0M\nL\t20\t+\t21\t+\t0M\nL\t21\t+\t2\t+\t0M\nL\t9\t+\t20\t+\t0M\nL\t21\t+\t10\t+\t0M\n" >> del-shared.gfa
for i in 1 2 3 4 5 6 7 8 9 10; do printf "P\tH$i\t1+,20+,21+,2+\t*\n"; done >> del-shared.gfa
printf "P\tX1\t9+,20+,21+,10+\t*\n" >> del-shared.gfa
is "$(vg clip del-shared.gfa -D 50 -P REF | vg paths -L -x - | grep -c '^H.*\[')" "0" "deletion edge clipping leaves an insertion shared with a jump to the haplotypes that stay local"

# BAD jumps from 1 to 10 through a hairpin at 20, whose two sides reach the same reference nodes
cp del-ref.gfa del-hairpin.gfa
printf "S\t20\tCCCCC\nL\t1\t+\t20\t+\t0M\nL\t10\t+\t20\t+\t0M\nL\t20\t+\t20\t-\t0M\nP\tBAD\t1+,20+,20-,10-\t*\n" >> del-hairpin.gfa
is "$(vg clip del-hairpin.gfa -D 50 -P REF | vg stats -E -)" "11" "deletion edge clipping cuts a jump through a hairpin"

# 16 to 2 through 100 is only a deletion once 101's edge into 100 is cut: one run has to find both
for i in $(seq 1 18); do printf "S\t$i\tAAAAAAAAAA\n"; done > del-twice.gfa
for i in $(seq 1 17); do printf "L\t$i\t+\t$((i+1))\t+\t0M\n"; done >> del-twice.gfa
printf "S\t100\tC\nS\t101\tG\nL\t16\t+\t100\t+\t0M\nL\t100\t+\t16\t+\t0M\nL\t100\t+\t2\t+\t0M\nL\t7\t+\t101\t+\t0M\n" >> del-twice.gfa
printf "L\t101\t+\t2\t+\t0M\nL\t101\t+\t100\t+\t0M\nP\tREF\t$(printf '%s+,' $(seq 1 18) | sed 's/,$//')\t*\n" >> del-twice.gfa
is "$(vg clip del-twice.gfa -D 60 -P REF | vg clip - -D 60 -P REF -v 2>&1 >/dev/null | grep -c 'Found deletion edge')" "0" "deletion edge clipping leaves nothing for a second run"

# J jumps from 2 to 9, turning around in a hairpin at 30 that L1 and L2 use. C1 and C2 jump from 4 to 10 through 30
cp del-ref.gfa del-hairpin-d.gfa
printf "S\t30\tCCCCC\nL\t30\t+\t2\t+\t0M\nL\t30\t+\t9\t+\t0M\nL\t30\t-\t30\t+\t0M\nL\t30\t-\t4\t+\t0M\nL\t30\t+\t10\t+\t0M\n" >> del-hairpin-d.gfa
printf "P\tJ\t2-,30-,30+,9+\t*\nP\tL1\t30-,30+\t*\nP\tL2\t30-,30+\t*\nP\tC1\t4-,30+,10+\t*\nP\tC2\t4-,30+,10+\t*\nP\tF\t30+,10+\t*\n" >> del-hairpin-d.gfa
is "$(vg clip del-hairpin-d.gfa -D 50 -d 3 -P REF | vg paths -L -x - | grep -c '^J$')" "0" "deletion edge clipping cuts a jump that turns around in a hairpin too deep to cut"

# 2 to 10 through 100 is the only jump. once its edge from 10 is cut, the one from 2 can go back, but then the
# search from 102 can reach 2 through 100, which no haplotype does: keep the cuts from before putting any back
for n in 1:5 2:12 3:10 4:8 5:7 6:12 7:5 8:15 9:11 10:14 11:8 12:7 13:10 14:7 15:8 16:11 17:12 18:9 100:4 101:10 102:8; do
    printf "S\t${n%:*}\t$(printf "%${n#*:}s" | tr ' ' A)\n"
done > del-putback.gfa
for i in $(seq 1 17); do printf "L\t$i\t+\t$((i+1))\t+\t0M\n"; done >> del-putback.gfa
printf "L\t2\t+\t100\t-\t0M\nL\t10\t+\t100\t+\t0M\nL\t100\t+\t102\t+\t0M\nL\t101\t+\t12\t+\t0M\nL\t102\t+\t17\t+\t0M\n" >> del-putback.gfa
printf "L\t102\t+\t101\t+\t0M\nL\t102\t-\t101\t+\t0M\nP\tREF\t$(printf '%s+,' $(seq 1 18) | sed 's/,$//')\t*\n" >> del-putback.gfa
is "$(vg clip del-putback.gfa -D 60 -c 1000 -P REF | vg stats -E -)" "22" "deletion edge clipping doesn't cut more after putting a cut back"

# a collapsed reference visits 20 twice, 70 bp apart
printf "S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t20\tGGGGGGGGGG\nS\t3\tTTTTTTTTTT\nS\t4\tAAAAAAAAAA\nS\t5\tCCCCCCCCCC\n" > del-collapse.gfa
printf "S\t6\tGGGGGGGGGG\nS\t7\tTTTTTTTTTT\nS\t8\tAAAAAAAAAA\nS\t9\tCCCCCCCCCC\nS\t30\tA\n" >> del-collapse.gfa
for e in "1 2" "2 20" "20 3" "3 4" "4 5" "5 6" "6 7" "7 8" "8 20" "20 9" "1 30" "30 20"; do set -- $e; printf "L\t$1\t+\t$2\t+\t0M\n"; done >> del-collapse.gfa
printf "P\tREF\t1+,2+,20+,3+,4+,5+,6+,7+,8+,20+,9+\t*\nP\tH1\t1+,30+,20+,3+,4+,5+,6+,7+,8+,20+,9+\t*\n" >> del-collapse.gfa
is "$(vg clip del-collapse.gfa -D 50 -P REF | vg stats -E -)" "12" "deletion edge clipping leaves a collapsed reference alone"

rm -f del-ref.gfa del-allele.gfa del-bubble.gfa del-collapse.gfa del-shared.gfa del-hairpin.gfa del-twice.gfa del-hairpin-d.gfa del-putback.gfa

vg clip graphs/snarl-clip.gfa -A 2 -d 2 -P x > sc-A2d2.gfa
is "$(vg find -x sc-A2d2.gfa -n 4 | wc -l)" "0" "Node 4 correctly clipped with snarl length and depth filter"
is "$(vg stats sc-A2d2.gfa -N)" "14" "No other nodes were clipped with snarl length and depth filter"

rm -f sc-A2d2.gfa

vg clip graphs/chain-clip.gfa -A 2 -d1 -P x -g > sc-A2d1g.gfa
is "$(vg stats -E sc-A2d1g.gfa)" "30" "One net edges clipped"
is "$(vg stats -N sc-A2d1g.gfa)" "22" "No other nodes were clipped with net edge filter"

rm -f sc-A2d1g.gfa

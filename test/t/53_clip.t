#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 36

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
is "$(vg clip rev-stub.gfa -S -P x | vg stats -E -)" "2" "path stubbification cuts the stub off a reverse-strand reference end"

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
is "$(vg clip ends.gfa -Sk -P x | vg view - | grep -c -P '^L\t10\t\+\t1\t')" "0" "connected path stubbification cuts a re-entering stub off the end"

# 20 hangs off 1 with a self loop and no tip: it stays attached, without the loop
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t20\tTTA\n" > loop-stub.gfa
printf "L\t1\t+\t2\t+\t0M\nL\t20\t+\t20\t+\t0M\nL\t20\t+\t1\t+\t0M\nP\tx\t1+,2+\t*\nP\th\t20+,20+,20+,1+,2+\t*\n" >> loop-stub.gfa
is "$(vg clip loop-stub.gfa -Sk -P x | vg stats -s - | wc -l)" "1" "connected path stubbification keeps a looped stub attached"
is "$(vg clip loop-stub.gfa -Sk -P x | vg stats -E -)" "2" "connected path stubbification cuts a tipless stub down to a tree"
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
is "$(vg clip two-ends.gfa -Sk -P x | vg view - | grep -c -P '^L\t4\t\+\t14\t')" "1" "connected path stubbification keeps a piece on the end most haplotypes use"

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

rm -f rev-stub.gfa ends.gfa loop-stub.gfa shared-end.gfa ref-ends.gfa two-ends.gfa ring.gfa bridge.gfa

vg clip graphs/snarl-clip.gfa -A 2 -d 2 -P x > sc-A2d2.gfa
is "$(vg find -x sc-A2d2.gfa -n 4 | wc -l)" "0" "Node 4 correctly clipped with snarl length and depth filter"
is "$(vg stats sc-A2d2.gfa -N)" "14" "No other nodes were clipped with snarl length and depth filter"

rm -f sc-A2d2.gfa

vg clip graphs/chain-clip.gfa -A 2 -d1 -P x -g > sc-A2d1g.gfa
is "$(vg stats -E sc-A2d1g.gfa)" "30" "One net edges clipped"
is "$(vg stats -N sc-A2d1g.gfa)" "22" "No other nodes were clipped with net edge filter"

rm -f sc-A2d1g.gfa

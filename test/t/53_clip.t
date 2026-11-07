#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 65

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

# the reference comes back through its last node, and y runs inside x
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGG\nS\t4\tTTTT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n" > ref-ends.gfa
printf "L\t3\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tx\t1+,2+,3+,3+,3+,4+\t*\nP\ty\t2+,3+\t*\n" >> ref-ends.gfa
is "$(vg clip ref-ends.gfa -S -P x -P y | vg paths -L -x - | sort | tr '\n' ' ')" "x y " "path stubbification never cuts an edge a reference path walks"

# -Sk cuts nothing: each selected path gets a private copy of its end base, and so ends in a tip
paths_seqs() { vg paths -F -x $1 | awk '/^>/ {if (n) print n, s; n = $0; s = ""; next} {s = s $0} END {print n, s}' | sort; }
paths_lens() { vg paths -E -x $1 | sort; }
sorted_view() { vg view $1 | sort | md5sum; }
vg clip ends.gfa -Sk -P x > ends-k.gfa
vg validate ends-k.gfa
is $? 0 "connected path stubbification makes a valid graph"
is "$(vg stats -N -E ends-k.gfa | tr '\n' ' ')" "13 14 " "connected path stubbification adds two nodes and two edges for each end of a node longer than 1 bp"
is "$(vg view ends-k.gfa | grep ^S | cut -f 2 | sort -n | head -9 | tr '\n' ' ')" "1 2 3 4 10 11 12 13 14 " "connected path stubbification removes no node"
is "$(paths_seqs ends-k.gfa)" "$(paths_seqs ends.gfa)" "connected path stubbification keeps every path's sequence"
is "$(paths_lens ends-k.gfa)" "$(paths_lens ends.gfa)" "connected path stubbification keeps every path's length"
is "$(vg paths -L -x ends-k.gfa | sort | tr '\n' ' ')" "h1 h2 h3 h4 h5 x " "connected path stubbification chops no path"
is "$(vg clip ends-k.gfa -S -P x | vg stats -E -)" "14" "after connected path stubbification, the reference ends are tips"
is "$(vg clip ends-k.gfa -Sk -P x | sorted_view -)" "$(sorted_view ends-k.gfa)" "connected path stubbification leaves nothing for a second run"
is "$(vg clip ends-k.gfa -Sk -P x | vg view - | grep -c $'^L\t[0-9]*\t+\t14\t+')" "1" "connected path stubbification keeps the loop from the last node back into the first"
# -S then finds nothing to cut, and -s keeps the loop through 14, which is attached at both ends (cutting it is
# left to -D or to whatever made the graph)
is "$(vg clip ends-k.gfa -sS -P x | vg view - | grep -c $'^L\t[0-9]*\t+\t14\t+') $(vg clip ends-k.gfa -sS -P x | vg paths -L -x - | grep -cx h5)" "1 1" "stub removal after connected path stubbification keeps a loop past the reference ends"
vg clip ends.gfa -k -P x > /dev/null 2>&1
is $? 1 "connected path stubbification needs -S"

# x's end nodes are 1 bp: they get a copy each, and nothing is divided
printf "H\tVN:Z:1.0\nS\t1\tA\nS\t2\tCCCC\nS\t3\tT\nS\t8\tGG\nS\t9\tTT\nL\t8\t+\t1\t+\t0M\nL\t1\t+\t2\t+\t0M\n" > short-ends.gfa
printf "L\t2\t+\t3\t+\t0M\nL\t3\t+\t9\t+\t0M\nP\tx\t1+,2+,3+\t*\nP\th\t8+,1+,2+,3+,9+\t*\n" >> short-ends.gfa
vg clip short-ends.gfa -Sk -P x > short-ends-k.gfa
is "$(vg stats -N -E short-ends-k.gfa | tr '\n' ' ')" "7 6 " "connected path stubbification adds one node and one edge for each end of a 1 bp node"
is "$(vg view short-ends-k.gfa | grep $'^P\th\t' | cut -f 3) $(paths_seqs short-ends-k.gfa | tr '\t\n' '  ')" "8+,1+,2+,3+,9+ $(paths_seqs short-ends.gfa | tr '\t\n' '  ')" "connected path stubbification leaves other paths over 1 bp end nodes as they are"

# reverse-strand ends, and a reference path whose end nodes are already tips
is "$(vg clip rev-stub.gfa -Sk -P x | vg view - | grep $'^P\tx\t' | cut -f 3) $(vg clip rev-stub.gfa -Sk -P x | vg paths -F -x - | tail -1)" "11+,1-,2+,3+ $(vg paths -F -x rev-stub.gfa | tail -1)" "connected path stubbification splits a reverse-strand reference start"
is "$(vg clip rev-end.gfa -Sk -P x | vg clip - -S -P x | vg stats -E -) $(vg clip rev-end.gfa -Sk -P x | vg paths -F -x - | tail -1)" "5 $(vg paths -F -x rev-end.gfa | tail -1)" "connected path stubbification splits a reverse-strand reference end"
is "$(vg clip ref-ends.gfa -Sk -P x | sorted_view -)" "$(sorted_view ref-ends.gfa)" "connected path stubbification leaves a path that ends in tips as it is"

# y runs inside x, so its ends are on x's nodes, and its last node is one x visits three times
vg clip ref-ends.gfa -Sk -P x -P y > ref-ends-k.gfa
vg validate ref-ends-k.gfa
is "$? $(paths_seqs ref-ends-k.gfa | tr '\t\n' '  ')" "0 $(paths_seqs ref-ends.gfa | tr '\t\n' '  ')" "connected path stubbification splits a path off the nodes another walks"
is "$(vg clip ref-ends-k.gfa -S -P x -P y | vg stats -E -) $(vg stats -N ref-ends-k.gfa)" "8 8" "connected path stubbification gives a path its own tips inside another"

# x and y both end on 4, and z visits its last node, 3, twice
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGG\nS\t4\tTTTT\nS\t5\tACGT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t3\t+\t0M\n" > shared-end.gfa
printf "L\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\nP\tx\t1+,2+,3+,4+\t*\nP\ty\t3+,4+\t*\nP\tz\t2+,3+,3+\t*\nP\th\t2+,3+,3+,4+,5+\t*\n" >> shared-end.gfa
vg clip shared-end.gfa -Sk -P x -P y -P z > shared-end-k.gfa
vg validate shared-end-k.gfa
is "$? $(paths_seqs shared-end-k.gfa | tr '\t\n' '  ')" "0 $(paths_seqs shared-end.gfa | tr '\t\n' '  ')" "connected path stubbification keeps every path's sequence with shared and revisited end nodes"
is "$(vg clip shared-end-k.gfa -S -P x -P y -P z | vg stats -E -) $(vg view shared-end-k.gfa | grep $'^P\t[xy]\t' | cut -f 3 | sed 's/.*,//' | sort -u | wc -l)" "$(vg stats -E shared-end-k.gfa) 2" "connected path stubbification gives two paths that end on one node a copy each"

# a path of one step on a 1 bp node has no piece to attach a copy to, and on a 2 bp node its copies are joined
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t5\tC\nS\t6\tACG\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n" > one-step.gfa
printf "L\t1\t+\t5\t+\t0M\nL\t5\t+\t3\t+\t0M\nL\t1\t+\t6\t+\t0M\nL\t6\t+\t3\t+\t0M\nP\ta\t2+\t*\nP\tb\t5+\t*\nP\tc\t6+\t*\nP\th\t1+,2+,3+,4+\t*\n" >> one-step.gfa
is "$(vg clip one-step.gfa -Sk -P a | vg stats -N -) $(vg clip one-step.gfa -Sk -P b | vg stats -N -) $(vg clip one-step.gfa -Sk -P c | vg stats -N -)" "9 6 10" "connected path stubbification splits a path of one step unless its node is 1 bp"
is "$(vg clip one-step.gfa -Sk -P a | vg clip - -Sk -P a | sorted_view -) $(vg clip one-step.gfa -Sk -P b | vg clip - -Sk -P b | sorted_view -)" "$(vg clip one-step.gfa -Sk -P a | sorted_view -) $(vg clip one-step.gfa -Sk -P b | sorted_view -)" "connected path stubbification of a path of one step leaves nothing for a second run"

# a circular reference path is left as it is
vg view -Fv ends.gfa > ends.vg 2> /dev/null
vg circularize -p x ends.vg > circular.vg 2> /dev/null
is "$(vg clip circular.vg -Sk -P x | sorted_view -)" "$(sorted_view circular.vg)" "connected path stubbification leaves a circular path as it is"
# circular paths c and d end with forward visits to x's last and first nodes, which are divided
printf "H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCCCT\nS\t3\tGGGG\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t1\t+\t0M\n" > circ-through.gfa
printf "P\tx\t1+,2+\t*\nP\tc\t3+,1+,2+\t*\nP\td\t2+,3+,1+\t*\n" >> circ-through.gfa
vg view -Fv circ-through.gfa 2> /dev/null | vg circularize -p c -p d - > circ-through.vg 2> /dev/null
vg convert -a circ-through.vg > circ-through.hg
vg convert -p circ-through.vg > circ-through.pg
is "$(vg clip circ-through.hg -Sk -P x | vg stats -N -) $(paths_seqs <(vg clip circ-through.hg -Sk -P x) | tr '\t\n' '  ')" "7 $(paths_seqs circ-through.hg | tr '\t\n' '  ')" "connected path stubbification keeps the walks of circular paths through divided end nodes (HashGraph)"
is "$(vg clip circ-through.pg -Sk -P x | vg stats -N -) $(paths_seqs <(vg clip circ-through.pg -Sk -P x) | tr '\t\n' '  ')" "7 $(paths_seqs circ-through.pg | tr '\t\n' '  ')" "connected path stubbification keeps the walks of circular paths through divided end nodes (PackedGraph)"

# past the end of REF, 6 runs back into its first node 1: the top-level chain is a loop, which haplotype sampling
# can't use, until REF ends in tips
printf 'H\tVN:Z:1.1\tRS:Z:REF\n' > ring.gfa
printf 'S\t%s\t%s\n' 1 GATTACAGATTACA 2 CATCATCATCAT 3 GTAGTAGTAGTA 4 TTAGTTAGTTAG 5 GGGATTTGGGATTT 6 ACACACACAC >> ring.gfa
printf 'L\t%s\t+\t%s\t+\t0M\n' 1 2 1 3 2 4 3 4 4 5 5 6 6 1 >> ring.gfa
printf 'W\t%s\t%s\t%s\t0\t*\t%s\n' REF 0 chr1 '>1>2>4>5' H1 1 chr1 '>1>3>4>5>6>1>2' H2 1 chr1 '>4>5>6>1>3>4' >> ring.gfa
vg clip ring.gfa -Sk -P REF > ring-k.gfa
for g in ring ring-k; do
    vg gbwt -G $g.gfa --gbz-format -g $g.gbz -r $g.ri
    vg index -j $g.dist -P REF $g.gbz
done
vg haplotypes --kmer-length 7 --window-length 3 --subchain-length 50 -H ring.hapl ring.gbz 2> /dev/null
is $? 1 "haplotype sampling can't use a graph whose reference ends are on a loop"
vg haplotypes --validate --kmer-length 7 --window-length 3 --subchain-length 50 -H ring-k.hapl ring-k.gbz
is $? 0 "haplotype sampling can use the graph after connected path stubbification"

rm -f rev-stub.gfa rev-end.gfa ends.gfa ref-ends.gfa ends-k.gfa short-ends.gfa short-ends-k.gfa ref-ends-k.gfa
rm -f shared-end.gfa shared-end-k.gfa one-step.gfa ends.vg circular.vg
rm -f circ-through.gfa circ-through.vg circ-through.hg circ-through.pg
rm -f ring.gfa ring-k.gfa ring.gbz ring.ri ring.dist ring.hapl ring-k.gbz ring-k.ri ring-k.dist ring-k.hapl

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

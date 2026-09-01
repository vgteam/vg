#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 30

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
#is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

nodes=$(vg stats -z z.vg | head -1 | cut -f 2)
real_nodes=$(vg view -j z.vg | jq -c '.node[]' | wc -l)
is $nodes $real_nodes "vg stats reports the expected number of nodes"

edges=$(vg stats -z z.vg | tail -1 | cut -f 2)
real_edges=$(vg view -j z.vg | jq -c '.edge[]' | wc -l)
is $edges $real_edges "vg stats reports the expected number of edges"

graph_length=$(vg stats -l z.vg | tail -1 | cut -f 2)
real_length=$(vg view -j z.vg | jq -r '.node[].sequence' | tr -d '\n' | wc -c)
is $graph_length $real_length "vg stats reports the expected graph length"

subgraph_count=$(vg stats -s z.vg | wc -l)
is $subgraph_count 1 "vg stats reports the correct number of subgraphs"

subgraph_length=$(vg stats -s z.vg | head -1 | cut -f 2)
is $subgraph_length $graph_length  "vg stats reports the correct subgraph length"

rm -f z.vg

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >t.vg
is $(vg stats -n 13 -d t.vg | cut -f 2) 38 "distance to head is correct"
is $(vg stats -n 13 -t t.vg | cut -f 2) 11 "distance to tail is correct"
rm -f t.vg

vg construct -m 1000 -r small/x.fa -a -f -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
vg map -x x.xg -g x.gcsa -T small/x-s1337-n100.reads >x.gam
is "$(vg stats -a x.gam x.vg | grep -v "^Speed" | grep -v "^Total time" | md5sum | cut -f 1 -d\ )" "$(md5sum correct/10_vg_stats/15.txt | cut -f 1 -d\ )" "aligned read stats are computed correctly"

is "$(vg stats -z x.vg)" "$(vg stats -z x.xg)" "basic stats agree between graph formats"

is "$(vg stats -a x.gam | grep 'Total alignments')" "Total alignments: 100" "stats can be computed for GAM files without graphs"
rm -f x.vg x.xg x.gcsa x.gam

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view -g - > tiny_names.gfa
printf "P\tref.1\t1+,3+,5+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
printf "P\talt1.1\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
vg view -Fv tiny_names.gfa > tiny_names.vg 
is $(vg stats -O tiny_names.vg | wc -l) 113 "a path overlap description of a test graph has the expected length"
rm -f tiny_names.gfa tiny_names.vg

is "$(vg stats -F graphs/atgc.vg)" "format: VG-Protobuf" "vg stats -F detects format of old protobuf graph"
vg convert graphs/atgc.vg -v > atgc.vg
is "$(vg stats -F atgc.vg)" "format: VG-Protobuf" "vg stats -F detects format of newer protobuf graph"
vg convert graphs/atgc.vg -a > atgc.hg
is "$(vg stats -F atgc.hg)" "format: HashGraph" "vg stats -F detects format of hash graph"
vg convert graphs/atgc.vg -p > atgc.pg
is "$(vg stats -F atgc.pg)" "format: PackedGraph" "vg stats -F detects format of packed graph"
vg index graphs/atgc.vg -x atgc.xg
is "$(vg stats -F atgc.xg)" "format: XG" "vg stats -F detects format of xg graph"
vg convert graphs/atgc.vg -f > atgc.gfa
is "$(vg stats -F atgc.gfa)" "format: GFA" "vg stats -F detects format of GFA graph"
is "$(vg stats -F tiny/tiny.gfaz)" "format: GFA" "vg stats -F classifies GFAZ input on the GFA graph path"
is "$(vg stats -z tiny/tiny.gfaz)" "$(vg stats -z tiny/tiny.gfa)" "basic stats agree between GFA and GFAZ"
rm -f  atgc.vg atgc.hg atgc.pg atgc.xg atgc.gfa
vg autoindex -v tiny/tiny.vcf.gz -r tiny/tiny.fa -w giraffe -p tiny
is "$(vg stats -F tiny.giraffe.gbz)" "format: GBZ" "vg stats -F detects format of GBZ graph"
rm -f tiny.giraffe.gbz tiny.dist tiny.min

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg stats -D - | head -4 | tail -1 > tiny.deg
printf "2\t12\t3\t9\t8\n" > tiny.true.deg
diff tiny.deg tiny.true.deg
is "$?" 0 "vg stats -D found correct degree distribution of tiny graph"
rm -f tiny.def tiny.true.deg


# List self-loops and nondeterministic edge sets.
vg convert -g -a graphs/special_cases.gfa > weird.vg

printf "self-loops\t3\n" > expected.txt
vg stats -L weird.vg > loops.txt
cmp loops.txt expected.txt
is $? 0 "vg stats -L counts self-loops correctly"

rm -f expected.txt
printf "nondeterministic\t1+\t2+\t3+\n" >> expected.txt
printf "nondeterministic\t4-\t2-\t3-\n" >> expected.txt
vg stats -e weird.vg > nondeterministic.txt
cmp nondeterministic.txt expected.txt
is $? 0 "vg stats -e finds the right nondeterministic edge sets"

rm -f weird.vg
rm -f expected.txt loops.txt nondeterministic.txt


# Reference coordinates for snarls (--snarl-sample).  All four graphs below share the
# topology 5(4bp) -> 1(10bp) -> {2(3bp), 3(1bp)} -> 4(7bp) -> 6(5bp), and a 29bp reference
# path visiting 5,1,2,4,6.  The reported coordinates are 0-based and half-open, and cover
# only the interior of each snarl: the boundary nodes are excluded, so the two trivial
# snarls (5,1) and (4,6) come out as empty intervals.

# Forward path: node 5 is at [0,4), 1 at [4,14), 2 at [14,17), 4 at [17,24), 6 at [24,29),
# so the bubble between 1 and 4 has its interior at [14,17).
printf "GRCh38#0#chr1#0\t4\t4\t5\t1\n" > expected.txt
printf "GRCh38#0#chr1#0\t14\t17\t1\t4\n" >> expected.txt
printf "GRCh38#0#chr1#0\t24\t24\t4\t6\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample gives half-open interior coordinates on a forward reference path"

# The same graph with the reference path walking every node in reverse.  The path is the
# reverse complement of the forward case, so each interval mirrors around the 29bp length:
# the bubble interior [14,17) becomes [29-17, 29-14) = [12,15).
printf "GRCh38#0#chr1#0\t5\t5\t4\t6\n" > expected.txt
printf "GRCh38#0#chr1#0\t12\t15\t1\t4\n" >> expected.txt
printf "GRCh38#0#chr1#0\t25\t25\t5\t1\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords_rev.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample gives half-open interior coordinates on a reverse reference path"

# A reference path that is a subpath starting at 1000.  Coordinates must be reported against
# the full contig, and the contig name must not keep the subrange suffix.
printf "GRCh38#0#chr1#0\t1004\t1004\t5\t1\n" > expected.txt
printf "GRCh38#0#chr1#0\t1014\t1017\t1\t4\n" >> expected.txt
printf "GRCh38#0#chr1#0\t1024\t1024\t4\t6\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords_sub.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample offsets coordinates by the subrange of a forward subpath"

# Both at once: a subpath starting at 1000 that also walks every node in reverse.
printf "GRCh38#0#chr1#0\t1005\t1005\t4\t6\n" > expected.txt
printf "GRCh38#0#chr1#0\t1012\t1015\t1\t4\n" >> expected.txt
printf "GRCh38#0#chr1#0\t1025\t1025\t5\t1\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords_sub_rev.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample offsets coordinates by the subrange of a reverse subpath"

# Snarls that only touch the reference at one end can't be bracketed, so they get no
# coordinates at all rather than a made-up empty interval.  These graphs bolt an off-reference
# bubble 6 -> 10 -> {11,12} -> 13 onto the end of the reference path, putting snarls (6,10)
# and (10,13) off the reference.  The answer must not depend on which way the path runs: the
# lone anchor is node 6 either way, and which side of it the bubble attaches to is exactly
# what we cannot tell.
printf ".\t.\t.\t10\t13\n" > expected.txt
printf ".\t.\t.\t6\t10\n" >> expected.txt
printf "GRCh38#0#chr1#0\t4\t4\t5\t1\n" >> expected.txt
printf "GRCh38#0#chr1#0\t14\t17\t1\t4\n" >> expected.txt
printf "GRCh38#0#chr1#0\t24\t24\t4\t6\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords_offref.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample reports no coordinates for snarls that only reach a forward reference at one end"

printf ".\t.\t.\t10\t13\n" > expected.txt
printf ".\t.\t.\t6\t10\n" >> expected.txt
printf "GRCh38#0#chr1#0\t5\t5\t4\t6\n" >> expected.txt
printf "GRCh38#0#chr1#0\t12\t15\t1\t4\n" >> expected.txt
printf "GRCh38#0#chr1#0\t25\t25\t5\t1\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords_offref_rev.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample reports no coordinates for snarls that only reach a reverse reference at one end"

# The same one-anchor situation arises when a contig is split into several subpaths, as in a
# GBZ from minigraph-cactus: snarl (1,4) straddles the junction between the two W-lines, so
# its boundaries land on different path handles and it cannot be bracketed.  The snarls that
# sit wholly within one subpath still get real coordinates, offset onto the full contig.
printf ".\t.\t.\t1\t4\n" > expected.txt
printf "GRCh38#0#chr1#0\t1004\t1004\t5\t1\n" >> expected.txt
printf "GRCh38#0#chr1#0\t1024\t1024\t4\t6\n" >> expected.txt
sort expected.txt > expected.sorted.txt
vg stats -R --snarl-sample GRCh38 graphs/snarl_ref_coords_chopped.gfa | tail -n +2 | cut -f1-4,6 | sort > coords.txt
cmp coords.txt expected.sorted.txt
is "$?" 0 "vg stats --snarl-sample reports no coordinates for a snarl straddling two subpaths of a contig"

rm -f expected.txt expected.sorted.txt coords.txt

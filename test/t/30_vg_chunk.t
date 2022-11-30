#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 30

# Construct a graph with alt paths so we can make a GBWT and a GBZ
vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg x.vg
vg gbwt -x x.vg -v small/x.vcf.gz -o x.haps.gbwt
vg gbwt -x x.vg -E -o x.paths.gbwt
vg gbwt -m x.haps.gbwt x.paths.gbwt -o x.gbwt
vg gbwt -x x.vg x.gbwt --gbz-format -g x.gbz

vg view -a small/x-l100-n1000-s10-e0.01-i0.01.gam > x.gam.json

# sanity check: does passing no options preserve input
is $(vg chunk -x x.xg -p x -c 10| vg stats - -N) 210 "vg chunk with no options preserves nodes"
is $(vg chunk -x x.xg -p x -c 10| vg stats - -E) 291 "vg chunk with no options preserves edges"

# check a small chunk
is $(vg chunk -x x.xg -p x:20-30 -c 0 | vg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == "9"))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk has path going through node 9"
is $(vg chunk -x x.gbz -p x:20-30 -c 0 | vg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == "9"))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk can be found from GBZ"

# check snarl chunking
vg snarls x.xg > x.snarls
is $(vg chunk -x x.xg -S x.snarls -p x:10-20 | vg view - | grep ^S | awk '{print $2}' | sort -n | awk 'BEGIN { ORS = "" } {print}') 6789 "snarl chunk works on simple example"
rm -f x.snarls

# check a small chunk, but using vg input and packed graph output
is $(vg chunk -x x.vg -p x:20-30 -c 0 -O pg | vg convert -v - | vg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == "9"))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk has path going through node 9"

# check no crash when using chunk_size, and filenames deterministic
rm -f _chunk_test*
vg chunk -x x.xg -p x -s 233 -o 50 -b _chunk_test -c 0 -t 2
vg chunk -x x.xg -p x -s 233 -o 50 -b _chunk_test -c 0 -t 1
is $(ls -l _chunk_test*.vg | wc -l) 6 "-s produces correct number of chunks"
rm -f _chunk_test*

#check that gam chunker runs through without crashing
vg gamsort small/x-l100-n1000-s10-e0.01-i0.01.gam -i x.sorted.gam.gai > x.sorted.gam
printf "x\t2\t200\nx\t500\t600\n" > _chunk_test_bed.bed
vg chunk -x x.xg -a x.sorted.gam -g -b _chunk_test -e _chunk_test_bed.bed -E _chunk_test_out.bed -c 0
is $(ls -l _chunk_test*.vg | wc -l) 2 "gam chunker produces correct number of graphs"
is $(ls -l _chunk_test*.gam | wc -l) 2 "gam chunker produces correct number of gams"
is $(grep x _chunk_test_out.bed | wc -l) 2 "gam chunker produces bed with correct number of chunks"
is "$(vg view -aj _chunk_test_0_x_0_199.gam | wc -l)" "$(vg view -aj _chunk_test_0_x_0_199.gam | sort | uniq | wc -l)" "gam chunker emits each matching read at most once"
is "$(vg view -aj _chunk_test_1_x_500_627.gam | wc -l)" "225" "chunk contains the expected number of alignments"
rm -f _chunk_test*

#check that we can chunk by read count
vg chunk -a small/x-l100-n1000-s10-e0.01-i0.01.gam -m 100 -b _chunk_test
is $(ls -l _chunk_test*.gam | wc -l) 10 "simple gam chunker produces correct number of gams"
is "$(vg view -aj _chunk_test000005.gam | wc -l)" "100" "simple chunk contains the expected number of alignments"

#check that id ranges work
is $(vg chunk -x x.xg -r 1:3 -c 0 | vg view - -j | jq .node | grep id |  wc -l) 3 "id chunker produces correct chunk size"
is $(vg chunk -x x.xg -r 1 -c 0 | vg view - -j | jq .node | grep id | wc -l) 1 "id chunker produces correct single chunk"

# Check that traces work on a GBWT
is $(vg chunk -x x.xg -G x.gbwt -r 1:1 -c 2 -T | vg view - -j | jq .node | grep id | wc -l) 5 "id chunker traces correct chunk size"
is "$(vg chunk -x x.xg -r 1:1 -c 2 -T | vg view - -j | jq -c '.path[] | select(.name != "x[0]")' | wc -l)" 0 "chunker extracts no threads from an empty gPBWT"
is "$(vg chunk -x x.xg -G x.haps.gbwt -r 1:1 -c 2 -T | vg view - -j | jq -c '.path[] | select(.name != "x[0]")' | wc -l)" 2 "chunker extracts 2 local threads from a gBWT with 2 locally distinct threads in it"
is "$(vg chunk -x x.xg -G x.gbwt -r 1:1 -c 2 -T | vg view - -j | jq -r '.path[] | select(.name == "thread_0") | .mapping | length')" 3 "chunker can extract a partial haplotype from a GBWT"
is "$(vg chunk -x x.gbz -r 1:1 -c 2 -T | vg view - -j | jq -r '.path[] | select(.name == "thread_0") | .mapping | length')" 3 "chunker can extract a partial haplotype from a GBZ"

vg chunk -x x.xg -G x.gbwt -p x:50-100 -c 0 -T | vg view - | grep ^S | sort > cnodes
vg find -x x.xg -p x:50-100 -c 0 | vg view - | grep ^S | sort > fnodes
diff cnodes fnodes
is "$?" 0 "trace consistent on path coordinates"
rm -f cnodes fnodes

#check that n-chunking works
# We know that it will drop _alt paths so we remake the graph without them for comparison.
vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
mkdir x.chunk
vg chunk -x x.xg -n 5 -b x.chunk/
# vg chunk no longer keeps ranks, so there's no way to stitch the path back together anymore!
# we grep it out of the comparison
is $(cat x.chunk/*vg | vg view -V - | grep -v P 2>/dev/null | sort |  md5sum | cut -f 1 -d\ ) $(vg view x.vg | grep -v P | sort  | md5sum | cut -f 1 -d\ ) "n-chunking works and chunks over the full graph"

rm -rf x.sorted.gam x.sorted.gam.gai _chunk_test_bed.bed _chunk_test* x.chunk
rm -f x.vg x.xg x.gbwt x.gbz x.gam.json filter_chunk*.gam chunks.bed
rm -f chunk_*.annotate.txt

vg construct -r small/xy.fa -v small/xy.vcf.gz > xy.vg
vg construct -r small/xy.fa -v small/xy.vcf.gz -R x > x.vg
vg construct -r small/xy.fa -v small/xy.vcf.gz -R y > y.vg
vg ids -j x.vg y.vg
vg sim -x x.vg -n 50 -a > x.gam
vg sim -x y.vg -n 100 -a > y.gam
cat x.gam y.gam > xy.gam
# test that exploding into components works
vg chunk -x xy.vg -M -b path_chunk -O hg -a xy.gam -g
vg view x.vg | grep "^S" | awk '{print $3}' | sort > x_nodes.txt
vg view y.vg | grep "^S" | awk '{print $3}' | sort > y_nodes.txt
vg convert path_chunk_x.hg -v | vg view - | grep "^S" | awk '{print $3}' | sort > pc_x_nodes.txt
vg convert path_chunk_y.hg -v | vg view - | grep "^S" | awk '{print $3}' | sort > pc_y_nodes.txt
diff x_nodes.txt pc_x_nodes.txt && diff y_nodes.txt pc_y_nodes.txt
is "$?" 0 "path-based components finds subgraphs"
is $(vg view -a path_chunk_x.gam | wc -l) $(vg view -a x.gam | wc -l) "x gam chunk has correct number of reads"
is $(vg view -a path_chunk_y.gam | wc -l) $(vg view -a y.gam | wc -l) "y gam chunk has correct number of reads"
vg chunk -x xy.vg -C -p x -b path_chunk_ind -O hg -a xy.gam -g > /dev/null
vg chunk -x xy.vg -C -p y -b path_chunk_ind -O hg -a xy.gam -g > /dev/null
is $(vg view -a path_chunk_ind_x.gam | wc -l) $(vg view -a x.gam | wc -l) "x gam chunk has correct number of reads with -p path"
is $(vg view -a path_chunk_ind_y.gam | wc -l) $(vg view -a y.gam | wc -l) "y gam chunk has correct number of reads with -p path"
vg paths -v x.vg -E > x_paths.txt
vg paths -v path_chunk_x.hg -E > pc_x_paths.txt
diff pc_x_paths.txt x_paths.txt
is "$?" 0 "path-based component contains correct path length"
vg chunk -x xy.vg -C -b components_chunk
vg view components_chunk_0.vg | grep "^S" | awk '{print $3}' > comp_0_nodes.txt
vg view components_chunk_1.vg | grep "^S" | awk '{print $3}' > comp_1_nodes.txt
cat comp_0_nodes.txt comp_1_nodes.txt | sort > comp_nodes.txt
cat x_nodes.txt y_nodes.txt | sort > nodes.txt
diff comp_nodes.txt nodes.txt
is "$?" 0 "components finds subgraphs"

rm -f xy.vg x.vg y.vg x_nodes.txt y_nodes.txt convert path_chunk_x.hg  convert path_chunk_y.hg pc_x_nodes.txt pc_y_nodes.txt x_paths.txt pc_x_paths.txt components_chunk_0.vg components_chunk_1.vg comp_0_nodes.txt comp_1_nodes.txt comp_nodes.txt nodes.txt x.gam y.gam xy.gam path_chunk_x.gam path_chunk_y.gam




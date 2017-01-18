#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 7

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg  x.vg
vg sim -x x.xg -l 100 -n 1000 -s 0 -e 0.01 -i 0.001 -a > x.gam
vg view -a x.gam > x.gam.json

# sanity check: does passing no options preserve input
is $(vg chunk -x x.xg -p x | vg stats - -N) 210 "vg chunk with no options preserves nodes"
is $(vg chunk -x x.xg -p x | vg stats - -E) 291 "vg chunk with no options preserves edges"

# check a small chunk
is $(vg chunk -x x.xg -p x:20-30 | vg view - -j | jq -c '.path[0].mapping[].position' | jq 'select ((.node_id == 9))' | grep node | sed s/,// | sort | uniq | wc -l) 1 "chunk has path going through node 9"

# check no crash when using chunk_size, and filenames deterministic
rm -f _chunk_test*
vg chunk -x x.xg -p x -s 233 -o 50 -b _chunk_test -t 2
vg chunk -x x.xg -p x -s 233 -o 50 -b _chunk_test -t 1
is $(ls -l _chunk_test*.vg | wc -l) 6 "-s produces correct number of chunks"
rm -f _chunk_test*

#check that gam chunker runs through without crashing
vg index -N x.gam -d x.gam.index
printf "x\t2\t200\nx\t500\t600\n" > _chunk_test_bed.bed
vg chunk -x x.xg -a x.gam.index -g -b _chunk_test -r _chunk_test_bed.bed -R _chunk_test_out.bed
is $(ls -l _chunk_test*.vg | wc -l) 2 "gam chunker produces correct number of graphs"
is $(ls -l _chunk_test*.gam | wc -l) 2 "gam chunker produces correct number of gams"
is $(grep x _chunk_test_out.bed | wc -l) 2 "gam chunker prodcues bed with correct number of chunks"

rm -rf x.gam.index _chunk_test_bed.bed _chunk_test*
rm -f x.vg x.xg x.gam x.gam.json filter_chunk*.gam chunks.bed

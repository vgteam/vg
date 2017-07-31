#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 49

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg

vg index -s -d x.idx x.vg
is $? 0 "indexing nodes and edges of graph"
rm -rf x.idx

vg index -x x.xg x.vg
is $? 0 "building an xg index of the graph"
rm -f x.xg

vg mod -D x.vg >y.vg
cp y.vg z.vg
vg ids -j y.vg z.vg
vg index -x x.xg y.vg z.vg
is $? 0 "building an xg index using multiple input files"
rm -rf x.xg y.vg z.vg

vg index -s -d x.idx x.vg bogus123.vg
is $? 134 "fail with nonexistent file"
rm -rf x.idx

vg index -k 11 -d x.idx x.vg
is $? 0 "indexing 11mers"

vg index -C -d x.idx
is $? 0 "index compaction"
rm -rf x.idx

vg index -g x.gcsa -k 16 x.vg
is $? 0 "building a GCSA2 index"
rm -f x.gcsa x.gcsa.lcp

vg kmers -k 16 -gB x.vg >x.graph
vg index -i x.graph -g x.gcsa
is $? 0 "a prebuilt deBruijn graph in GCSA2 format may be used"
rm -f x.gcsa x.gcsa.lcp x.graph

vg index -s -k 11 -d x.idx x.vg
num_records=$(vg index -D -d x.idx | wc -l)
is $? 0 "dumping graph index"
is $num_records 3208 "correct number of records in graph index"

vg index -x x.xg -g x.gcsa -k 11 x.vg
vg map -T <(vg sim -s 1337 -n 100 -x x.xg) -d x > x1337.gam
vg index -a x1337.gam -d x.vg.aln
is $(vg index -D -d x.vg.aln | wc -l) 101 "index can store alignments"
is $(vg index -A -d x.vg.aln | vg view -a - | wc -l) 100 "index can dump alignments"

# repeat with an unmapped read (sequence from phiX)
rm -rf x.vg.aln
(vg map -s "CTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTC" -d x -P 0.9; vg map -s "CTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTC" -d x ) | vg index -a - -d x.vg.aln
is $(vg index -A -d x.vg.aln | vg view -a - | wc -l) 2 "alignment index stores unmapped reads"

# load the same data again & see that the records are duplicated.
# that's not really a wise use-case, but it tests that we don't
# overwrite anything in an existing alignment index, by virtue
# of keying each entry with a unqiue nonce.
rm -rf x.vg.aln
vg index -a x1337.gam -d x.vg.aln
vg index -a x1337.gam -d x.vg.aln
is $(vg index -D -d x.vg.aln | wc -l) 201 "alignment index can be loaded using sequential invocations; next_nonce persistence"

vg map -T <(vg sim -s 1337 -n 100 -x x.xg) -d x | vg index -m - -d x.vg.map
is $(vg index -D -d x.vg.map | wc -l) 1477 "index stores all mappings"

rm -rf x.idx x.vg.map x.vg.aln x1337.gam

vg construct -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg -v small/x.vcf.gz x.vg
is $? 0 "building an xg index containing a gPBWT"

xg -i x.xg -x > part.vg
is "$(cat x.vg part.vg | vg view -j - | jq '.path[].name' | grep '_thread' | wc -l)" 4 "the gPBWT can be queried for two threads for each haplotype"

is $(vg find -x x.xg -q _thread_1_x_0 | vg paths -L - | wc -l) 1 "a specific thread may be pulled from the graph by name"

rm -f x.vg x.xg part.vg x.gcsa

vg construct -r small/xy.fa -v small/xy.vcf.gz -a >xy.vg
vg index -x xy.xg -v small/xy.vcf.gz xy.vg
is $(vg find -x xy.xg -t | vg paths -L - | wc -l) 4 "a thread is stored per haplotype, sample, and reference sequence"
is $(vg find -x xy.xg -q _thread_1_y | vg paths -L - | wc -l) 2 "we have the expected number of threads per chromosome"
rm -f xy.vg xy.xg

vg construct -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg -v small/x.vcf.gz -H haps.bin x.vg
is $(du -b haps.bin | cut -f 1) 2264 "threads may be exported to binary for use in GBWT construction"

rm -f x.vg x.xg part.vg x.gcsa haps.bin

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg concat x.vg y.vg z.vg >q.vg

vg index -s -d q.idx q.vg
vg index -s -d x.idx x.vg

single=$(vg index -D -d x.idx | wc -l)
triple=$(vg index -D -d q.idx | wc -l)

# subtract four for metadata lines about paths that aren't duplicated in the merged
is $triple $(echo "$single * 3 - 4" | bc) "storage of multiple graphs in an index succeeds"

vg ids -j x.vg q.vg
vg index -k 2 -g qx.vg.gcsa q.vg x.vg
is $? 0 "building a GCSA2 index of two graphs"

rm x.vg y.vg z.vg q.vg
rm -rf x.idx q.vg.index qx.vg.gcsa qx.vg.gcsa.lcp

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s -d x.idx x.vg
is $(vg index -D -d x.idx | grep +g | grep +p | wc -l) $(vg view x.vg | grep ^P | cut -f 3 | grep -o "[0-9]\+" | wc -l) "correct number of elements in path index"
is $(vg index -D -d x.idx | grep +path_id | wc -l) 1 "path id recorded"
is $(vg index -D -d x.idx | grep +path_name | wc -l) 1 "path name recorded"
rm -rf x.idx x.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -v small/x.vcf.gz -r small/x.fa | vg view - | sed s/x/y/ | vg view -Fv - >y.vg
vg ids -j x.vg y.vg

vg index -s -d q.idx x.vg y.vg
## Seems this can sometimes be -420... is that a permissible value?
is $(vg index -L -d q.idx | tail -1 | awk '{ print $3 }' ) 420 "end of the second path found correctly"

rm -rf q.idx x.vg y.vg

# Now test backward nodes
vg index -s -d r.idx reversing/reversing_x.vg
is $? 0 "can index backward nodes"
ls r.idx/*.sst
is $? 0 "backward node index contains data"

vg index -k 16 -d r.idx reversing/reversing_x.vg
is $? 0 "can index kmers for backward nodes"

vg index -x r.xg -g r.gcsa -k 11 reversing/reversing_x.vg

is $(vg index -D -d r.idx | grep "TATTAGCCATGTGACT" | wc -l) 1 "kmers crossing reversing edges are in index"

is $(vg index -D -d r.idx | grep '"from": 55' | grep '"from_start": true' | wc -l) 2 "from_start edges in index"

is $(vg index -D -d r.idx | grep '"to": 55' | grep '"to_end": true' | wc -l) 2 "to_end edges in index"

vg map -T <(vg sim -s 1338 -n 100 -x r.xg) -d r | vg index -a - -d r.aln.idx
is $(vg index -D -d r.aln.idx | wc -l) 101 "index can store alignments to backward nodes"

rm -rf r.idx r.aln.idx r.xg r.gcsa

vg index -k 16 -s -d c.idx cyclic/all.vg
is $? 0 "index can store a cyclic graph"

vg index -x c.xg -g c.gcsa -k 11 cyclic/all.vg
vg map -T <(vg sim -s 1337 -n 100 -x c.xg) -d c | vg index -a - -d all.vg.aln
is $(vg index -D -d all.vg.aln | wc -l) 101 "index can store alignments to cyclic graphs"

rm -rf c.idx all.vg.aln c.xg c.gcsa

is $(vg index -g x.gcsa -k 16 -V <(vg view -Fv cyclic/two_node.gfa) 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with heads and tails"

is $(vg index -g x.gcsa -k 16 -V cyclic/no_heads.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with no heads or tails"

is $(vg index -g x.gcsa -k 16 -V cyclic/self_loops.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with self loops"

is $(vg index -g x.gcsa -k 16 -V cyclic/all.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on general cyclic graphs"

is $(vg index -g x.gcsa -k 16 -V -F cyclic/no_heads.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 forward-only indexing works on cyclic graphs with no heads or tails"

rm -f x.gcsa x.gcsa.lcp

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg index -g t.gcsa -k 16 -V - 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing of a tiny graph works"

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg index -g t.gcsa -k 16 -V - 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 forward-only indexing of a tiny graph works"

is $(vg construct -r tiny/tiny.fa | vg index -g t.gcsa -k 16 -V - 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing succeeds on a single-node graph"

is $(vg construct -r tiny/tiny.fa | vg index -g t.gcsa -k 16 -V -F - 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 forward-only indexing succeeds on a single-node graph"

is $(vg index -g t.gcsa reversing/cactus.vg -k 16 -V 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing succeeds on graph with heads but no tails"

is $(vg index -g t.gcsa reversing/cactus.vg -k 16 -V -F 2>&1 | grep 'Index verification complete' | wc -l) 0 "GCSA2 forward-only indexing fails due to impossibility on graph with heads but no tails"

vg construct -r ins_and_del/ins_and_del.fa -v ins_and_del/ins_and_del.vcf.gz -a >ins_and_del.vg
is $(vg index -x ins_and_del.vg.xg -v ins_and_del/ins_and_del.vcf.gz ins_and_del.vg 2>&1 | wc -l) 0 "indexing with allele paths handles combination insert-and-deletes"

vg construct -m 1025 -r 1mb1kgp/z.fa > big.vg

is "$(vg index -g big.gcsa big.vg -k 16 2>&1 | head -n10 | grep "vg mod -X" | wc -l)" "1" "a useful error message is produced when nodes are too large"

rm -f big.vg

rm -f t.gcsa
rm -f x.vg

rm -f r.gcsa.lcp c.gcsa.lcp t.gcsa.lcp ins_and_del.vg ins_and_del.vg.xg

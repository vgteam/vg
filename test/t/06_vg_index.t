#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 63

# Single graph without haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg

vg index -x x.xg x.vg
is $? 0 "building an XG index of a graph"

vg index -g x.gcsa x.vg
is $? 0 "building a GCSA index of a graph"

vg index -x x2.xg -g x2.gcsa x.vg
is $? 0 "building both indexes at once"

cmp x.xg x2.xg && cmp x.gcsa x2.gcsa && cmp x.gcsa.lcp x2.gcsa.lcp
is $? 0 "the indexes are identical when built one at a time and together"

vg index -x x.xg -g x3.gcsa
is $? 0 "building GCSA from XG"

cmp x.gcsa x3.gcsa && cmp x.gcsa.lcp x3.gcsa.lcp
is $? 0 "the GCSA indexes are identical when built from vg and from xg"

rm -f x.vg
rm -f x.xg x.gcsa x.gcsa.lcp
rm -f x2.xg x2.gcsa x2.gcsa.lcp
rm -f x3.gcsa x3.gcsa.lcp


# Single graph with haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg

vg index -G x.gbwt -v small/x.vcf.gz x.vg
is $? 0 "building a GBWT index of a graph with haplotypes"

vg index -x x.xg x.vg
is $? 0 "building an XG index of a graph with haplotypes"

is $(vg paths -x x.xg -L | wc -l) 1 "xg index does not contain alt paths by default"

vg index -x x-ap.xg x.vg -L
is $? 0 "building an XG index of a graph with haplotypes and alt paths included"

is $(vg paths -x x-ap.xg -L | wc -l) $(vg paths -v x.vg -L | wc -l) "xg index does contains alt paths with index -L"

vg index -g x.gcsa x.vg
is $? 0 "building a GCSA index of a graph with haplotypes"

vg index -x x2.xg -G x2.gbwt -v small/x.vcf.gz -g x2.gcsa x.vg
is $? 0 "building all indexes at once"

cmp x.xg x2.xg && cmp x.gbwt x2.gbwt && cmp x.gcsa x2.gcsa && cmp x.gcsa.lcp x2.gcsa.lcp
is $? 0 "the indexes are identical"

vg index -x x2-ap.xg -G x2-ap.gbwt -v small/x.vcf.gz -g x2-ap.gcsa x.vg -L
is $? 0 "building all indexes at once, while leaving alt paths in xg"

cmp x.gbwt x2-ap.gbwt && cmp x.gcsa x2-ap.gcsa && cmp x.gcsa.lcp x2-ap.gcsa.lcp
is $? 0 "the indexes are identical with -L"

is $(vg paths -x x2-ap.xg -L | wc -l) $(vg paths -v x.vg -L | wc -l) "xg index does contains alt paths with index -L all at once"

# Build the same GBWT indirectly from a VCF parse
vg index -v small/x.vcf.gz -e parse x.vg
is $? 0 "storing a VCF parse for a graph with haplotypes"

../deps/gbwt/build_gbwt -p -r parse_x > /dev/null 2> /dev/null
is $? 0 "building a GBWT index from the VCF parse"

# Dump the tagged stream and compare the indexes built with build_gbwt and vg index
vg view --extract-tag GBWT x.gbwt > x.bare.gbwt
cmp parse_x.gbwt x.bare.gbwt
is $? 0 "the indexes are identical"

# Exclude a sample from the GBWT index
vg index -G empty.gbwt -v small/x.vcf.gz --exclude 1 x.vg
is $? 0 "samples can be excluded from haplotype indexing"
is $(vg gbwt -c empty.gbwt) 0 "excluded samples were not included in the GBWT index"

# Make GBWT from GAM
vg paths -v x.vg -X -Q _alt > x-alts.gam
vg index x.vg -M x-alts.gam -G x-gam.gbwt
# Make GBWT from GAF
vg convert x.vg -G x-alts.gam > x-alts.gaf
vg index x.vg -F x-alts.gaf -G x-gaf.gbwt
cmp x-gaf.gbwt x-gam.gbwt
is $? 0 "GBWT from GAF same as from GAM"

rm -f x.vg
rm -f x.xg x-ap.xg x.gbwtx.gcsa x.gcsa.lcp
rm -f x2.xg x2.gbwt x2.gcsa x2.gcsa.lcp
rm -f x2-ap.xg x2-ap.gbwt x2-ap.gcsa x2-ap.gcsa.lcp
rm -f parse_x parse_x_0_1 parse_x.gbwt x.bare.gbwt
rm -f empty.gbwt
rm -f x-alts.gam x-alts.gaf x-gam.gbwt x-gaf.gbwt


# Subregion graph with haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a --region x:100-200 > x.part.vg

vg index -x x.part.xg -G x.part.gbwt --region x:100-200 -v small/x.vcf.gz x.part.vg 2>log.txt
is $? 0 "building GBWT index for a regional graph"

is "$(cat log.txt | wc -c)" "0" "no warnings about missing variants produced"

rm -f x.part.vg x.part.xg x.part.gbwt log.txt


# Multiple graphs without haplotypes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C > y.vg 2> /dev/null
vg ids -j x.vg y.vg

vg index -x xy.xg x.vg y.vg
is $? 0 "building an XG index of multiple graphs"

vg index -g xy.gcsa -k 2 x.vg y.vg
is $? 0 "building a GCSA index of multiple graphs"

vg index -x xy2.xg -g xy2.gcsa -k 2 x.vg y.vg
is $? 0 "building both indexes at once"

cmp xy.xg xy2.xg && cmp xy.gcsa xy2.gcsa && cmp xy.gcsa.lcp xy2.gcsa.lcp
is $? 0 "the indexes are identical"

rm -f x.vg y.vg
rm -f xy.xg xy.gcsa xy.gcsa.lcp
rm -f xy2.xg xy2.gcsa xy2.gcsa.lcp


# Multiple graphs with haplotypes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg

vg index -G x.gbwt -v small/xy2.vcf.gz x.vg && vg index -G y.gbwt -v small/xy2.vcf.gz y.vg && vg gbwt -m -f -o xy.gbwt x.gbwt y.gbwt
is $? 0 "building a GBWT index of multiple graphs with haplotypes"

vg index -x xy.xg x.vg y.vg
is $? 0 "building an XG index of multiple graphs with haplotypes"

vg index -g xy.gcsa -k 2 x.vg y.vg
is $? 0 "building a GCSA index of multiple graphs with haplotypes"

vg index -x xy2.xg -g xy2.gcsa -k 2 -G xy2.gbwt -v small/xy2.vcf.gz x.vg y.vg
is $? 0 "building all three indexes at once"

cmp xy.xg xy2.xg && cmp xy.gcsa xy2.gcsa && cmp xy.gcsa.lcp xy2.gcsa.lcp && cmp xy.gbwt xy2.gbwt
is $? 0 "the indexes are identical"

# Build the same GBWT indirectly from a VCF parse
vg index -v small/xy2.vcf.gz -e parse x.vg && vg index -v small/xy2.vcf.gz -e parse y.vg
is $? 0 "storing a VCF parse for multiple graphs with haplotypes"

../deps/gbwt/build_gbwt -p -r -o parse_xy parse_x parse_y > /dev/null 2> /dev/null
is $? 0 "building a GBWT index from the VCF parses"

# Dump the tagged stream and compare the indexes built with build_gbwt and vg index
vg view --extract-tag GBWT xy.gbwt > xy.bare.gbwt
cmp parse_xy.gbwt xy.bare.gbwt
is $? 0 "the indexes are identical"

rm -f x.vg y.vg
rm -f x.gbwt y.gbwt
rm -f xy.xg xy.gbwt xy.gcsa xy.gcsa.lcp
rm -f xy2.xg xy2.gbwt xy2.gcsa xy2.gcsa.lcp
rm -f parse_x parse_x_0_1 parse_y parse_y_0_1 parse_xy.gbwt xy.bare.gbwt


# GBWT construction options
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null

vg index -G x_ref.gbwt -T x.vg
is $? 0 "GBWT can be built for paths"

vg index -G x_both.gbwt -T -v small/xy2.vcf.gz x.vg
is $? 0 "GBWT can be built for both paths and haplotypes"

rm -f x_ref.gbwt x_both.gbwt

# We do not test GBWT construction parameters (-B, -u, -n) because they matter only for large inputs.
# We do not test chromosome-length path generation (-P, -o) for the same reason.


# Other tests
vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg bogus123.vg 2>/dev/null
is $? 1 "fail with nonexistent file"
rm -rf x.idx

vg kmers -k 16 -gB x.vg >x.graph
vg index -i x.graph -g x.gcsa
is $? 0 "a prebuilt deBruijn graph in GCSA2 format may be used"
rm -f x.gcsa x.gcsa.lcp x.graph

vg index -x x.xg -g x.gcsa -k 11 x.vg
vg map -T <(vg sim -n 100 -x x.xg) -d x > x1337.gam
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

vg map -T small/x-s1337-n100-v2.reads -d x | vg index -m - -d x.vg.map
is $(vg index -D -d x.vg.map | wc -l) 1477 "index stores all mappings"

# Compare our indexes against gamsort's
vg gamsort -i x1337.sorted.gam.gai2 x1337.gam > x1337.sorted.gam
vg index --index-sorted-gam x1337.sorted.gam

is "$(md5sum <x1337.sorted.gam.gai)" "$(md5sum <x1337.sorted.gam.gai2)" "vg index and vg gamsort produce identical sorted GAM indexes"

rm -rf x.idx x.vg.map x.vg.aln x1337.gam x1337.sorted.gam.gai2 x1337.sorted.gam.gai x1337.sorted.gam


vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg -v small/x.vcf.gz -H haps.bin x.vg
is $(du -b haps.bin | cut -f 1) 329 "threads may be exported to binary for use in GBWT construction"

rm -f x.vg x.xg part.vg x.gcsa haps.bin x.gbwt

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg concat x.vg y.vg z.vg >q.vg

vg index -x q.xg q.vg

is $? 0 "storage of multiple graphs in an index succeeds"

rm x.vg y.vg z.vg q.vg
rm -rf x.idx q.xg

# Now test backward nodes
vg index -x r.xg reversing/reversing_x.vg
is $? 0 "can index backward nodes"

vg index -k 11 -g r.gcsa reversing/reversing_x.vg
is $? 0 "can index kmers for backward nodes"

vg map -T <(vg sim -n 100 -x r.xg) -d r | vg index -a - -d r.aln.idx
is $(vg index -D -d r.aln.idx | wc -l) 101 "index can store alignments to backward nodes"

rm -rf r.aln.idx r.xg r.gcsa

vg index -x c.xg -g c.gcsa -k 11 cyclic/all.vg
vg map -T <(vg sim -n 100 -x c.xg) -d c | vg index -a - -d all.vg.aln
is $(vg index -D -d all.vg.aln | wc -l) 101 "index can store alignments to cyclic graphs"

rm -rf all.vg.aln c.xg c.gcsa

is $(vg index -g x.gcsa -k 16 -V <(vg view -Fv cyclic/two_node.gfa) 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with heads and tails"

is $(vg index -g x.gcsa -k 16 -V cyclic/no_heads.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with no heads or tails"

is $(vg index -g x.gcsa -k 16 -V cyclic/self_loops.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with self loops"

is $(vg index -g x.gcsa -k 16 -V cyclic/all.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on general cyclic graphs"

rm -f x.gcsa x.gcsa.lcp

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg index -g t.gcsa -k 16 -V - 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing of a tiny graph works"

is $(vg construct -r tiny/tiny.fa | vg index -g t.gcsa -k 16 -V - 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing succeeds on a single-node graph"

is $(vg index -g t.gcsa reversing/cactus.vg -k 16 -V 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing succeeds on graph with heads but no tails"

vg construct -m 1025 -r 1mb1kgp/z.fa > big.vg

is $(vg index -g big.gcsa big.vg -k 16 2>&1 | head -n10 | grep 'Found kmer with offset' | wc -l) 1 "a useful error message is produced when nodes are too large"

rm -f big.vg

rm -f t.gcsa
rm -f x.vg

rm -f r.gcsa.lcp c.gcsa.lcp t.gcsa.lcp

# Test distance index 
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg snarls -t x.vg > snarls.pb

vg index -s snarls.pb -j distIndex -w 100 x.vg
is $? 0 "building a distance index of a graph"

vg index -s snarls.pb -j distIndex x.vg
is $? 0 "building a distance index of a graph without maximum index"

rm -f x.vg distIndex snarls.pb

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order

plan tests 45

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

vg index -g x3.gcsa x.xg
is $? 0 "building GCSA from XG"

cmp x.gcsa x3.gcsa && cmp x.gcsa.lcp x3.gcsa.lcp
is $? 0 "the GCSA indexes are identical when built from vg and from xg"

rm -f x.vg
rm -f x.xg x.gcsa x.gcsa.lcp
rm -f x2.xg x2.gcsa x2.gcsa.lcp
rm -f x3.gcsa x3.gcsa.lcp


# Single graph with haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg

vg index -x x.xg x.vg
is $? 0 "building an XG index of a graph with haplotypes"

is $(vg paths -x x.xg -L | wc -l) 1 "xg index does not contain alt paths by default"

vg index -x x-ap.xg x.vg -L
is $? 0 "building an XG index of a graph with haplotypes and alt paths included"

is $(vg paths -x x-ap.xg -L | wc -l) $(vg paths -v x.vg -L | wc -l) "xg index does contains alt paths with index -L"

vg index -g x.gcsa x.vg
is $? 0 "building a GCSA index of a graph with haplotypes"

vg index -x x2.xg -g x2.gcsa x.vg
is $? 0 "building all indexes at once"

cmp x.xg x2.xg && cmp x.gcsa x2.gcsa && cmp x.gcsa.lcp x2.gcsa.lcp
is $? 0 "the indexes are identical"

vg index -x x2-ap.xg -g x2-ap.gcsa x.vg -L
is $? 0 "building all indexes at once, while leaving alt paths in xg"

cmp x.gcsa x2-ap.gcsa && cmp x.gcsa.lcp x2-ap.gcsa.lcp
is $? 0 "the indexes are identical with -L"

is $(vg paths -x x2-ap.xg -L | wc -l) $(vg paths -v x.vg -L | wc -l) "xg index does contains alt paths with index -L all at once"

rm -f x.vg
rm -f x.xg x-ap.xg x.gcsa x.gcsa.lcp
rm -f x2.xg x2.gcsa x2.gcsa.lcp
rm -f x2-ap.xg x2-ap.gcsa x2-ap.gcsa.lcp
rm -f x-alts.gam x-alts.gaf


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

vg index -x xy.xg x.vg y.vg
is $? 0 "building an XG index of multiple graphs with haplotypes"

vg index -g xy.gcsa -k 2 x.vg y.vg
is $? 0 "building a GCSA index of multiple graphs with haplotypes"

vg index -x xy2.xg -g xy2.gcsa -k 2 x.vg y.vg
is $? 0 "building XG and GCSA indexes at once"

vg index -x xy-alt.xg -L x.vg y.vg
is $? 0 "building an XG index with alt paths"

cmp xy.xg xy2.xg && cmp xy.gcsa xy2.gcsa && cmp xy.gcsa.lcp xy2.gcsa.lcp
is $? 0 "the indexes are identical"

rm -f x.vg y.vg
rm -f xy.xg xy.gcsa xy.gcsa.lcp
rm -f xy2.xg xy2.gcsa xy2.gcsa.lcp
rm -f xy-alt.xg

# Other tests
vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg bogus123.vg 2>/dev/null
is $? 1 "fail with nonexistent file"

vg kmers -k 16 -gB x.vg >x.graph
vg index -i x.graph -g x.gcsa
is $? 0 "a prebuilt deBruijn graph in GCSA2 format may be used"
rm -f x.gcsa x.gcsa.lcp x.graph

vg index -x x.xg -g x.gcsa -k 11 x.vg
vg map -T <(vg sim -n 100 -x x.xg) -d x > x1337.gam

# Compare our indexes against gamsort's
vg gamsort -i x1337.sorted.gam.gai2 x1337.gam > x1337.sorted.gam
vg index --index-sorted-gam x1337.sorted.gam

is "$(md5sum <x1337.sorted.gam.gai)" "$(md5sum <x1337.sorted.gam.gai2)" "vg index and vg gamsort produce identical sorted GAM indexes"

rm -rf x1337.gam x1337.sorted.gam.gai2 x1337.sorted.gam.gai x1337.sorted.gam


vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg concat x.vg y.vg z.vg >q.vg

vg index -x q.xg q.vg

is $? 0 "storage of multiple graphs in an index succeeds"

rm x.vg y.vg z.vg q.vg
rm -rf q.xg

# Now test backward nodes
vg index -x r.xg reversing/reversing_x.vg
is $? 0 "can index backward nodes"

vg index -k 11 -g r.gcsa reversing/reversing_x.vg
is $? 0 "can index kmers for backward nodes"

rm -rf r.xg r.gcsa


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
rm -f x.vg x.xg

rm -f r.gcsa.lcp c.gcsa.lcp t.gcsa.lcp


# Test distance index
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg

vg snarls -T x.vg > snarls.pb
is $? 0 "snarl finding with trivial snarls"

vg index -j distIndex x.vg
is $? 0 "building a distance index of a graph"

vg convert -v x.vg | vg index -j distIndex  -
is $? 0 "building a distance index of a piped vg graph"

vg convert -p x.vg | vg index -j distIndex -
is $? 0 "building a distance index of a piped PackedGraph"

rm -f x.vg distIndex snarls.pb

# Test distance index with GBZ
vg gbwt -g graph.gbz --gbz-format -G graphs/components_walks.gfa
vg snarls -T graph.gbz > graph.snarls
vg index -j graph.dist graph.gbz
is $? 0 "distance index construction from GBZ"

cat graph.gbz | vg index -j graph.dist -
is $? 0 "distance index construction from piped GBZ"

rm -f graph.gbz graph.snarls graph.dist

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 17


# Indexing a single graph
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg

# Default construction
vg minimizer -i x.mi x.xg
is $? 0 "default parameters"

# Single-threaded
vg minimizer -t 1 -i x.mi x.xg
is $? 0 "single-threaded construction"
is $(md5sum x.mi | cut -f 1 -d\ ) b31ba27e61452bae522a6fc296965ed4 "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -i x.mi x.xg
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) 727a7990d534c13498833e9da63d6360 "setting -k -w works correctly"

# Max occs (-k 7 -w 3 -m 2)
vg minimizer -t 1 -k 7 -w 3 -m 2 -i x.mi x.xg
is $? 0 "max occurrences"
is $(md5sum x.mi | cut -f 1 -d\ ) f03f948f9ed0428d188f8c491a44e32f "frequent minimizers can be excluded"

# Haplotype-consistent minimizers
vg minimizer -t 1 -g x.gbwt -i x.mi x.xg
is $? 0 "haplotype-consistent minimizers"
is $(md5sum x.mi | cut -f 1 -d\ ) cc839bb44d8047e1ba60b1e0f4d193e8 "construction is deterministic"

rm -f x.vg x.xg x.gbwt x.mi


# Indexing two graphs
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg
vg index -x y.xg -G y.gbwt -v small/xy2.vcf.gz y.vg
vg index -x xy.xg -G xy.gbwt -v small/xy2.vcf.gz x.vg y.vg

# All minimizers
vg minimizer -t 1 -i x.mi x.xg
is $? 0 "all minimizers: first graph"
vg minimizer -t 1 -l x.mi -i xy.mi y.xg
is $? 0 "all minimizers: second graph"
vg minimizer -t 1 -i xy2.mi xy.xg
is $? 0 "all minimizers: both graphs"
cmp xy.mi xy2.mi
is $? 0 "the indexes are identical"

# Haplotype-consistent minimizers
vg minimizer -t 1 -g x.gbwt -i x.mi x.xg
is $? 0 "all minimizers: first graph"
vg minimizer -t 1 -g y.gbwt -l x.mi -i xy.mi y.xg
is $? 0 "all minimizers: second graph"
vg minimizer -t 1 -g xy.gbwt -i xy2.mi xy.xg
is $? 0 "all minimizers: both graphs"
cmp xy.mi xy2.mi
is $? 0 "the indexes are identical"

rm -f x.vg y.vg
rm -f x.xg y.xg xy.xg x.gbwt y.gbwt xy.gbwt
rm -f x.mi xy.mi xy2.mi

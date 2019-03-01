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
is $(md5sum x.mi | cut -f 1 -d\ ) 66abe4e3023de99363d54903b2b39f21 "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -i x.mi x.xg
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) 672db7e6be34f41c10001b1b9abcda2f "setting -k -w works correctly"

# Max occs (-k 7 -w 3 -m 2)
vg minimizer -t 1 -k 7 -w 3 -m 2 -i x.mi x.xg
is $? 0 "max occurrences"
is $(md5sum x.mi | cut -f 1 -d\ ) d1f001dd33828c5df4263f430f508e4c "frequent minimizers can be excluded"

# Haplotype-consistent minimizers
vg minimizer -t 1 -g x.gbwt -i x.mi x.xg
is $? 0 "haplotype-consistent minimizers"
is $(md5sum x.mi | cut -f 1 -d\ ) 7b66504fd7d0001e56613d23de5a937e "construction is deterministic"

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

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10


# Indexing a single graph
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg

# Default construction
vg minimizer -i x.mi x.xg
is $? 0 "default parameters"

# Single-threaded
vg minimizer -t 1 -i x.mi x.xg
is $? 0 "single-threaded construction"
is $(md5sum x.mi | cut -f 1 -d\ ) 5e32f24011133d89da496fbc9a2f7b29 "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -i x.mi x.xg
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) 1a5b8ab49a403adb3cadcd83833e1798 "setting -k -w works correctly"

# Haplotype-consistent minimizers
vg minimizer -t 1 -g x.gbwt -i x.mi x.xg
is $? 0 "haplotype-consistent minimizers"
is $(md5sum x.mi | cut -f 1 -d\ ) 6bedb0d725cfcab531a7938cdc8120a6 "construction is deterministic"

rm -f x.vg x.xg x.gbwt x.mi


# Indexing two graphs
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg
vg index -x y.xg -G y.gbwt -v small/xy2.vcf.gz y.vg

# Appending to the index
vg minimizer -t 1 -i x.mi x.xg
is $? 0 "multiple graphs: first"
vg minimizer -t 1 -l x.mi -i xy.mi y.xg
is $? 0 "multiple graphs: second"
is $(md5sum xy.mi | cut -f 1 -d\ ) 6597a86c79bdd63bc344c262267c3f7e "construction is deterministic"

rm -f x.vg y.vg
rm -f x.xg y.xg
rm -f x.mi xy.mi

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 12


# Indexing a single graph
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg

# Default construction
vg minimizer -i x.mi x.xg
is $? 0 "default parameters"

# Single-threaded
vg minimizer -t 1 -i x.mi x.xg
is $? 0 "single-threaded construction"
is $(md5sum x.mi | cut -f 1 -d\ ) 23216d7454a79dbc5e1cceb6f4613670 "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -i x.mi x.xg
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) 5a44c5a0b900bb4f5f32f76462f1eec9 "setting -k -w works correctly"

# Max occs (-k 7 -w 3 -m 2)
vg minimizer -t 1 -k 7 -w 3 -m 2 -i x.mi x.xg
is $? 0 "max occurrences"
is $(md5sum x.mi | cut -f 1 -d\ ) b3f091d939826322601c4ef570711da0 "frequent minimizers can be excluded"

# Haplotype-consistent minimizers
vg minimizer -t 1 -g x.gbwt -i x.mi x.xg
is $? 0 "haplotype-consistent minimizers"
is $(md5sum x.mi | cut -f 1 -d\ ) c786cf3fdad7356a55a2894d874c50da "construction is deterministic"

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
is $(md5sum xy.mi | cut -f 1 -d\ ) a7d670354123f9b990e1a492c51b3c29 "construction is deterministic"

rm -f x.vg y.vg
rm -f x.xg y.xg
rm -f x.mi xy.mi

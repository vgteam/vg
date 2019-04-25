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
is $(md5sum x.mi | cut -f 1 -d\ ) 859da428d4f90d10247cf3a0c8f6cafb "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -i x.mi x.xg
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) 08756cda9cd89b923ff4c8bf018f3257 "setting -k -w works correctly"

# Max occs (-k 7 -w 3 -m 2)
vg minimizer -t 1 -k 7 -w 3 -m 2 -i x.mi x.xg
is $? 0 "max occurrences"
is $(md5sum x.mi | cut -f 1 -d\ ) e917beec607f67269322107679c7d415 "frequent minimizers can be excluded"

# Haplotype-consistent minimizers
vg minimizer -t 1 -g x.gbwt -i x.mi x.xg
is $? 0 "haplotype-consistent minimizers"
is $(md5sum x.mi | cut -f 1 -d\ ) b53302959d910faa8b13b1a9e627e1d8 "construction is deterministic"

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
is $(md5sum xy.mi | cut -f 1 -d\ ) 45ccf9727fb3147b6133da12de2019f9 "construction is deterministic"

rm -f x.vg y.vg
rm -f x.xg y.xg
rm -f x.mi xy.mi

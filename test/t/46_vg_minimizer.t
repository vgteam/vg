#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10


# Indexing a single graph
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg

# Default construction
vg minimizer -i x.mi -g x.gbwt x.xg
is $? 0 "default parameters"

# Single-threaded for deterministic results
vg minimizer -t 1 -i x.mi -g x.gbwt x.xg
is $? 0 "single-threaded construction"
vg view --extract-tag MinimizerIndex x.mi > x.extracted.mi
is $(md5sum x.extracted.mi | cut -f 1 -d\ ) 6d173d906b6eead5074f1d7a26464e7c "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -i x.mi -g x.gbwt x.xg
is $? 0 "minimizer parameters"
vg view --extract-tag MinimizerIndex x.mi > x.extracted.mi
is $(md5sum x.extracted.mi | cut -f 1 -d\ ) fd3cea04319fe9ff7eb623788be4b31e "setting -k -w works correctly"
# FIXME

# Construction from GBWTGraph
vg gbwt -x x.xg -g x.gg x.gbwt
vg minimizer -t 1 -g x.gbwt -G -i x.mi x.gg
is $? 0 "construction from GBWTGraph"
vg view --extract-tag MinimizerIndex x.mi > x.extracted.mi
is $(md5sum x.extracted.mi | cut -f 1 -d\ ) 6d173d906b6eead5074f1d7a26464e7c "construction is deterministic"

rm -f x.vg x.xg x.gbwt x.mi x.extracted.mi x.gg


# Indexing two graphs
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg
vg index -x y.xg -G y.gbwt -v small/xy2.vcf.gz y.vg

# Appending to the index
vg minimizer -t 1 -i x.mi -g x.gbwt x.xg
is $? 0 "multiple graphs: first"
vg minimizer -t 1 -l x.mi -i xy.mi -g y.gbwt y.xg
is $? 0 "multiple graphs: second"
vg view --extract-tag MinimizerIndex xy.mi > xy.extracted.mi
is $(md5sum xy.extracted.mi | cut -f 1 -d\ ) f64285d9aeca03551c1a8a48e4c3f654 "construction is deterministic"

rm -f x.vg y.vg
rm -f x.xg y.xg
rm -f x.gbwt y.gbwt
rm -f x.mi xy.mi xy.extracted.mi

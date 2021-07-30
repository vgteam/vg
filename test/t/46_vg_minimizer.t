#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 16


# Indexing a single graph
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -x x.xg x.vg
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
vg snarls -T x.vg > x.snarls
vg index -s x.snarls -j x.dist x.xg

# Default construction
vg minimizer -o x.mi -g x.gbwt x.xg
is $? 0 "default parameters"

# Single-threaded for deterministic results
vg minimizer -t 1 -o x.mi -g x.gbwt x.xg
is $? 0 "single-threaded construction"
is $(md5sum x.mi | cut -f 1 -d\ ) 58b2bd98902df9acbe416bcfde649571 "construction is deterministic"

# Indexing syncmers
vg minimizer -t 1 -o x.mi -c -g x.gbwt x.xg
is $? 0 "syncmer index"
is $(md5sum x.mi | cut -f 1 -d\ ) 111bb0658db34d88a3a075f269513a36 "construction is deterministic"

# Minimizer parameters
vg minimizer -t 1 -k 7 -w 3 -o x.mi -g x.gbwt x.xg
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) c69f4f2dfab192fc66055dcf2fa8a7da "setting -k -w works correctly"

# Construction from GBWTGraph
vg gbwt -x x.xg -g x.gg x.gbwt
vg minimizer -t 1 -g x.gbwt -o x.mi x.gg
is $? 0 "construction from GBWTGraph"
is $(md5sum x.mi | cut -f 1 -d\ ) 58b2bd98902df9acbe416bcfde649571 "construction is deterministic"

# Construction from GBZ
vg gbwt -x x.xg -g x.gbz --gbz-format x.gbwt
vg minimizer -t 1 -o x.mi x.gbz
is $? 0 "construction from GBZ"
is $(md5sum x.mi | cut -f 1 -d\ ) 58b2bd98902df9acbe416bcfde649571 "construction is deterministic"

# Store payload in the index
vg minimizer -t 1 -o x.mi -g x.gbwt -d x.dist x.gg
is $? 0 "construction with payload"
is $(md5sum x.mi | cut -f 1 -d\ ) 58a6780c18921e4f6701b57fdb9c2e44 "construction is deterministic"

rm -f x.vg x.xg x.gbwt x.snarls x.dist x.mi x.gg x.gbz


# Indexing two graphs
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg
vg index -x y.xg -G y.gbwt -v small/xy2.vcf.gz y.vg

# Appending to the index
vg minimizer -t 1 -o x.mi -g x.gbwt x.xg
is $? 0 "multiple graphs: first"
vg minimizer -t 1 -l x.mi -o xy.mi -g y.gbwt y.xg
is $? 0 "multiple graphs: second"
is $(md5sum xy.mi | cut -f 1 -d\ ) 7500d106741f2dd8aa34aa178c917832 "construction is deterministic"

rm -f x.vg y.vg
rm -f x.xg y.xg
rm -f x.gbwt y.gbwt
rm -f x.mi xy.mi

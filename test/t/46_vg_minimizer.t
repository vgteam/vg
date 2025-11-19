#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 10


# Indexing a single graph
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg gbwt -x x.vg -v small/xy2.vcf.gz -g x.gbz --gbz-format
vg index -j x.dist x.gbz

# Default construction
vg minimizer --no-dist -o x.mi x.gbz
is $? 0 "default parameters"

# Construction fails without a distance index
vg minimizer -o x.mi x.gbz 2> /dev/null
is $? 1 "distance index or --no-dist is required"

# Single-threaded for deterministic results
vg minimizer --no-dist -t 1 -o x.mi x.gbz
is $? 0 "single-threaded construction"
is $(md5sum x.mi | cut -f 1 -d\ ) 8939d819aafa2183b393ee20ae6c4cc1 "construction is deterministic"

# Indexing syncmers
vg minimizer --no-dist -t 1 -o x.mi -c x.gbz
is $? 0 "syncmer index"
is $(md5sum x.mi | cut -f 1 -d\ ) ca5b291ba8096fb12ba0400a3dc833b2 "construction is deterministic"

# Minimizer parameters
vg minimizer --no-dist -t 1 -k 7 -w 3 -o x.mi x.gbz
is $? 0 "minimizer parameters"
is $(md5sum x.mi | cut -f 1 -d\ ) a7ee5b99f0bc0c0c7974681492548750 "setting -k -w works correctly"

# Store zipcode payload in the index
# Construction will not be deterministic because the snarls are not deterministic
vg minimizer -t 1 -o x.mi -z x.zipcodes -d x.dist x.gbz
is $? 0 "construction with zipcode payload"

# Store zipcode payload and path information in the index
# Construction will not be deterministic because the snarls are not deterministic
vg minimizer -t 1 -o x.mi -z x.zipcodes --rec-mode -d x.dist x.gbz
is $? 0 "construction with zipcode payload and path information"


rm -f x.vg x.dist x.mi x.zipcodes x.gbz

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 7

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg
vg construct -r small/x.fa -v small/x.vcf.gz > x2.vg
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg

# List path/thread names from various input formats
is "$(vg paths --list -v x2.vg)" "x" "path listing works from vg"
is "$(vg paths --list -x x.xg)" "x" "path listing works from XG"
is $(vg paths --list -g x.gbwt | wc -l) 2 "thread listing works from GBWT"

# Select threads by name
is $(vg paths --list -Q _thread_1_x_0 -g x.gbwt | wc -l) 1 "thread selection by name prefix works correctly"
is $(vg paths --list -S 1 -g x.gbwt | wc -l) 2 "thread selection by sample name works correctly"
is $(vg paths --list -S 2 -g x.gbwt | wc -l) 0 "no threads are reported for invalid samples"

# Extract threads as alignments
is $(vg paths -x x.xg -g x.gbwt -X | vg view -a -  | wc -l) 2 "vg paths may be used to extract threads"

rm -f x.xg x.gbwt x.vg x2.vg

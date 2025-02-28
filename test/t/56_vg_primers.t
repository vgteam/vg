#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 11

# make graph and snarl dist index
vg construct -r small/y.fa -v small/y.vcf.gz > y.vg

vg convert -x y.vg > y.xg

vg index -j y.dist y.vg

is $(vg primers primers/y.primer3_with_ref_pos.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz        | wc -l) 6 "Get the expected number of primer pairs"
is $(vg primers primers/y.primer3_with_ref_pos.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -a     | wc -l) 6 "Get the expected number of primer pairs using --all-primers tag"
is $(vg primers primers/y.primer3_with_ref_pos.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -l 2   | wc -l) 3 "Get the expected number of primer pairs using --tolerance tag"
is $(vg primers primers/y.primer3_with_ref_pos.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -n 137 | wc -l) 4 "Get the expected number of primer pairs using --minimum-size tag"
is $(vg primers primers/y.primer3_with_ref_pos.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -m 140 | wc -l) 4 "Get the expected number of primer pairs using --maximum-size tag"

is $(vg primers primers/y.split.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz        | wc -l) 9  "Get the expected number of primer pairs"
is $(vg primers primers/y.split.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -a     | wc -l) 11 "Get the expected number of primer pairs using --all-primers tag"
is $(vg primers primers/y.split.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -l 2   | wc -l) 6  "Get the expected number of primer pairs using --tolerance tag"
is $(vg primers primers/y.split.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -n 137 | wc -l) 4  "Get the expected number of primer pairs using --minimum-size tag"
is $(vg primers primers/y.split.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz -m 140 | wc -l) 7  "Get the expected number of primer pairs using --maximum-size tag"

vg primers primers/y.primer3_with_ref_pos.out    -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz> y.ref_pos_0.out
vg primers primers/y.primer3_with_ref_pos_11.out -x y.xg -d y.dist -r primers/y.ri -g primers/y.giraffe.gbz > y.ref_pos_11.out
diff -q <(awk '{$2=$3=$6=$7=""; print $0}' y.ref_pos_0.out) <(awk '{$2=$3=$6=$7=""; print $0}' y.ref_pos_11.out) > diff_0_11
is $(cat diff_0_11 | wc -l) 0 "These two output files should have identical primers except for their positions on template"

# clean up
rm diff_0_11
rm y.vg y.xg y.dist
rm y.ref_pos_0.out y.ref_pos_11.out

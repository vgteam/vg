#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 2

num_biallelic_variants=$(zcat < tiny/tiny.vcf.gz | grep -v "^#" | wc -l)
num_paths=$(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg paths -s -l 100 - | sort | uniq | wc -l)

is $num_paths $(echo "2^"$num_biallelic_variants | bc) "number of paths as expected for tiny graph"

is "$(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg paths --list -)" "x" "path listing works"

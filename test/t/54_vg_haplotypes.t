#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 4

# Build the indexes for the two-chromosome (1+1 kbp) test case
vg construct -r small/xy.fa -v small/xy.vcf.gz -a > small.vg 2> /dev/null
vg gbwt -v small/xy.vcf.gz -x small.vg --gbz-format -g small.gbz -r small.ri
vg index -j small.dist small.gbz

# Build the haplotype information and sample the haplotypes separately
vg haplotypes --validate --kmer-length 21 --subchain-length 300 -H small.hapl small.gbz
is $? 0 "generating haplotype information"
vg haplotypes --validate --num-haplotypes 1 --coverage 200 -i small.hapl -k small/xy-k21-coverage200.kff -g indirect.gbz small.gbz
is $? 0 "sampling from existing haplotype information"

# Sample the haplotypes directly
vg haplotypes --validate --kmer-length 21 --subchain-length 300 --num-haplotypes 1 --coverage 200 -k small/xy-k21-coverage200.kff -g direct.gbz small.gbz
is $? 0 "sampling the haplotypes directly"
cmp indirect.gbz direct.gbz
is $? 0 "the outputs are identical"

# FIXME: Test --include-reference with both named and reference paths.

# Cleanup
rm -r small.vg small.gbz small.ri small.dist
rm -f small.hapl indirect.gbz direct.gbz

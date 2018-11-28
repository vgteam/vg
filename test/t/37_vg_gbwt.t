#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 27


# Build vg graphs for two chromosomes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg


# Chromosome X
vg index -x x.xg -G x.gbwt -v small/xy2.vcf.gz x.vg
is $(vg gbwt -c x.gbwt) 2 "chromosome x: 2 threads"
is $(vg gbwt -C x.gbwt) 1 "chromosome x: 1 contig"
is $(vg gbwt -H x.gbwt) 2 "chromosome x: 2 haplotypes"
is $(vg gbwt -S x.gbwt) 1 "chromosome x: 1 sample"

# Full extraction of threads
is $(vg paths -x x.xg -g x.gbwt -X -T | vg view -a -  | wc -l) 2 "vg paths may be used to extract threads"

# Query test
is $(vg paths -x x.xg -g x.gbwt -X -Q _thread_1_x_0 | vg view -a -  | wc -l) 1 "vg paths can extract one thread by name prefix"

# Chromosome Y
vg index -G y.gbwt -v small/xy2.vcf.gz y.vg
is $(vg gbwt -c y.gbwt) 2 "chromosome y: 2 threads"
is $(vg gbwt -C y.gbwt) 1 "chromosome y: 1 contig"
is $(vg gbwt -H y.gbwt) 2 "chromosome y: 2 haplotypes"
is $(vg gbwt -S y.gbwt) 1 "chromosome y: 1 sample"

# Normal merging
vg gbwt -m -o xy.gbwt x.gbwt y.gbwt
is $? 0 "GBWT indexes can be merged"
is $(vg gbwt -c xy.gbwt) 4 "merge: 4 threads"
is $(vg gbwt -C xy.gbwt) 1 "merge: 1 contig"
is $(vg gbwt -H xy.gbwt) 4 "merge: 4 haplotypes"
is $(vg gbwt -S xy.gbwt) 2 "merge: 2 samples"

# Fast merging
vg gbwt -f -o xy2.gbwt x.gbwt y.gbwt
is $? 0 "GBWT indexes can be merged with the fast algorithm"
is $(vg gbwt -c xy2.gbwt) 4 "fast merge: 4 threads"
is $(vg gbwt -C xy2.gbwt) 2 "fast merge: 2 contigs"
is $(vg gbwt -H xy2.gbwt) 2 "fast merge: 2 haplotypes"
is $(vg gbwt -S xy2.gbwt) 1 "fast merge: 1 sample"

# Remove metadata from the merged indexes and compare them
../deps/gbwt/metadata -r xy > /dev/null
../deps/gbwt/metadata -r xy2 > /dev/null
cmp xy.gbwt xy2.gbwt
is $? 0 "the merged indexes are identical"

# Remove threads from a GBWT
vg gbwt -r 1 -r 2 xy.gbwt
is $? 0 "threads can be removed from a GBWT index"
is $(vg gbwt -c xy.gbwt) 2 "remove: 2 threads"

rm -f x.gbwt y.gbwt xy.gbwt xy2.gbwt x.xg


# Build a GBWT for paths
vg index -G x_ref.gbwt -T x.vg
is $(vg gbwt -c x_ref.gbwt) 1 "there is 1 thread in the index"

# Build a GBWT for both paths and threads
vg index -G x_both.gbwt -T -v small/xy2.vcf.gz x.vg
is $(vg gbwt -c x_both.gbwt) 3 "there are 3 threads in the index"

# Use vg paths to extract the gbwt


rm x_ref.gbwt x_both.gbwt


# Store haplotypes in both GBWT and a binary file
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg
vg index -H x.bin -v small/xy2.vcf.gz x.vg

# Extract threads from GBWT
vg gbwt -e x.extract x.gbwt
is $? 0 "threads can be extracted from GBWT"
cmp x.bin x.extract
is $? 0 "the thread files are identical"

rm -f x.gbwt x.bin x.extract


rm -f x.vg y.vg

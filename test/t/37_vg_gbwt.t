#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 56


# Build vg graphs for two chromosomes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg

vg index -x x.xg x.vg
vg index -x xy.xg x.vg y.vg
vg index -x xy-alt.xg -L x.vg y.vg


# Chromosome X
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg
is $(vg gbwt -c x.gbwt) 2 "chromosome x: 2 threads"
is $(vg gbwt -C x.gbwt) 1 "chromosome x: 1 contig"
is $(vg gbwt -H x.gbwt) 2 "chromosome x: 2 haplotypes"
is $(vg gbwt -S x.gbwt) 1 "chromosome x: 1 sample"

# Thread / contig / sample names
is $(vg gbwt -T x.gbwt | wc -l) 2 "chromosome x: 2 thread names"
is $(vg gbwt -C -L x.gbwt | wc -l) 1 "chromosome x: 1 contig name"
is $(vg gbwt -S -L x.gbwt | wc -l) 1 "chromosome x: 1 sample name"

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
is $(vg gbwt -C xy.gbwt) 2 "merge: 2 contigs"
is $(vg gbwt -H xy.gbwt) 2 "merge: 2 haplotypes"
is $(vg gbwt -S xy.gbwt) 1 "merge: 1 sample"

# Fast merging
vg gbwt -f -o xy2.gbwt x.gbwt y.gbwt
is $? 0 "GBWT indexes can be merged with the fast algorithm"
is $(vg gbwt -c xy2.gbwt) 4 "fast merge: 4 threads"
is $(vg gbwt -C xy2.gbwt) 2 "fast merge: 2 contigs"
is $(vg gbwt -H xy2.gbwt) 2 "fast merge: 2 haplotypes"
is $(vg gbwt -S xy2.gbwt) 1 "fast merge: 1 sample"

# Both merging algorithms should produce the same results, because metadata merging is based on names
cmp xy.gbwt xy2.gbwt
is $? 0 "the merged indexes are identical"

rm -f x.gbwt y.gbwt xy.gbwt xy2.gbwt


# Build a GBWT for paths
vg index -G x_ref.gbwt -T x.vg
is $(vg gbwt -c x_ref.gbwt) 1 "there is 1 thread in the index"

# Build a GBWT for both paths and threads
vg index -G x_haplo.gbwt -v small/xy2.vcf.gz x.vg
vg gbwt -m -o x_both.gbwt x_ref.gbwt x_haplo.gbwt
is $(vg gbwt -c x_both.gbwt) 3 "there are 3 threads in the index"

# Remove a sample (actually the reference) from a GBWT
vg gbwt -R ref x_both.gbwt -o removed.gbwt
is $? 0 "samples can be removed from a GBWT index"
is $(vg gbwt -c removed.gbwt) 2 "the sample was removed"

# Build a GBWT with paths as contigs
vg index -G xy_contigs.gbwt -T xy.xg
is $(vg gbwt -c xy_contigs.gbwt) 2 "paths as contigs: 2 threads"
is $(vg gbwt -C xy_contigs.gbwt) 2 "paths as contigs: 2 contigs"
is $(vg gbwt -H xy_contigs.gbwt) 1 "paths as contigs: 1 haplotype"
is $(vg gbwt -S xy_contigs.gbwt) 1 "paths as contigs: 1 sample"

# Build a GBWT with paths as samples
vg index -G xy_samples.gbwt -T --paths-as-samples xy.xg
is $(vg gbwt -c xy_samples.gbwt) 2 "paths as samples: 2 threads"
is $(vg gbwt -C xy_samples.gbwt) 1 "paths as samples: 1 contig"
is $(vg gbwt -H xy_samples.gbwt) 2 "paths as samples: 2 haplotypes"
is $(vg gbwt -S xy_samples.gbwt) 2 "paths as samples: 2 samples"

rm x_ref.gbwt x_haplo.gbwt x_both.gbwt removed.gbwt xy_contigs.gbwt xy_samples.gbwt


# Extract threads from GBWT
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg
vg gbwt -e x.extract x.gbwt
is $? 0 "threads can be extracted from GBWT"
is $(cat x.extract | wc -c) 121 "correct size for the thread file"

rm -f x.gbwt x.extract


# Build and serialize GBWTGraph
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg
vg gbwt -g x.gg -x x.xg x.gbwt
is $? 0 "GBWTGraph construction was successful"
vg view --extract-tag GBWTGraph x.gg > x.extracted.gg
is $(md5sum x.extracted.gg | cut -f 1 -d\ ) 62d451917c5076d7e84a6837dfb836cb "GBWTGraph was correctly serialized"

rm -f x.gbwt x.gg x.extracted.gg


# Build both GBWT and GBWTGraph from a 16-path cover
vg gbwt -P -n 16 -x xy.xg -g xy.gg -o xy.gbwt
is $? 0 "GBWT/GBWTGraph construction from path cover was successful"
vg view --extract-tag GBWTGraph xy.gg > xy.extracted.gg
is $(md5sum xy.extracted.gg | cut -f 1 -d\ ) ae6ba365e7e5fac6456f9a5a130aa98f "GBWTGraph was correctly serialized"
is $(vg gbwt -c xy.gbwt) 32 "path cover: 32 threads"
is $(vg gbwt -C xy.gbwt) 2 "path cover: 2 contigs"
is $(vg gbwt -H xy.gbwt) 16 "path cover: 16 haplotypes"
is $(vg gbwt -S xy.gbwt) 16 "path cover: 16 samples"

rm -f xy.gg xy.gbwt xy.extracted.gg


# Build both GBWT and GBWTGraph from 16 paths of local haplotypes
vg index -G xy.gbwt -v small/xy2.vcf.gz xy-alt.xg
vg gbwt -l -n 16 -x xy.xg -g xy.gg -o xy.local.gbwt xy.gbwt
is $? 0 "GBWT/GBWTGraph construction from local haplotypes was successful"
vg view --extract-tag GBWTGraph xy.gg > xy.extracted.gg
is $(md5sum xy.extracted.gg | cut -f 1 -d\ ) b7b40fb5296ded80cc659cd2300015af "GBWTGraph was correctly serialized"
is $(vg gbwt -c xy.local.gbwt) 32 "local haplotypes: 32 threads"
is $(vg gbwt -C xy.local.gbwt) 2 "local haplotypes: 2 contigs"
is $(vg gbwt -H xy.local.gbwt) 16 "local haplotypes: 16 haplotypes"
is $(vg gbwt -S xy.local.gbwt) 16 "local haplotypes: 16 samples"

rm -f xy.gg xy.gbwt xy.local.gbwt xy.extracted.gg


# Build GBWTGraph from an augmented GBWT
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg
vg gbwt -a -n 16 -x xy.xg -g augmented.gg -o augmented.gbwt x.gbwt
is $? 0 "GBWT/GBWTGraph construction by augmenting GBWT was successful"
vg view --extract-tag GBWTGraph augmented.gg > augmented.extracted.gg
is $(md5sum augmented.extracted.gg | cut -f 1 -d\ ) b7b40fb5296ded80cc659cd2300015af "GBWTGraph was correctly serialized"
is $(vg gbwt -c augmented.gbwt) 18 "augmented: 18 threads"
is $(vg gbwt -C augmented.gbwt) 2 "augmented: 2 contigs"
is $(vg gbwt -H augmented.gbwt) 2 "augmented: 2 haplotypes"
is $(vg gbwt -S augmented.gbwt) 17 "augmented: 17 samples"

rm -f x.gbwt augmented.gg augmented.gbwt augmented.extracted.gg


rm -f x.vg y.vg x.xg xy.xg xy-alt.xg

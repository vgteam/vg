#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 91


# Build vg graphs for two chromosomes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg

vg index -x x.xg x.vg
vg index -x xy.xg x.vg y.vg
vg index -x xy-alt.xg -L x.vg y.vg


# Single chromosome: haplotypes
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
is $? 0 "chromosome X GBWT with vg gbwt"
vg index -G x2.gbwt -v small/xy2.vcf.gz x.vg
is $? 0 "chromosome X GBWT with vg index"
cmp x.gbwt x2.gbwt
is $? 0 "identical construction results with vg gbwt and vg index"
vg gbwt -x x.vg -o parse --parse-only -v small/xy2.vcf.gz
is $? 0 "chromosome X VCF parse"
../deps/gbwt/build_gbwt -p -r parse_x > /dev/null 2> /dev/null
is $? 0 "chromosome X GBWT from VCF parse"
vg view -x GBWT x.gbwt > x.bare.gbwt
cmp x.bare.gbwt parse_x.gbwt
is $? 0 "identical construction results with vg gbwt and from VCF parse"

# Single chromosome: metadata for haplotypes
is $(vg gbwt -c x.gbwt) 2 "chromosome X: 2 threads"
is $(vg gbwt -C x.gbwt) 1 "chromosome X: 1 contig"
is $(vg gbwt -H x.gbwt) 2 "chromosome X: 2 haplotypes"
is $(vg gbwt -S x.gbwt) 1 "chromosome X: 1 sample"
is $(vg gbwt -T x.gbwt | wc -l) 2 "chromosome X: 2 thread names"
is $(vg gbwt -C -L x.gbwt | wc -l) 1 "chromosome X: 1 contig name"
is $(vg gbwt -S -L x.gbwt | wc -l) 1 "chromosome X: 1 sample name"

rm -f x.gbwt x2.gbwt x.bare.gbwt parse_x.gbwt
rm -f parse_x parse_x_0_1


# Single chromosome: paths
vg gbwt -E -o x.ref.gbwt -x x.vg
is $? 0 "chromosome X reference GBWT with vg gbwt"
vg index -G x2.ref.gbwt -T x.vg
is $? 0 "chromosome X reference GBWT with vg index"
cmp x.ref.gbwt x2.ref.gbwt
is $? 0 "identical construction results with vg gbwt and vg index"
is $(vg gbwt -c x.ref.gbwt) 1 "chromosome X reference: 1 thread"

rm -f x.ref.gbwt x2.ref.gbwt


# Single chromosome: alignments
vg paths -v x.vg -X -Q _alt > x.alts.gam
vg convert -G x.alts.gam x.vg > x.alts.gaf
vg gbwt -A -o x.alts.gaf.gbwt -x x.vg x.alts.gaf
is $? 0 "chromosome X GAF with vg gbwt"
vg index -F x.alts.gaf -G x2.alts.gaf.gbwt x.vg
is $? 0 "chromosome X GAF with vg index"
cmp x.alts.gaf.gbwt x2.alts.gaf.gbwt
is $? 0 "identical construction results with vg gbwt and vg index"
vg gbwt -A --gam-format -o x.alts.gam.gbwt -x x.vg x.alts.gam
is $? 0 "chromosome X GAM with vg gbwt"
cmp x.alts.gaf.gbwt x.alts.gaf.gbwt
is $? 0 "identical construction results from GAF and GAM"

rm -f x.alts.gam x.alts.gaf
rm -f x.alts.gaf.gbwt x2.alts.gaf.gbwt x.alts.gam.gbwt


# Graph region: haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a --region x:100-200 > x.part.vg
vg gbwt -x x.part.vg -o x.part.gbwt --vcf-region x:100-200 -v small/x.vcf.gz 2> log.txt
is $? 0 "chromosome X subgraph GBWT with vg gbwt"
is "$(cat log.txt | wc -c)" 0 "no warnings about missing variants"
vg index -G x2.part.gbwt --region x:100-200 -v small/x.vcf.gz x.part.vg 2> log.txt
is $? 0 "chromosome X subgraph GBWT with vg index"
cmp x.part.gbwt x2.part.gbwt
is $? 0 "identical construction results with vg gbwt and vg index"

rm -f x.part.vg x.part.gbwt x2.part.gbwt log.txt


# Multiple chromosomes: haplotypes
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
vg gbwt -x y.vg -o y.gbwt -v small/xy2.vcf.gz
vg gbwt -m -o xy.merge.gbwt x.gbwt y.gbwt
is $? 0 "insertion merging with multiple chromosomes"
vg gbwt -f -o xy.fast.gbwt x.gbwt y.gbwt
is $? 0 "fast merging with multiple chromosomes"
cmp xy.merge.gbwt xy.fast.gbwt
is $? 0 "identical merging results with the insertion and fast algorithms"
vg gbwt -b --merge-jobs 1 -o xy.parallel.gbwt x.gbwt y.gbwt
is $? 0 "parallel merging with multiple chromosomes"
cmp xy.merge.gbwt xy.parallel.gbwt
is $? 0 "identical merging results with the insertion and parallel algorithms"
vg gbwt -x xy-alt.xg -o xy.direct.gbwt -v small/xy2.vcf.gz
is $? 0 "direct construction with multiple chromosomes and a single VCF"
cmp xy.direct.gbwt xy.merge.gbwt
is $? 0 "identical results with direct construction and merging"
vg gbwt -x xy-alt.xg -o xy.multi.gbwt -v --inputs-as-jobs small/xy2_x.vcf.gz small/xy2_y.vcf.gz
is $? 0 "direct construction with multiple chromosomes and multiple VCFs"
cmp xy.direct.gbwt xy.multi.gbwt
is $? 0 "identical construction results with a single VCF and multiple VCFs"

# Multiple chromosomes: haplotypes with presets
vg gbwt -x xy-alt.xg -o xy.1000gp.gbwt --preset 1000gp -v small/xy2.vcf.gz
is $? 0 "construction preset: 1000gp"

# Multiple chromosomes: metadata for haplotypes
is $(vg gbwt -c xy.merge.gbwt) 4 "multiple chromosomes: 4 threads"
is $(vg gbwt -C xy.merge.gbwt) 2 "multiple chromosomes: 2 contigs"
is $(vg gbwt -H xy.merge.gbwt) 2 "multiple chromosomes: 2 haplotypes"
is $(vg gbwt -S xy.merge.gbwt) 1 "multiple chromosomes: 1 sample"

rm -f x.gbwt y.gbwt xy.merge.gbwt xy.fast.gbwt xy.parallel.gbwt xy.direct.gbwt xy.multi.gbwt
rm -f xy.1000gp.gbwt


# Multiple chromosomes: paths as contigs
vg gbwt -E -o xy.contigs.gbwt -x xy.xg
is $? 0 "paths as contigs with vg gbwt"
vg index -G xy2.contigs.gbwt -T xy.xg
is $? 0 "paths as contigs with vg index"
cmp xy.contigs.gbwt xy2.contigs.gbwt
is $? 0 "identical construction results with vg gbwt and vg index"
is $(vg gbwt -c xy.contigs.gbwt) 2 "paths as contigs: 2 threads"
is $(vg gbwt -C xy.contigs.gbwt) 2 "paths as contigs: 2 contigs"
is $(vg gbwt -H xy.contigs.gbwt) 1 "paths as contigs: 1 haplotype"
is $(vg gbwt -S xy.contigs.gbwt) 1 "paths as contigs: 1 sample"

# Multiple chromosomes: paths as samples
vg gbwt -E --paths-as-samples -o xy.samples.gbwt -x xy.xg
is $? 0 "paths as samples with vg gbwt"
vg index -G xy2.samples.gbwt -T --paths-as-samples xy.xg
is $? 0 "paths as samples with vg index"
cmp xy.samples.gbwt xy2.samples.gbwt
is $(vg gbwt -c xy.samples.gbwt) 2 "paths as samples: 2 threads"
is $(vg gbwt -C xy.samples.gbwt) 1 "paths as samples: 1 contig"
is $(vg gbwt -H xy.samples.gbwt) 2 "paths as samples: 2 haplotypes"
is $(vg gbwt -S xy.samples.gbwt) 2 "paths as samples: 2 samples"

rm -f xy.contigs.gbwt xy2.contigs.gbwt xy.samples.gbwt xy2.samples.gbwt


# Build an r-index
# TODO: More tests once we can use the r-index for something
vg gbwt -x xy-alt.xg -o xy.gbwt -r xy.ri --num-threads 1 -v small/xy2.vcf.gz
is $? 0 "r-index construction"

rm -f xy.gbwt xy.ri


# Non-empty and empty GBWTs
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
../deps/gbwt/build_gbwt -e empty > /dev/null

# Normal merging, x + empty
vg gbwt -m -o x2.gbwt x.gbwt empty.gbwt
is $? 0 "normal merging: non-empty + empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Normal merging, empty + x; silence the warning about the merging order
vg gbwt -m -o x2.gbwt empty.gbwt x.gbwt 2> /dev/null
is $? 0 "normal merging: empty + non-empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Fast merging, x + empty
vg gbwt -f -o x2.gbwt x.gbwt empty.gbwt
is $? 0 "fast merging: non-empty + empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Fast merging, empty + x
vg gbwt -f -o x2.gbwt empty.gbwt x.gbwt
is $? 0 "fast merging: empty + non-empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Parallel merging, x + empty
vg gbwt -b --merge-jobs 1 -o x2.gbwt x.gbwt empty.gbwt
is $? 0 "parallel merging: non-empty + empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

# Parallel merging, empty + x; silence the warning about the merging order
vg gbwt -b --merge-jobs 1 -o x2.gbwt empty.gbwt x.gbwt 2> /dev/null
is $? 0 "parallel merging: empty + non-empty"
cmp x.gbwt x2.gbwt
is $? 0 "the index remains unchanged"

rm -f x.gbwt empty.gbwt x2.gbwt


# Build a GBWT with both haplotypes and embedded paths
vg gbwt -x xy-alt.xg -o xy.gbwt -v small/xy2.vcf.gz
vg gbwt -E -o xy.ref.gbwt -x xy.xg
vg gbwt -m -o xy.both.gbwt xy.gbwt xy.ref.gbwt
is $(vg gbwt -c xy.both.gbwt) 6 "haplotypes and paths: 6 threads"

# Remove the reference
vg gbwt -R ref -o xy.removed.gbwt xy.both.gbwt
is $? 0 "samples can be removed from a GBWT index"
is $(vg gbwt -c xy.removed.gbwt) 4 "haplotypes only: 4 threads"

rm -f xy.gbwt xy.ref.gbwt xy.both.gbwt xy.removed.gbwt


# Extract threads from GBWT
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
vg gbwt -e x.extract x.gbwt
is $? 0 "threads can be extracted from GBWT"
is $(cat x.extract | wc -c) 121 "correct size for the thread file"

rm -f x.gbwt x.extract


# Build and serialize GBWTGraph
vg gbwt -x x.vg -g x.gg -o x.gbwt -v small/xy2.vcf.gz
is $? 0 "GBWTGraph construction was successful"
vg view --extract-tag GBWTGraph x.gg > x.extracted.gg
is $(md5sum x.extracted.gg | cut -f 1 -d\ ) 62d451917c5076d7e84a6837dfb836cb "GBWTGraph was correctly serialized"

rm -f x.gbwt x.gg x.extracted.gg


# Build both GBWT and GBWTGraph from a 16-path cover
vg gbwt -P -n 16 -x xy.xg -g xy.cover.gg -o xy.cover.gbwt
is $? 0 "Path cover GBWTGraph construction was successful"
vg view --extract-tag GBWTGraph xy.cover.gg > xy.extracted.gg
is $(md5sum xy.extracted.gg | cut -f 1 -d\ ) ae6ba365e7e5fac6456f9a5a130aa98f "GBWTGraph was correctly serialized"
is $(vg gbwt -c xy.cover.gbwt) 32 "path cover: 32 threads"
is $(vg gbwt -C xy.cover.gbwt) 2 "path cover: 2 contigs"
is $(vg gbwt -H xy.cover.gbwt) 16 "path cover: 16 haplotypes"
is $(vg gbwt -S xy.cover.gbwt) 16 "path cover: 16 samples"

rm -f xy.cover.gg xy.cover.gbwt xy.extracted.gg


# Build both GBWT and GBWTGraph from 16 paths of local haplotypes
vg gbwt -x xy-alt.xg -g xy.local.gg -l -n 16 -o xy.local.gbwt -v small/xy2.vcf.gz
is $? 0 "Local haplotypes GBWTGraph construction was successful"
vg view --extract-tag GBWTGraph xy.local.gg > xy.extracted.gg
is $(md5sum xy.extracted.gg | cut -f 1 -d\ ) b7b40fb5296ded80cc659cd2300015af "GBWTGraph was correctly serialized"
is $(vg gbwt -c xy.local.gbwt) 32 "local haplotypes: 32 threads"
is $(vg gbwt -C xy.local.gbwt) 2 "local haplotypes: 2 contigs"
is $(vg gbwt -H xy.local.gbwt) 16 "local haplotypes: 16 haplotypes"
is $(vg gbwt -S xy.local.gbwt) 16 "local haplotypes: 16 samples"

rm -f xy.local.gg xy.local.gbwt xy.extracted.gg


# Build GBWTGraph from an augmented GBWT
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
vg gbwt -a -n 16 -x xy.xg -g augmented.gg -o augmented.gbwt x.gbwt
is $? 0 "Augmented GBWTGraph construction was successful"
vg view --extract-tag GBWTGraph augmented.gg > augmented.extracted.gg
is $(md5sum augmented.extracted.gg | cut -f 1 -d\ ) b7b40fb5296ded80cc659cd2300015af "GBWTGraph was correctly serialized"
is $(vg gbwt -c augmented.gbwt) 18 "augmented: 18 threads"
is $(vg gbwt -C augmented.gbwt) 2 "augmented: 2 contigs"
is $(vg gbwt -H augmented.gbwt) 2 "augmented: 2 haplotypes"
is $(vg gbwt -S augmented.gbwt) 17 "augmented: 17 samples"

rm -f x.gbwt augmented.gg augmented.gbwt augmented.extracted.gg


rm -f x.vg y.vg x.xg xy.xg xy-alt.xg

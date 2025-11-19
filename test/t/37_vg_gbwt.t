#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 149


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
vg gbwt -x x.vg -o parse --parse-only -v small/xy2.vcf.gz
is $? 0 "chromosome X VCF parse"
../deps/gbwt/bin/build_gbwt -p -r parse_x > /dev/null 2> /dev/null
is $? 0 "chromosome X GBWT from VCF parse"
cmp x.gbwt parse_x.gbwt
is $? 0 "identical construction results with vg gbwt and from VCF parse"

# Single chromosome: metadata for haplotypes
is $(vg gbwt -c x.gbwt) 2 "chromosome X: 2 paths"
is $(vg gbwt -C x.gbwt) 1 "chromosome X: 1 contig"
is $(vg gbwt -H x.gbwt) 2 "chromosome X: 2 haplotypes"
is $(vg gbwt -S x.gbwt) 1 "chromosome X: 1 sample"
is $(vg gbwt -T x.gbwt | wc -l) 2 "chromosome X: 2 path names"
is $(vg gbwt -C -L x.gbwt | wc -l) 1 "chromosome X: 1 contig name"
is $(vg gbwt -S -L x.gbwt | wc -l) 1 "chromosome X: 1 sample name"

rm -f x.gbwt parse parse_x.gbwt
rm -f parse_x parse_x_0_1


# Single chromosome: paths
vg gbwt -E -o x.ref.gbwt -x x.vg
is $? 0 "chromosome X reference GBWT with vg gbwt"
is $(vg gbwt -c x.ref.gbwt) 1 "chromosome X reference: 1 path"

rm -f x.ref.gbwt


# Single chromosome: alignments
vg paths -x x.vg -X -Q _alt > x.alts.gam
vg convert -G x.alts.gam x.vg > x.alts.gaf
vg gbwt -A --num-jobs 1 -o x.alts.gaf.gbwt -x x.vg x.alts.gaf
is $? 0 "chromosome X GAF with vg gbwt"
vg gbwt -A --num-jobs 1 --gam-format -o x.alts.gam.gbwt -x x.vg x.alts.gam
is $? 0 "chromosome X GAM with vg gbwt"
cmp x.alts.gaf.gbwt x.alts.gaf.gbwt
is $? 0 "identical construction results from GAF and GAM"

rm -f x.alts.gam x.alts.gaf
rm -f x.alts.gaf.gbwt x.alts.gam.gbwt


# Graph region: haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a --region x:100-200 > x.part.vg
vg gbwt -x x.part.vg -o x.part.gbwt --vcf-region x:100-200 -v small/x.vcf.gz 2> log.txt
is $? 0 "chromosome X subgraph GBWT with vg gbwt"
is "$(cat log.txt | wc -c)" 0 "no warnings about missing variants"

rm -f x.part.vg x.part.gbwt log.txt


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
is $(vg gbwt -c xy.merge.gbwt) 4 "multiple chromosomes: 4 paths"
is $(vg gbwt -C xy.merge.gbwt) 2 "multiple chromosomes: 2 contigs"
is $(vg gbwt -H xy.merge.gbwt) 2 "multiple chromosomes: 2 haplotypes"
is $(vg gbwt -S xy.merge.gbwt) 1 "multiple chromosomes: 1 sample"

rm -f x.gbwt y.gbwt xy.merge.gbwt xy.fast.gbwt xy.parallel.gbwt xy.direct.gbwt xy.multi.gbwt
rm -f xy.1000gp.gbwt


# Multiple chromosomes: alignments
vg paths -x xy-alt.xg -X -Q _alt > xy.alts.gam
vg convert -G xy.alts.gam xy.xg > xy.alts.gaf
vg gbwt -A --num-jobs 1 -o xy.alts.gaf.gbwt -x xy.xg xy.alts.gaf
is $? 0 "multi-chromosome GAF with vg gbwt"
vg gbwt -A --num-jobs 1 --gam-format -o xy.alts.gam.gbwt -x xy.xg xy.alts.gam
is $? 0 "multi-chromosome GAM with vg gbwt"
cmp xy.alts.gaf.gbwt xy.alts.gaf.gbwt
is $? 0 "identical construction results from GAF and GAM"

vg gbwt -A --num-jobs 2 -o multi.gbwt -x xy.xg xy.alts.gaf
is $? 0 "multi-chromosome GAF with vg gbwt using multiple jobs"
is $(vg gbwt -c xy.alts.gaf.gbwt) 58 "single job: 58 paths"
is $(vg gbwt -c multi.gbwt) 58 "multiple jobs: 58 paths"

rm -f xy.alts.gam xy.alts.gaf
rm -f xy.alts.gaf.gbwt xy.alts.gam.gbwt multi.gbwt


# Multiple chromosomes: paths as contigs
vg gbwt -E -o xy.contigs.gbwt -x xy.xg
is $? 0 "paths as contigs with vg gbwt"
is $(vg gbwt -c xy.contigs.gbwt) 2 "paths as contigs: 2 paths"
is $(vg gbwt -C xy.contigs.gbwt) 2 "paths as contigs: 2 contigs"
is $(vg gbwt -H xy.contigs.gbwt) 1 "paths as contigs: 1 haplotype"
is $(vg gbwt -S xy.contigs.gbwt) 1 "paths as contigs: 1 sample"

rm -f xy.contigs.gbwt


# Build an r-index
# TODO: More tests once we can use the r-index for something
vg gbwt -x xy-alt.xg -o xy.gbwt -r xy.ri --num-threads 1 -v small/xy2.vcf.gz
is $? 0 "r-index construction"

rm -f xy.gbwt xy.ri


# Non-empty and empty GBWTs
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
../deps/gbwt/bin/build_gbwt -e empty > /dev/null

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
is $(vg gbwt -c xy.both.gbwt) 6 "haplotypes and paths: 6 paths"

# Remove the reference sample that GBWTs use for paths
vg gbwt -R _gbwt_ref -o xy.removed.gbwt xy.both.gbwt
is $? 0 "samples can be removed from a GBWT index"
is $(vg gbwt -c xy.removed.gbwt) 4 "haplotypes only: 4 paths"

rm -f xy.gbwt xy.ref.gbwt xy.both.gbwt xy.removed.gbwt


# Build a three-sample GBWT from a simple GFA
vg gbwt -o all.gbwt -G graphs/three_samples.gfa
is $(vg gbwt -c all.gbwt) 12 "all samples: 12 paths"
is $(vg gbwt -H all.gbwt) 6 "all samples: 6 haplotypes"

# Remove samples 1 and 3
vg gbwt -R sample1 -R sample3 -o removed.gbwt all.gbwt
is $? 0 "multiple samples can be removed from a GBWT index"
is $(vg gbwt -c removed.gbwt) 4 "sample 2: 4 paths"
is $(vg gbwt -H removed.gbwt) 2 "sample 2: 2 haplotypes"

rm -f all.gbwt removed.gbwt


# Extract paths from GBWT
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
vg gbwt -e x.extract x.gbwt
is $? 0 "paths can be extracted from GBWT"
is $(cat x.extract | wc -c) 121 "correct size for the paths file"

rm -f x.extract


# Build and serialize GBZ from an existing GBWT
vg gbwt -x x.vg -g x.gbz x.gbwt
is $? 0 "GBZ construction from GBWT"
is $(md5sum x.gbz | cut -f 1 -d\ ) c25b16aa5315e98be8f6cead13bc567c "GBZ was serialized correctly"

# Build and serialize GBZ from VCF
vg gbwt -x x.vg -g x2.gbz -v small/xy2.vcf.gz
is $? 0 "GBZ construction from VCF"
cmp x.gbz x2.gbz
is $? 0 "Identical construction results from GBWT and VCF"

# Extract GBWT from GBZ
vg gbwt -o extracted.gbwt -Z x.gbz
is $? 0 "GBWT extraction from GBZ"
cmp x.gbwt extracted.gbwt
is $? 0 "Identical GBWT indexes"

rm -f x.gbwt x.gbz x2.gbz
rm -f extracted.gbwt


# Build GBZ from a 16-path cover
vg gbwt -P -n 16 -x xy.xg -g xy.cover.gbz
is $? 0 "Path cover GBZ construction"
is $(md5sum xy.cover.gbz | cut -f 1 -d\ ) ffd93d2a6e603817bcba701db43a5ee1 "GBZ was serialized correctly"
is $(vg gbwt -c -Z xy.cover.gbz) 32 "path cover: 32 paths"
is $(vg gbwt -C -Z xy.cover.gbz) 2 "path cover: 2 contigs"
is $(vg gbwt -H -Z xy.cover.gbz) 16 "path cover: 16 haplotypes"
is $(vg gbwt -S -Z xy.cover.gbz) 16 "path cover: 16 samples"

rm -f xy.cover.gbz

# Build GBZ from a 16-path cover, passing through named paths
vg gbwt -P -n 16 -x xy.xg -g xy.cover.gbz --pass-paths
is $? 0 "Path cover GBZ construction"
is $(md5sum xy.cover.gbz | cut -f 1 -d\ ) f6e725bacf2fde5e8a69ae9c1384874c "GBZ was serialized correctly"
is $(vg gbwt -c -Z xy.cover.gbz) 34 "path cover w/ reference paths: 34 paths"
is $(vg gbwt -C -Z xy.cover.gbz) 2 "path cover w/ reference paths: 2 contigs"
is $(vg gbwt -H -Z xy.cover.gbz) 17 "path cover w/ reference paths: 17 haplotypes"
is $(vg gbwt -S -Z xy.cover.gbz) 17 "path cover w/ reference paths: 17 samples"

rm -f xy.cover.gbz

# Build GBZ from 16 paths of local haplotypes with a single job
vg gbwt -x xy-alt.xg -g xy.local.gbz -l -n 16 --num-jobs 1 -v small/xy2.vcf.gz
is $? 0 "Local haplotypes GBZ construction"
is $(md5sum xy.local.gbz | cut -f 1 -d\ ) e0b32538b533a69d26afb266e227da52 "GBZ was serialized correctly"
is $(vg gbwt -c -Z xy.local.gbz) 32 "local haplotypes: 32 paths"
is $(vg gbwt -C -Z xy.local.gbz) 2 "local haplotypes: 2 contigs"
is $(vg gbwt -H -Z xy.local.gbz) 16 "local haplotypes: 16 haplotypes"
is $(vg gbwt -S -Z xy.local.gbz) 16 "local haplotypes: 16 samples"

rm -f xy.local.gbz

# As above, but with two parallel jobs
vg gbwt -x xy-alt.xg -g xy.local.gbz -l -n 16 --num-jobs 2 -v small/xy2.vcf.gz
is $? 0 "Local haplotypes GBZ construction (two jobs)"
is $(md5sum xy.local.gbz | cut -f 1 -d\ ) e0b32538b533a69d26afb266e227da52 "GBZ was serialized correctly"

rm -f xy.local.gbz

# Build GBZ from 16 paths of local haplotypes, passing through named paths
vg gbwt -x xy-alt.xg -g xy.local.gbz -l -n 16 -v small/xy2.vcf.gz --pass-paths
is $? 0 "Local haplotypes GBZ construction"
is $(md5sum xy.local.gbz | cut -f 1 -d\ ) a98013fb492ee0a8422a73a15f7e65c3 "GBZ was serialized correctly"
is $(vg gbwt -c -Z xy.local.gbz) 34 "local haplotypes w/ reference paths: 34 paths"
is $(vg gbwt -C -Z xy.local.gbz) 2 "local haplotypes w/ reference paths: 2 contigs"
is $(vg gbwt -H -Z xy.local.gbz) 17 "local haplotypes w/ reference paths: 17 haplotypes"
is $(vg gbwt -S -Z xy.local.gbz) 17 "local haplotypes w/ reference paths: 17 samples"

rm -f xy.local.gbz

# Build GBZ from a GFA and then build a local haplotype cover with reference paths from the GBZ
vg gbwt -G haplotype-sampling/micb-kir3dl1.gfa -g large.gbz
vg gbwt -Z large.gbz -l -n 16 --pass-paths -o large.local.gbwt
is $? 0 "Local haplotypes with reference paths from a larger GBZ"
is $(vg gbwt -c large.local.gbwt) 36 "local haplotypes w/ reference paths: 36 paths"
is $(vg gbwt -C large.local.gbwt) 2 "local haplotypes w/ reference paths: 2 contigs"
is $(vg gbwt -H large.local.gbwt) 18 "local haplotypes w/ reference paths: 18 haplotypes"
is $(vg gbwt -S large.local.gbwt) 18 "local haplotypes w/ reference paths: 18 samples"
is $(vg gbwt --tags large.local.gbwt | grep -c reference_samples) 1 "local haplotypes w/ reference paths: reference_samples set"

# Now the same but with a specified reference sample
vg gbwt -Z large.gbz -l -n 16 --pass-paths --set-reference GRCh38 -o large.local.GCRh38.gbwt
is $? 0 "Local haplotypes with a specified reference sample"
is $(vg gbwt -c large.local.GCRh38.gbwt) 34 "local haplotypes w/ reference paths: 34 paths"
is $(vg gbwt -C large.local.GCRh38.gbwt) 2 "local haplotypes w/ reference paths: 2 contigs"
is $(vg gbwt -H large.local.GCRh38.gbwt) 17 "local haplotypes w/ reference paths: 17 haplotypes"
is $(vg gbwt -S large.local.GCRh38.gbwt) 17 "local haplotypes w/ reference paths: 17 samples"
is $(vg gbwt --tags large.local.GCRh38.gbwt | grep -c reference_samples) 1 "local haplotypes w/ reference paths: reference_samples set"

rm -f large.gbz large.local.gbwt large.local.GCRh38.gbwt


# Build GBZ from an augmented GBWT
vg gbwt -x x.vg -o x.gbwt -v small/xy2.vcf.gz
vg gbwt -a -n 16 -x xy.xg -g augmented.gbz x.gbwt
is $? 0 "Augmented GBZ construction"
is $(md5sum augmented.gbz | cut -f 1 -d\ ) 6333d0929717a14bf49941b3ec7ab87e "GBZ was serialized correctly"
is $(vg gbwt -c -Z augmented.gbz) 18 "augmented: 18 paths"
is $(vg gbwt -C -Z augmented.gbz) 2 "augmented: 2 contigs"
is $(vg gbwt -H -Z augmented.gbz) 18 "augmented: 18 haplotypes"
is $(vg gbwt -S -Z augmented.gbz) 17 "augmented: 17 samples"

rm -f x.gbwt augmented.gbz


# Remove the graphs
rm -f x.vg y.vg x.xg xy.xg xy-alt.xg


# Build GBWT from GFA
vg gbwt -o gfa.gbwt -G graphs/components_walks.gfa
is $? 0 "GBWT construction from GFA"
is $(md5sum gfa.gbwt | cut -f 1 -d\ ) 44c27c37c7af6911c26aea2a41008460 "GBWT was serialized correctly"
is $(vg gbwt -c gfa.gbwt) 4 "gfa: 4 paths"
is $(vg gbwt -C gfa.gbwt) 2 "gfa: 2 contigs"
is $(vg gbwt -H gfa.gbwt) 2 "gfa: 2 haplotypes"
is $(vg gbwt -S gfa.gbwt) 1 "gfa: 1 sample"

# Build GBZ from GFA
vg gbwt -g gfa2.gbz -G graphs/components_walks.gfa
is $? 0 "GBZ construction from GFA"
is $(md5sum gfa2.gbz | cut -f 1 -d\ ) 125428508faa78b96c4e752374a534b1 "GBZ was serialized correctly"
vg gbwt -Z gfa2.gbz -o /dev/null --translation gfa2.trans
is $(wc -l < gfa2.trans) 0 "no chopping: 0 translations"

# Build GBZ from GFA with node chopping
vg gbwt -g chopping.gbz --max-node 2 --translation chopping.trans -G graphs/chopping_walks.gfa
is $? 0 "GBZ construction from GFA with chopping"
is $(vg gbwt -c -Z chopping.gbz) 3 "chopping: 3 paths"
is $(vg gbwt -C -Z chopping.gbz) 1 "chopping: 1 contig"
is $(vg gbwt -H -Z chopping.gbz) 3 "chopping: 3 haplotypes"
is $(vg gbwt -S -Z chopping.gbz) 2 "chopping: 2 samples"
is $(wc -l < chopping.trans) 9 "chopping: 9 translations"
vg gbwt -Z chopping.gbz -o /dev/null --translation from_gbz.trans
is $? 0 "Translation can be extracted from GBZ"
is $(wc -l < from_gbz.trans) 8 "from GBZ: 8 translations"

# Build GBZ from GFA with both paths and walks
vg gbwt -g ref_paths.gbz --translation ref_paths.trans -G graphs/components_paths_walks.gfa
is $? 0 "GBZ construction from GFA with reference paths"
is $(vg gbwt -c -Z ref_paths.gbz) 6 "ref paths: 6 paths"
is $(vg gbwt -C -Z ref_paths.gbz) 2 "ref paths: 2 contigs"
is $(vg gbwt -H -Z ref_paths.gbz) 3 "ref paths: 3 haplotypes"
is $(vg gbwt -S -Z ref_paths.gbz) 2 "ref paths: 2 samples"
is $(wc -l < ref_paths.trans) 0 "ref paths: 0 translations"

rm -f gfa.gbwt
rm -f gfa2.trans gfa2.gbz
rm -f ref_paths.gbz ref_paths.trans
rm -f chopping.gbz chopping.trans from_gbz.trans

# Build a GBZ from a graph with a reference
vg gbwt -g gfa.gbz -G graphs/gfa_with_reference.gfa
is $? 0 "GBZ construction from GFA with reference"
vg gbwt -Z gfa.gbz --tags >tags.tsv
is $? 0 "GBZ GBWT tag extraction works"
is "$(grep reference_samples tags.tsv | cut -f2 | tr ' ' '\\n' | sort | tr '\\n' ' ')" "GRCh37 GRCh38" "GBWT tags contain the correct reference samples"
vg gbwt -g gfa2.gbz -Z gfa.gbz --set-tag "reference_samples=GRCh38 CHM13"
is $? 0 "GBZ GBWT tag modification works"
is "$(vg paths -M -S GRCh37 -x gfa2.gbz | grep -v "^#" | cut -f2 | grep HAPLOTYPE | wc -l)" "1" "Changing reference_samples tag can make a reference a haplotype"
is "$(vg paths -M -S CHM13 -x gfa2.gbz | grep -v "^#" | cut -f2 | grep REFERENCE | wc -l)" "1" "Changing reference_samples tag can make a haplotype a reference"
vg gbwt -g gfa2.gbz -Z gfa.gbz --set-tag "reference_samples=GRCh38#1 CHM13" 2>/dev/null
is $? 1 "GBZ GBWT tag modification validation works"
vg gbwt -g gfa3.gbz --set-reference GRCh37 --set-reference CHM13 -Z gfa2.gbz
is $? 0 "Samples can be direcly set as references"
is "$(vg gbwt --tags -Z gfa3.gbz | grep reference_samples | cut -f 2)" "CHM13 GRCh37" "Direct reference assignment works"

rm -f gfa.gbz gfa2.gbz gfa3.gbz tags.tsv

# Build a GBZ from a graph with a reference but no haplotype phase number
vg gbwt -g gfa.gbz -G graphs/gfa_two_part_reference.gfa
is "$(vg paths -M --reference-paths -x gfa.gbz | grep -v "^#" | cut -f4 | grep NO_HAPLOTYPE | wc -l)" "2" "GBZ can represent reference paths without haplotype numbers"

rm -f gfa.gbz



#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 27

# The test graph consists of two subgraphs of the HPRC Minigraph-Cactus v1.1 graph:
# - GRCh38#chr6:31498145-31511124 (micb)
# - GRCh38#chr19:54816468-54830778 (kir3dl1)
# The reads are 50x novaseq reads for HG003 mapping to those regions.

# Build indexes for the full graph.
vg gbwt -G haplotype-sampling/micb-kir3dl1.gfa --gbz-format -g full.gbz -r full.ri
vg index -j full.dist full.gbz

# Build the haplotype information and sample the haplotypes separately
vg haplotypes --validate --subchain-length 300 -H full.hapl full.gbz
is $? 0 "generating haplotype information"
vg haplotypes --validate -i full.hapl -k haplotype-sampling/HG003.kff --include-reference -g indirect.gbz full.gbz
is $? 0 "sampling from existing haplotype information"
is $(vg gbwt -S -Z indirect.gbz) 3 "1 generated + 2 reference samples"
is $(vg gbwt -C -Z indirect.gbz) 2 "2 contigs"
is $(vg gbwt -H -Z indirect.gbz) 6 "4 generated + 2 reference haplotypes"

# Sample the haplotypes directly
vg haplotypes --validate --subchain-length 300 -k haplotype-sampling/HG003.kff --include-reference -g direct.gbz full.gbz
is $? 0 "sampling the haplotypes directly"
cmp indirect.gbz direct.gbz
is $? 0 "the outputs are identical"

# Sample without reference
vg haplotypes --validate -i full.hapl -k haplotype-sampling/HG003.kff -g no_ref.gbz full.gbz
is $? 0 "sampling without adding reference"
is $(vg gbwt -S -Z no_ref.gbz) 1 "1 sample"
is $(vg gbwt -C -Z no_ref.gbz) 2 "2 contigs"
is $(vg gbwt -H -Z no_ref.gbz) 4 "4 haplotypes"

# Diploid sampling
vg haplotypes --validate -i full.hapl -k haplotype-sampling/HG003.kff --include-reference --diploid-sampling -g diploid.gbz full.gbz
is $? 0 "diploid sampling"
is $(vg gbwt -S -Z diploid.gbz) 3 "1 generated + 2 reference samples"
is $(vg gbwt -C -Z diploid.gbz) 2 "2 contigs"
is $(vg gbwt -H -Z diploid.gbz) 4 "2 generated + 2 reference haplotypes"

# Diploid sampling using a preset
vg haplotypes -i full.hapl -k haplotype-sampling/HG003.kff --preset diploid -g diploid2.gbz full.gbz
is $? 0 "diploid sampling using a preset"
cmp diploid.gbz diploid2.gbz
is $? 0 "the outputs are identical"

# Diploid sampling with a specified reference sample
vg haplotypes -i full.hapl -k haplotype-sampling/HG003.kff --include-reference --diploid-sampling --set-reference GRCh38 -g diploid3.gbz full.gbz
is $? 0 "diploid sampling with a specified reference sample"
is $(vg gbwt -S -Z diploid3.gbz) 2 "1 generated + 1 reference samples"
is $(vg gbwt -C -Z diploid3.gbz) 2 "2 contigs"
is $(vg gbwt -H -Z diploid3.gbz) 3 "2 generated + 1 reference haplotypes"

# Giraffe integration, guessed output name
vg giraffe -Z full.gbz --haplotype-name full.hapl --kff-name haplotype-sampling/HG003.kff \
    -f haplotype-sampling/HG003.fq.gz > default.gam 2> /dev/null
is $? 0 "Giraffe integration with a guessed output name"
cmp diploid.gbz full.HG003.gbz
is $? 0 "the sampled graph is identical to a manually sampled one"

# Giraffe integration, specified output name
vg giraffe -Z full.gbz --haplotype-name full.hapl --kff-name haplotype-sampling/HG003.kff \
    --index-basename sampled -N 003HG \
    -f haplotype-sampling/HG003.fq.gz > specified.gam 2> /dev/null
is $? 0 "Giraffe integration with a specified output name"
cmp full.HG003.gbz sampled.003HG.gbz
is $? 0 "the sampled graphs are identical"

# Giraffe integration, specified reference sample
vg giraffe -Z full.gbz --haplotype-name full.hapl --kff-name haplotype-sampling/HG003.kff \
    --index-basename GRCh38 -N HG003 --set-reference GRCh38 \
    -f haplotype-sampling/HG003.fq.gz > HG003_GRCh38.gam 2> /dev/null
is $? 0 "Giraffe integration with a specified reference sample"
cmp diploid3.gbz GRCh38.HG003.gbz
is $? 0 "the sampled graph is identical to a manually sampled one"

# Cleanup
rm -r full.gbz full.ri full.dist full.hapl
rm -f indirect.gbz direct.gbz no_ref.gbz
rm -f diploid.gbz diploid2.gbz diploid3.gbz
rm -f full.HG003.* default.gam
rm -f sampled.003HG.* specified.gam
rm -f GRCh38.HG003.* HG003_GRCh38.gam

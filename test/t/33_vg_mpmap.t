#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 20


# Exercise the GBWT
# Index a couple of nearly identical contigs
vg construct -m 1000 -a -r small/xy.fa -v small/xy2.vcf.gz >xy2.vg
vg index -x xy2.xg -g xy2.gcsa -v small/xy2.vcf.gz --gbwt-name xy2.gbwt -k 16 xy2.vg

# We turn off the background model calibration with -B and ignore it with -P 1

# This read is part ref and part alt which matches a haplotype on X, but is possible on Y as well.
is "$(vg mpmap --suppress-mismapping -B -P 1 -n dna -x xy2.xg -g xy2.gcsa -f reads/xy2.match.fq -F GAM | vg view -aj - | jq '.mapping_quality')" "3" "MAPQ is 50% without haplotype info"
is "$(vg mpmap --suppress-mismapping -B -P 1 -n dna -x xy2.xg -g xy2.gcsa --gbwt-name xy2.gbwt -f reads/xy2.match.fq -F GAM | vg view -aj - | jq '.mapping_quality')" "4" "haplotype match can disambiguate"

# For paired end, don't do any fragment length estimation

is "$(vg mpmap -B -P 1 -n dna -I 200 -D 200 -x xy2.xg -g xy2.gcsa -f reads/xy2.matchpaired.fq -i -F GAM -t 1 | vg view -aj - | head -n1 | jq '.mapping_quality')" "3" "MAPQ is 50% when paired without haplotype info"
vg mpmap -B -P 1 -n dna -I 200 -D 200 -x xy2.xg -g xy2.gcsa --gbwt-name xy2.gbwt -f reads/xy2.matchpaired.fq -i -F GAM -t 1 > gbwt.gam
is "$(vg view -aj gbwt.gam | head -n1 | jq '.mapping_quality')" "4" "haplotype match can disambiguate paired"
is "$(vg view -aj gbwt.gam | head -n1 | jq '.annotation.haplotype_score_used')" "true" "use of haplotype-aware mapping is recorded"

rm -f gbwt.gam

# Now we map a read with genotype 0,1,0,1 against haplotypes 1,1,1,1|0,1,0,0 and 1,1,0,1|0,0,1,0 on X and Y
# First 2 variants are inserts; second 2 are SNPs
# The right place is haplotype 2 on X with a SNP mismatch. But unless you do multiple tracebacks it gets placed on Y at MAPQ 0

# We need to use snarl cutting to actually consider the alignment as not just a single subpath
vg snarls xy2.vg > xy2.snarls

is "$(vg mpmap -B -P 1 -n dna -x xy2.xg -g xy2.gcsa -s xy2.snarls -f reads/xy2.discordant.fq -F GAM -t 1 | vg view -aj - | jq -r '.path.mapping[0].position.node_id')" "50" "Haplotype-oblivious mapping places read on the wrong contig"
is "$(vg mpmap -B -P 1 -n dna -x xy2.xg -g xy2.gcsa -s xy2.snarls -f reads/xy2.discordant.fq -F GAM -t 1| vg view -aj - | jq '.mapping_quality')" "3" "Haplotype-oblivious mapping places read with MAPQ of 50%"
is "$(vg mpmap -B -P 1 -n dna -x xy2.xg -g xy2.gcsa --gbwt-name xy2.gbwt -s xy2.snarls -f reads/xy2.discordant.fq -t 1 -F GAM | vg view -aj - | jq -r '.path.mapping[0].position.node_id')" "1" "Haplotype-aware mapping places read on the right contig"
is "$(vg mpmap -B -P 1 -n dna -x xy2.xg -g xy2.gcsa --gbwt-name xy2.gbwt -s xy2.snarls -f reads/xy2.discordant.fq -t 1 -F GAM | vg view -aj - | jq '.mapping_quality')" "6" "Haplotype-aware mapping places read with MAPQ > 50%"

rm -f xy2.vg xy2.xg xy2.gcsa xy2.gcsa.lcp xy2.gbwt xy2.snarls

# Do a larger-scale test

vg index -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -k 16 graphs/refonly-lrc_kir.vg

vg mpmap -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -n dna -B -i -I 10 -D 50 -F GAM -t 1 | vg view -aj -  > temp_paired_alignment.json
vg mpmap -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -n dna -B -i -I 100000 -D 5 -F GAM -t 1 | vg view -aj -  > temp_distant_alignment.json
vg mpmap -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -n dna -B -i -F GAM -t 1 | vg view -aj - > temp_independent_alignment.json
paired_score=$(jq -r ".score" < temp_paired_alignment.json | awk '{ sum+=$1} END {print sum}')
independent_score=$(jq -r ".score" < temp_independent_alignment.json | awk '{ sum+=$1} END {print sum}')
is $(printf "%s\t%s\n" $paired_score $independent_score | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent have lower score than unrestricted alignments"

paired_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_paired_alignment.json | sort | rs -T | awk '{print ($2 - $1)}')
distant_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_distant_alignment.json | sort | rs -T | awk '{print ($2 - $1)}')
independent_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_independent_alignment.json | sort | rs -T | awk '{print ($2 - $1)}')
is $(printf "%s\t%s\n" $paired_range $independent_range | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent are closer together in node id space than unrestricted alignments"
is $(printf "%s\t%s\n" $paired_range $distant_range | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be near each other are closer together in node id space than those forced to be far apart"

rm -f temp_paired_alignment.json temp_distant_alignment.json temp_independent_alignment.json

vg sim -x graphs/refonly-lrc_kir.vg.xg -n 1000 -p 500 -l 100 -a > input.gam

vg mpmap -B -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -G input.gam -I 500 -D 100 -n dna -i -F GAM --no-qual-adjust > output.gam
is "$(vg view -aj output.gam | jq -c 'select(.fragment_next == null and .fragment_prev == null)' | wc -l)" "0" "small batches are still all paired in the output"

rm -f graphs/refonly-lrc_kir.vg.xg graphs/refonly-lrc_kir.vg.gcsa graphs/refonly-lrc_kir.vg.gcsa.lcp input.gam output.gam

# Test the anchor trimming

vg construct -m 1000 -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
vg index -x t.xg -g t.gcsa -k 16 t.vg

echo "@read1
CAAATAAGG
+
HHHHHHHHH
@read2
AAAATTTTCT
+
HHHHHHHHHH
@read3
CAAATAAGGT
+
HHHHHHHHHH" > t.fq

is "$(vg mpmap -B -n dna -x t.xg -g t.gcsa -f t.fq | vg view -Kj - | wc -l)" "3" "multipath mapping works in scenarios that trigger branch point trimming"

rm t.vg t.xg t.gcsa t.gcsa.lcp t.fq

# test spliced alignment

vg construct -m 32 -r small/x.fa -v small/x.vcf.gz > s.vg
vg snarls -T s.vg > s.snarls
vg index s.vg -x s.xg -g s.gcsa
vg index s.vg -j s.dist -s s.snarls

echo "@read
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAGATGA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" > s.fq 

is "$(vg mpmap -x s.xg -d s.dist -g s.gcsa -B -n rna -f s.fq | vg view -KG - | vg view -aj - | jq .score)" "105" "spliced alignments can be found when aligning RNA"

echo "@read1
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAG
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read2
TCATCTCCCCGCACCTTTGTCCTTAGTCCACAGGAAACTCTGCTGTCAGTAGTATCATCTCCATATTAGAGATA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" > s.fq 


is "$(vg mpmap -x s.xg -d s.dist -g s.gcsa -B -n rna -f s.fq -i -I 350 -D 10 | vg view -Kj - | grep connection | wc -l)" "1" "paired mapping can identify a splice junction on read 1"

echo "@read1
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read2
TCATCTCCCCGCACCTTTGTCCTTAGTCCACAGGAAACTTAAAATGGCCTGGGTAGCTTTGGATGTGGAGTAGTTAAAGGCCCTTGAGG
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" > s.fq 

is "$(vg mpmap -x s.xg -d s.dist -g s.gcsa -B -n rna -f s.fq -i -I 530 -D 10 | vg view -Kj - | grep connection | wc -l)" "1" "paired mapping can identify a splice junction on read 2"

is "$(vg mpmap -x s.xg -d s.dist -g s.gcsa -B -n rna -f s.fq -i -I 530 -D 10 -F SAM --report-group-mapq | grep GM:i: | wc -l)" "2" "HTS output contains group mapq annotation"

rm s.vg s.xg s.gcsa s.gcsa.lcp s.dist s.snarls s.fq

# Now make sure we randomly choose between equivalent mappings
vg construct -r small/x.fa > x.vg
vg sim -a -p 200 -v 10 -l 50 -n 1000 -s 12345 -x x.vg >x.gam

vg construct -r small/xy.fa >xy.vg
vg snarls -T xy.vg > xy.snarls
vg index xy.vg -x xy.xg -g xy.gcsa
vg index xy.vg -j xy.dist -s xy.snarls

vg mpmap -x xy.xg -d xy.dist -g xy.gcsa -G x.gam -F SAM -i --frag-mean 50 --frag-stddev 10 >xy.sam
X_HITS="$(cat xy.sam | grep -v "^@" | cut -f3 | grep x | wc -l)"
if [ "${X_HITS}" -lt 1200 ] && [ "${X_HITS}" -gt 800 ] ; then
    IN_RANGE="1"
else
    IN_RANGE="0"
fi
is "${IN_RANGE}" "1" "paired reads are evenly split between equivalent mappings"

vg mpmap -x xy.xg -d xy.dist -g xy.gcsa -G x.gam -F SAM >xy.sam
X_HITS="$(cat xy.sam | grep -v "^@" | cut -f3 | grep x | wc -l)"
if [ "${X_HITS}" -lt 1200 ] && [ "${X_HITS}" -gt 800 ] ; then
    IN_RANGE="1"
else
    IN_RANGE="0"
fi
is "${IN_RANGE}" "1" "unpaired reads are evenly split between equivalent mappings"

rm x.vg x.gam xy.vg xy.xg xy.gcsa xy.snarls xy.dist xy.sam



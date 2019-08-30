#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 4

# Toy example of hand-made pileup (and hand inspected truth) to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg view -J -v call/tiny.json > tiny.vg

# With an empty pileup and loci mode we should assert the primery path.
true > empty.gam
vg augment tiny.vg empty.gam -A empty_aug.gam > tiny_aug.vg
vg index tiny_aug.vg -x tiny_aug.xg
vg pack -x tiny_aug.xg -g empty_aug.gam -o tiny_aug.pack
vg call tiny_aug.xg -k tiny_aug.pack > tiny_aug.vcf

is $(grep -v '#' tiny_aug.vcf | wc -l) 0 "calling empty gam gives empty VCF"

rm -f tiny.vg tiny_aug.vg tiny_aug.xg empty_aug.gam tiny_aug.pack tiny_aug.vcf

echo '{"node": [{"id": 1, "sequence": "CGTAGCGTGGTCGCATAAGTACAGTAGATCCTCCCCGCGCATCCTATTTATTAAGTTAAT"}]}' | vg view -Jv - > test.vg
vg index -x test.xg -g test.gcsa -k 16 test.vg
true >reads.txt
for REP in seq 1 5; do
    echo 'CGTAGCGTGGTCGCATAAGTACAGTANATCCTCCCCGCGCATCCTATTTATTAAGTTAAT' >>reads.txt
done
vg map -x test.xg -g test.gcsa --reads reads.txt > test.gam
vg augment test.vg test.gam -A test_aug.gam > test_aug.vg
vg index test_aug.vg -x test_aug.xg
vg pack -x test_aug.xg -g test_aug.gam -o test_aug.pack
vg call test_aug.xg -k test_aug.pack > /dev/null

N_COUNT=$(vg view -j test.aug.vg | grep "N" | wc -l)

is "${N_COUNT}" "0" "N bases are not augmented into the graph"

rm -rf reads.txt test.vg test.gam test_aug.gam test_aug.vg test_aug.xg test_aug.pack

vg construct -r inverting/miniFasta.fa -v inverting/miniFasta_VCFinversion.vcf.gz -S > miniFastaGraph.vg
vg index -x miniFastaGraph.xg -g miniFastaGraph.gcsa miniFastaGraph.vg
vg sim -x miniFastaGraph.xg -n 1000 -l 30 -a > miniFasta.gam
vg map -G miniFasta.gam -g miniFastaGraph.gcsa -x miniFastaGraph.xg > miniFastaGraph.gam
vg augment  miniFastaGraph.vg miniFastaGraph.gam -A mappedminitest_aug.gam > mappedminitest_aug.vg
vg index mappedminitest_aug.vg -x mappedminitest_aug.xg
vg pack -x mappedminitest_aug.xg -g mappedminitest_aug.gam -o mappedminitest_aug.pack
vg call  mappedminitest_aug.xg -k mappedminitest_aug.pack > calledminitest.vcf

L_COUNT=$(cat calledminitest.vcf | grep "#" -v | wc -l)
is "${L_COUNT}" "1" "Called microinversion"
 
rm -f miniFastaGraph.vg miniFasta.gam miniFastaGraph.gam mappedminitest.aug.vg calledminitest.vcf mappedminitest.trans mappedminitest.support mappedminitest.pileup miniFastaGraph.xg miniFastaGraph.gcsa

## SV Genotyping test
# augment the graph with the alt paths
vg augment -i call/HGSVC_chr22_17119590_17880307.vg call/HGSVC_chr22_17119590_17880307_alts.gam > HGSVC_alts.vg
# compute the xg-index with alts!
vg index HGSVC_alts.vg -x HGSVC_alts.xg -L
# get the support
vg pack -x HGSVC_alts.xg -g call/HGSVC_chr22_17119590_17880307.gam -o HGSVC_alts.pack
# genotype the VCF
vg call HGSVC_alts.xg -k HGSVC_alts.pack -v call/HGSVC_chr22_17200000_17800000.vcf.gz -s HG00514 > HGSVC.vcf
# extract the "true" calls
gzip -dc call/HGSVC_chr22_17200000_17800000.vcf.gz | grep -v '#' | awk '{print $10}' | awk -F ':' '{print $1}' > baseline_gts.txt
# extract the called genotypes
grep -v '#' HGSVC.vcf | sort -k1,1d -k2,2n | awk '{print $10}' | awk -F ':' '{print $1}' | sed 's/\//\|/g' > gts.txt
DIFF_COUNT=$(diff -U 0 baseline_gts.txt gts.txt | grep ^@ | wc -l)
LESS_SIX=$(if (( $DIFF_COUNT < 8 )); then echo 1; else echo 0; fi)
is "${LESS_SIX}" "1" "Fewer than 6 differences between called and true SV genotypes" 

rm -f HGSVC_alts.vg HGSVC_alts.xg HGSVC_alts.pack HGSVC.vcf baseline_gts.txt gts.txt




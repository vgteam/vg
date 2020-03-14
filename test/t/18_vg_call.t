#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 6

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

rm -f tiny.vg tiny_aug.vg tiny_aug.xg empty_aug.gam tiny_aug.pack tiny_aug.vcf empty.gam

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
 
rm -f miniFastaGraph.vg miniFasta.gam miniFastaGraph.gam calledminitest.vcf mappedminitest.trans mappedminitest.support mappedminitest.pileup miniFastaGraph.xg miniFastaGraph.gcsa mappedminitest_aug.vg mappedminitest_aug.gam mappedminitest_aug.xg mappedminitest_aug.pack miniFastaGraph.gcsa.lcp

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

vg construct -a -r small/x.fa -v small/x.vcf.gz > x.vg
vg index -x x.xg x.vg -L
vg sim -s 1 -n 1000 -l 150 -x x.xg -a > sim.gam
vg pack -x x.xg -g sim.gam -o x.xg.cx
vg pack -x x.vg -g sim.gam -o x.vg.cx
vg snarls x.xg > x.snarls
vg call x.xg -k x.xg.cx -r x.snarls -t 1 > x.xg.vcf
vg call x.vg -k x.vg.cx -r x.snarls -t 1 > x.vg.vcf
diff x.xg.vcf x.vg.vcf
is "$?" 0 "call output same on vg as xg"

vg call x.xg -k x.xg.cx -r x.snarls -t 1 -v tiny/tiny.vcf.gz > x.xg.gt.vcf
vg call x.vg -k x.vg.cx -r x.snarls -t 1 -v tiny/tiny.vcf.gz > x.vg.gt.vcf

diff x.xg.gt.vcf x.vg.gt.vcf
is "$?" 0 "call output same on vg as xg"

rm -f x.vg x.xg sim.gam x.xg.cx x.vg.cx x.xg.vcf x.vg.vcf x.xg.gt.vcf x.vg.gt.vcf x.snarls

vg msga -f msgas/cycle.fa -b s1 -w 64 -t 1 >c.vg
vg index -x c.xg -g c.gcsa c.vg
# True alignment has 3 variants:
# TCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTC
# TCCCTCCTCAAAGGCTTCTCACTACTCCA-ATCAAAGCTACCCAGGCCATTTTAAGTTTC
#            *       *
cat msgas/cycle.fa | sed s/TCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTC/TCCCTCCTCAAAGGCTTCTCACTACTCCAATCAAAGCTACCCAGGCCATTTTAAGTTTC/ >m.fa
vg construct -r m.fa >m.vg
vg index -x m.xg m.vg
vg sim -n 200 -s 23823 -l 50 -x m.xg -a >m.sim
vg map -x c.xg -g c.gcsa -G m.sim >m.gam
vg augment c.vg m.gam -A m.aug.gam >c.aug.vg
vg index -x c.aug.xg c.aug.vg
vg pack -x c.aug.xg -g m.aug.gam -o m.aug.pack
vg call c.aug.xg -k m.aug.pack >m.vcf
is $(cat m.vcf | grep -v "^#" | grep -v "0/0" | wc -l) 3 "vg call finds true homozygous variants in a cyclic graph"
rm -f c.vg c.xg c.gcsa c.gcsa.lcp m.fa m.vg m.xg m.sim m.gam m.aug.gam c.aug.vg c.aug.xg m.aug.pack m.vcf


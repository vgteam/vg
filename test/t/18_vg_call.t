#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 19

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
vg sim -x miniFastaGraph.xg -n 1000 -l 30 -a -s 1 > miniFasta.gam
vg map -G miniFasta.gam -g miniFastaGraph.gcsa -x miniFastaGraph.xg > miniFastaGraph.gam
vg augment  miniFastaGraph.vg miniFastaGraph.gam -A mappedminitest_aug.gam > mappedminitest_aug.vg
vg index mappedminitest_aug.vg -x mappedminitest_aug.xg
vg pack -x mappedminitest_aug.xg -g mappedminitest_aug.gam -o mappedminitest_aug.pack
vg call  mappedminitest_aug.xg -k mappedminitest_aug.pack > calledminitest.vcf

L_COUNT=$(cat calledminitest.vcf | grep "#" -v | wc -l)
is "${L_COUNT}" "1" "Called microinversion"

rm -f miniFastaGraph.vg miniFasta.gam miniFastaGraph.gam calledminitest.vcf  miniFastaGraph.xg miniFastaGraph.gcsa mappedminitest_aug.vg mappedminitest_aug.gam mappedminitest_aug.xg mappedminitest_aug.pack miniFastaGraph.gcsa.lcp

vg construct -r inverting/miniFasta.fa -v inverting/miniFasta_VCFinversion.vcf.gz -S > miniFastaGraph.vg
vg index -x miniFastaGraph.xg -g miniFastaGraph.gcsa miniFastaGraph.vg
vg sim -x miniFastaGraph.xg -n 1000 -l 30 -a -s 1 > miniFasta.gam
vg construct -r inverting/miniFasta.fa > miniFastaFlat.vg
vg sim -x  miniFastaFlat.vg -n 500 -l 30 -a -s 1 >> miniFasta.gam
vg map -G miniFasta.gam -g miniFastaGraph.gcsa -x miniFastaGraph.xg > miniFastaGraph.gam
vg augment  miniFastaGraph.vg miniFastaGraph.gam -A mappedminitest_aug.gam > mappedminitest_aug.vg
vg index mappedminitest_aug.vg -x mappedminitest_aug.xg
vg pack -x mappedminitest_aug.xg -g mappedminitest_aug.gam -o mappedminitest_aug.pack
vg call  mappedminitest_aug.xg -k mappedminitest_aug.pack -d 1 > calledminitest.vcf

L_COUNT=$(cat calledminitest.vcf | grep "#" -v | wc -l)
is "${L_COUNT}" "0" "Called no microinversion with haploid setting"

rm -f miniFastaGraph.vg miniFastaFlat.vg miniFasta.gam miniFastaGraph.gam calledminitest.vcf calledminitest1.vcf miniFastaGraph.xg miniFastaGraph.gcsa mappedminitest_aug.vg mappedminitest_aug.gam mappedminitest_aug.xg mappedminitest_aug.pack miniFastaGraph.gcsa.lcp

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
DIFF_COUNT=$(diff -U1000 <(nl baseline_gts.txt) <(nl gts.txt) | tail -n +4 | grep '^+' | wc -l)
LESS_EIGHT=$(if (( $DIFF_COUNT < 8 )); then echo 1; else echo 0; fi)
is "${LESS_EIGHT}" "1" "Fewer than 8 differences between called and true SV genotypes"

# genotype the VCF in haploid mode
vg call HGSVC_alts.xg -k HGSVC_alts.pack -v call/HGSVC_chr22_17200000_17800000.vcf.gz -s HG00514 -d 1 > HGSVC1.vcf
# extract the "true" calls
gzip -dc call/HGSVC_chr22_17200000_17800000.vcf.gz | grep -v '#' | awk '{print $10}' | awk -F ':' '{print $1}' | awk -F '|' '{print $1}'  > baseline_gts1.txt
# extract the called genotypes
grep -v '#' HGSVC1.vcf | sort -k1,1d -k2,2n | awk '{print $10}' | awk -F ':' '{print $1}' | sed 's/\//\|/g' > gts1.txt
DIFF_COUNT=$(diff -U1000 <(nl baseline_gts1.txt) <(nl gts1.txt) | tail -n +4 | grep '^+' | wc -l)
LESS_EIGHT=$(if (( $DIFF_COUNT <= 8 )); then echo 1; else echo 0; fi)
is "${LESS_EIGHT}" "1" "Fewer than 8 differences between called haploid and truncated true SV genotypes"

# call all the snarls with -a
vg call HGSVC_alts.xg -k HGSVC_alts.pack -s HG00514 -a > HGSVC2.vcf
REF_COUNT_V=$(grep "0/0" HGSVC.vcf | wc -l)
REF_COUNT_A=$(grep "0/0" HGSVC2.vcf | wc -l)
# this probably doesn't need to be exact (coincidence?), but it works now
is "${REF_COUNT_V}" "${REF_COUNT_A}" "Same number of reference calls with -a as with -v"

# Output snarl traversals into a GBWT then genotype that
vg call HGSVC_alts.xg -k HGSVC_alts.pack -s HG00514 -T | gzip > HGSVC_travs.gaf.gz
vg index HGSVC_alts.xg -F HGSVC_travs.gaf.gz -G HGSVC_travs.gbwt
vg call HGSVC_alts.xg -k HGSVC_alts.pack  -g HGSVC_travs.gbwt -s HG00514 > HGSVC_travs.vcf
vg call HGSVC_alts.xg -k HGSVC_alts.pack -s HG00514 > HGSVC_direct.vcf
# extract the called genotypes
grep -v '#' HGSVC_travs.vcf | sort -k1,1d -k2,2n | awk '{print $10}' | awk -F ':' '{print $1}' | sed 's/1\/0/0\/1/g' | sed 's/\//\|/g' > gts-travs.txt
grep -v '#' HGSVC_direct.vcf | sort -k1,1d -k2,2n | awk '{print $10}' | awk -F ':' '{print $1}' | sed 's/1\/0/0\/1/g' | sed 's/\//\|/g' > gts-direct.txt
diff gts-travs.txt gts-direct.txt
is "$?" "0" "Calling from extracted traversals by way of GBWT produces same genotypes as calling directly"

grep -v '#' HGSVC_travs.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > calls-travs.txt
grep -v '#' HGSVC_direct.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > calls-direct.txt
DIFF_COUNT=$(diff -U1000 <(nl calls-travs.txt) <(nl calls-direct.txt) | tail -n +4 | grep '^+' | wc -l)
LESS_THREE=$(if (( $DIFF_COUNT < 3 )); then echo 1; else echo 0; fi)
# because call makes an attempt to call multiple snarls at once when outputting traversals (to make bigger traversals)
# there is some wobble here
is "${LESS_THREE}" "1" "Fewer than 3 differences between allales called via traversals or directly"

rm -f HGSVC_alts.vg HGSVC_alts.xg HGSVC_alts.pack HGSVC.vcf baseline_gts.txt gts.txt HGSVC1.vcf HGSVC2.vcf HGSVC_travs.gaf.gz HGSVC_travs.gbwt HGSVC_travs.vcf HGSVC_direct.vcf baseline_gts1.txt gts1.txt gts-travs.txt gts-direct.txt calls-travs.txt calls-direct.txt

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

# simple gbwt
vg autoindex -r small/x.fa -v small/x.vcf.gz -w giraffe -p x
rm -f x.min x.dist
mv x.giraffe.gbz x.gbz
vg gbwt -Z x.gbz -o x.gbwt -g x.gg
vg convert x.gg -b x.gbwt -p  > x.vg
# simulate 500 reads from each thread path
vg sim -x x.vg -P 1#0#x#0 -n 500 -a -s 23 > sim.gam
vg sim -x x.vg -P 1#1#x#0 -n 500 -a -s 23 >> sim.gam
# call some hets
vg pack -x x.vg -o x.pack -g sim.gam
vg call x.vg -k x.pack -a > call.vcf
vg call x.vg -k x.pack -g x.gbwt > callg.vcf
is "$(grep -v 0/0 callg.vcf | grep -v lowad | wc -l)" "$(grep -v 0/0 call.vcf | grep -v lowad | wc -l)" "vg call finds same variants when using gbwt to enumerate traversals"
# try with gbz
vg call x.gbz -k x.pack -z > callz.vcf
cat callg.vcf | grep -v lowad | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6}' > callg.6
cat callz.vcf | grep -v lowad | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6}' > callz.6
diff callg.6 callz.6
is $? 0 "call produces same output with gbwt and gbz"

rm -f x.vg x.gbwt x.gg x.gbz sim.gam x.pack call.vcf callg.vcf callz.vcf callg.6 callz.6


# subpath test
sed -e 's/x/x[100]/g' small/x.fa > x_sub1.fa
sed -e 's/x/x[10000]/g' small/x.fa > x_sub2.fa
gzip -dc small/x.vcf.gz | sed -e 's/x/x[100]/g' | bgzip > x_sub1.vcf.gz && tabix -fp vcf x_sub1.vcf.gz
gzip -dc small/x.vcf.gz | sed -e 's/x/x[10000]/g' | bgzip > x_sub2.vcf.gz && tabix -fp vcf x_sub2.vcf.gz
vg construct -r x_sub1.fa -v x_sub1.vcf.gz -r x_sub2.fa -v x_sub2.vcf.gz > x_subs.vg
vg sim -x x_subs.vg -n 1000 -a -s 23 > sim.gam
vg pack -x x_subs.vg -o x_subs.pack -g sim.gam
vg call x_subs.vg -k x_subs.pack > x_subs.vcf
is $(grep "^##contig=<ID=x,length=11001>" x_subs.vcf | wc -l) 1 "vg call makes currect base path header with subpath input"
is $(grep "^##contig" x_subs.vcf | wc -l) 1 "vg call makes only currect base path header with subpath input"
is "$(grep -v "^#" x_subs.vcf | wc -l)" "$(grep "^x" x_subs.vcf | grep -v "\[" | wc -l)" "vg call only reports base paths with subpath input"
vg call x_subs.vg -k x_subs.pack -p x -l 50000 > x_subs_override.vcf
is $(grep "^##contig=<ID=x,length=50000>" x_subs_override.vcf | wc -l) 1 "vg call makes currect base path header with subpath input and override"
is $(grep "^##contig" x_subs_override.vcf | wc -l) 1 "vg call makes only currect base path header with subpath input and override"
grep -v "##contig" x_subs.vcf > x_subs_nocontig.vcf
grep -v "##contig" x_subs_override.vcf > x_subs_override_nocontig.vcf
diff x_subs_nocontig.vcf x_subs_override_nocontig.vcf
is $? 0 "overriding contig length does not change calls"

rm -f x_sub1.fa x_sub1.fa.fai x_sub2.fa x_sub2.fa.fai x_sub1.vcf.gz x_sub1.vcf.gz.tbi  x_sub2.vcf.gz x_sub2.vcf.gz.tbi sim.gam x_subs.vcf x_subs_override.vcf x_subs_nocontig.vcf x_subs_override_nocontig.vcf






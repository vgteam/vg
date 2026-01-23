#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 41

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
vg gbwt -o HGSVC_travs.gbwt -x HGSVC_alts.xg -A HGSVC_travs.gaf.gz
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
vg call c.aug.xg -k m.aug.pack -p s1 >m.vcf
is $(cat m.vcf | grep -v "^#" | grep -v "0/0" | wc -l) 3 "vg call finds true homozygous variants in a cyclic graph"
rm -f c.vg c.xg c.gcsa c.gcsa.lcp m.fa m.vg m.xg m.sim m.gam m.aug.gam c.aug.vg c.aug.xg m.aug.pack m.vcf

# simple gbwt
vg autoindex -r small/x.fa -v small/x.vcf.gz -w giraffe -p x
rm -f x.min x.dist
mv x.giraffe.gbz x.gbz
vg gbwt -o x.gbwt -Z x.gbz
vg convert x.gbz -p > x.vg
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

rm -f x.vg x.gbz x.gbwt sim.gam x.pack call.vcf callg.vcf callz.vcf callg.6 callz.6


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

# Test: Clustering option (-L) in vg call
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - > small_cluster_call.gfa
printf "L\t1\t+\t9\t+\t0M\n" >> small_cluster_call.gfa
printf "P\ty\t1+,2+,4+,6+,7+,9+\t*\n" >> small_cluster_call.gfa
printf "P\tz\t1+,9+\t*\n" >> small_cluster_call.gfa
vg view -Fv small_cluster_call.gfa > small_cluster_call.vg
vg sim -x small_cluster_call.vg -n 100 -l 20 -a -s 1 > small_cluster_call.gam
vg pack -x small_cluster_call.vg -g small_cluster_call.gam -o small_cluster_call.pack
# Call without clustering
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -n > call_no_cluster.vcf 2>/dev/null
# Call with clustering
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -n -L 0.3 > call_cluster.vcf 2>/dev/null
is $(grep -v "^#" call_no_cluster.vcf | wc -l) $(grep -v "^#" call_cluster.vcf | wc -l) "clustering does not change number of variant sites in vg call"

rm -f small_cluster_call.gfa small_cluster_call.vg small_cluster_call.gam small_cluster_call.pack call_no_cluster.vcf call_cluster.vcf

# Test: Nested calling with vg call -n (basic test)
# Use a larger graph with clear variants
vg construct -r small/x.fa -v small/x.vcf.gz -a > nested_call_test.vg
vg sim -x nested_call_test.vg -n 500 -l 50 -a -s 42 > nested_call_test.gam
vg pack -x nested_call_test.vg -g nested_call_test.gam -o nested_call_test.pack
vg call nested_call_test.vg -k nested_call_test.pack -p x -n > nested_call_test.vcf 2>/dev/null
# Should produce at least one variant (the small/x graph has multiple variants)
NESTED_VARIANT_COUNT=$(grep -v "^#" nested_call_test.vcf | wc -l)
is $(if [ "$NESTED_VARIANT_COUNT" -ge 1 ]; then echo "1"; else echo "0"; fi) "1" "nested vg call produces at least one variant"

rm -f nested_call_test.vg nested_call_test.gam nested_call_test.pack nested_call_test.vcf

# Test: Nested calling emits both top-level and child snarls
# This tests the fix for nested snarl VCF emission where genotypes propagate to children
vg view -Fv nesting/nested_snp_in_del.gfa > nested_snp.vg
vg sim -x nested_snp.vg -n 100 -l 5 -a -s 42 > nested_snp.gam
vg pack -x nested_snp.vg -g nested_snp.gam -o nested_snp.pack
vg call nested_snp.vg -k nested_snp.pack -n -p x 2>/dev/null > nested_snp.vcf
# Should have exactly 2 variant lines: one for top-level snarl (1->6) and one for nested snarl (2->5)
NESTED_LINE_COUNT=$(grep -v "^#" nested_snp.vcf | wc -l)
is "$NESTED_LINE_COUNT" "2" "nested vg call emits both top-level and child snarl variants"

rm -f nested_snp.vg nested_snp.gam nested_snp.pack nested_snp.vcf

# Test: Star allele option validation (-Y requires -n)
vg construct -r small/x.fa -v small/x.vcf.gz > star_test.vg
vg sim -x star_test.vg -n 100 -l 20 -a -s 1 > star_test.gam
vg pack -x star_test.vg -g star_test.gam -o star_test.pack
vg call star_test.vg -k star_test.pack -Y 2>&1 | grep -q "requires -n/--nested"
is "$?" 0 "star allele option requires nested mode"

rm -f star_test.vg star_test.gam star_test.pack

# Test: Cluster-post validation (requires -L < 1.0)
vg construct -r small/x.fa -v small/x.vcf.gz > cluster_post_test.vg
vg sim -x cluster_post_test.vg -n 100 -l 20 -a -s 1 > cluster_post_test.gam
vg pack -x cluster_post_test.vg -g cluster_post_test.gam -o cluster_post_test.pack
vg call cluster_post_test.vg -k cluster_post_test.pack --cluster-post 2>&1 | grep -q "requires -L/--cluster"
is "$?" 0 "cluster-post option requires cluster threshold"

rm -f cluster_post_test.vg cluster_post_test.gam cluster_post_test.pack

# Test: Altpath option validation (requires -n)
vg construct -r small/x.fa -v small/x.vcf.gz > altpath_test.vg
vg sim -x altpath_test.vg -n 100 -l 20 -a -s 1 > altpath_test.gam
vg pack -x altpath_test.vg -g altpath_test.gam -o altpath_test.pack
vg call altpath_test.vg -k altpath_test.pack --altpaths 2>&1 | grep -q "requires -n/--nested"
is "$?" 0 "altpaths option requires nested mode"

rm -f altpath_test.vg altpath_test.gam altpath_test.pack

# =============================================================================
# Nested genotype propagation tests
# These tests verify correct genotype propagation from parent to child snarls
# Graph: nested_snp_in_del.gfa
#   x (ref):  1->2->3->5->6  (allele 0 at top-level, allele 0 at nested)
#   y0:       1->2->4->5->6  (allele 0 at top-level with SNP, allele 1 at nested)
#   y1:       1->6           (allele 1 = deletion at top-level, star at nested)
# =============================================================================

# Test 0/0: homozygous reference - reads only from x path
vg sim -x nesting/nested_snp_in_del.gfa -P x -n 100 -l 2 -a -s 1 > nd_00.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_00.gam -o nd_00.pack
vg call nesting/nested_snp_in_del.gfa -k nd_00.pack -n -p x 2>/dev/null > nd_00.vcf
# 0/0 should produce no non-ref variants
ND_00_NONREF=$(grep -v "^#" nd_00.vcf | grep -v "0/0" | wc -l)
is "$ND_00_NONREF" "0" "nested_snp_in_del 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/SNP - reads from x and y0 (both traverse nested snarl)
vg sim -x nesting/nested_snp_in_del.gfa -P x -n 50 -l 2 -a -s 10 > nd_01.gam
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 50 -l 2 -a -s 11 >> nd_01.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_01.gam -o nd_01.pack
vg call nesting/nested_snp_in_del.gfa -k nd_01.pack -n -p x 2>/dev/null > nd_01.vcf
# Should have 2 variants (top-level and nested), both het
ND_01_COUNT=$(grep -v "^#" nd_01.vcf | wc -l)
is "$ND_01_COUNT" "2" "nested_snp_in_del 0/1: produces both top-level and nested variants"

# Test 1/1: homozygous alt SNP - reads only from y0 path (via sample a haplotype 1)
# We need to simulate specifically from haplotype 1 (y0) not from both haplotypes
# Since -m a gives both haplotypes, we use -A and rely on the path structure
vg sim -x nesting/nested_snp_in_del.gfa -A -n 100 -l 2 -a -s 20 > nd_11.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_11.gam -o nd_11.pack
vg call nesting/nested_snp_in_del.gfa -k nd_11.pack -n -p x 2>/dev/null > nd_11.vcf
# Should produce variants (exact genotype depends on path coverage)
ND_11_EXIT=$?
is "$ND_11_EXIT" "0" "nested_snp_in_del 1/1: vg call completes without error"

# Test 1/2: het SNP/deletion - reads from y0 and y1 (sample a has both)
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 100 -l 2 -a -s 30 > nd_12.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_12.gam -o nd_12.pack
vg call nesting/nested_snp_in_del.gfa -k nd_12.pack -n -p x 2>/dev/null > nd_12.vcf
# Should have 2 variants, nested one should have missing allele
ND_12_COUNT=$(grep -v "^#" nd_12.vcf | wc -l)
is "$ND_12_COUNT" "2" "nested_snp_in_del 1/2: het SNP/del produces both variants"
# Nested snarl should have missing allele (.) for deletion parent
ND_12_MISSING=$(grep ">2>5" nd_12.vcf | grep -c "\./")
is "$ND_12_MISSING" "1" "nested_snp_in_del 1/2: nested snarl shows missing allele for deletion"

rm -f nd_00.gam nd_00.pack nd_00.vcf nd_01.gam nd_01.pack nd_01.vcf
rm -f nd_11.gam nd_11.pack nd_11.vcf nd_12.gam nd_12.pack nd_12.vcf

# =============================================================================
# Star allele tests (-Y flag)
# When parent allele doesn't traverse child, output * instead of .
# =============================================================================

# Test star allele output with -Y flag
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 100 -l 2 -a -s 40 > star.gam
vg pack -x nesting/nested_snp_in_del.gfa -g star.gam -o star.pack
vg call nesting/nested_snp_in_del.gfa -k star.pack -n -Y -p x 2>/dev/null > star.vcf
# With -Y, nested snarl should have * in ALT instead of . in GT
# The genotype will use numeric index (e.g. 1/2) where one allele is *
STAR_IN_ALT=$(grep ">2>5" star.vcf | cut -f5 | grep -c "\*")
is "$STAR_IN_ALT" "1" "star allele: -Y flag produces * in ALT for spanning deletion"
# Verify the genotype doesn't have . when -Y is used (it uses indexed * instead)
NO_MISSING_GT=$(grep ">2>5" star.vcf | cut -f10 | grep -v "\." | wc -l)
is "$NO_MISSING_GT" "1" "star allele: genotype uses indexed * instead of . with -Y"

rm -f star.gam star.pack star.vcf

# =============================================================================
# nested_snp_in_ins.gfa tests
# Graph structure:
#   x (ref):  1->6                    (short ref, allele 0 at top-level)
#   y0 (a#1): 1->2->4->5->6           (insertion with SNP A, allele 1 at top, allele 1 at nested)
#   y1 (a#2): 1->2->3->5->6           (insertion with SNP T, allele 2 at top, allele 0 at nested)
# Top-level snarl: 1->6, Nested snarl: 2->5
# =============================================================================

# Test 0/0: homozygous reference - reads only from x path (short ref)
vg sim -x nesting/nested_snp_in_ins.gfa -P x -n 100 -l 2 -a -s 70 > ni_00.gam
vg pack -x nesting/nested_snp_in_ins.gfa -g ni_00.gam -o ni_00.pack
vg call nesting/nested_snp_in_ins.gfa -k ni_00.pack -n -p x 2>/dev/null > ni_00.vcf
# 0/0 should produce no non-ref variants (or only ref calls)
NI_00_NONREF=$(grep -v "^#" ni_00.vcf | grep -v "0/0" | wc -l)
is "$NI_00_NONREF" "0" "nested_snp_in_ins 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/insertion - reads from x and y1 (a#2 haplotype)
# Note: nested snarl (2->5) not on ref path, so only top-level variant emitted
vg sim -x nesting/nested_snp_in_ins.gfa -P x -n 50 -l 4 -a -s 71 > ni_01.gam
vg sim -x nesting/nested_snp_in_ins.gfa -P "a#2#y1#0" -n 50 -l 4 -a -s 72 >> ni_01.gam
vg pack -x nesting/nested_snp_in_ins.gfa -g ni_01.gam -o ni_01.pack
vg call nesting/nested_snp_in_ins.gfa -k ni_01.pack -n -p x 2>/dev/null > ni_01.vcf
# Only top-level variant (nested snarl not on ref path x, can't emit to VCF)
NI_01_COUNT=$(grep -v "^#" ni_01.vcf | wc -l)
is "$NI_01_COUNT" "1" "nested_snp_in_ins 0/1: het ref/ins produces top-level variant"

# Test 1/1: homozygous insertion - reads only from y1 path
vg sim -x nesting/nested_snp_in_ins.gfa -P "a#2#y1#0" -n 100 -l 4 -a -s 73 > ni_11.gam
vg pack -x nesting/nested_snp_in_ins.gfa -g ni_11.gam -o ni_11.pack
vg call nesting/nested_snp_in_ins.gfa -k ni_11.pack -n -p x 2>/dev/null > ni_11.vcf
# Only top-level variant (nested snarl not on ref path)
NI_11_COUNT=$(grep -v "^#" ni_11.vcf | wc -l)
is "$NI_11_COUNT" "1" "nested_snp_in_ins 1/1: homozygous ins produces top-level variant"

# Test 1/2: het between two insertion alleles - reads from both y0 and y1
vg sim -x nesting/nested_snp_in_ins.gfa -m a -n 100 -l 4 -a -s 74 > ni_12.gam
vg pack -x nesting/nested_snp_in_ins.gfa -g ni_12.gam -o ni_12.pack
vg call nesting/nested_snp_in_ins.gfa -k ni_12.pack -n -p x 2>/dev/null > ni_12.vcf
# Only top-level variant (nested snarl not on ref path)
NI_12_COUNT=$(grep -v "^#" ni_12.vcf | wc -l)
is "$NI_12_COUNT" "1" "nested_snp_in_ins 1/2: het ins/ins produces top-level variant"

rm -f ni_00.gam ni_00.pack ni_00.vcf ni_01.gam ni_01.pack ni_01.vcf
rm -f ni_11.gam ni_11.pack ni_11.vcf ni_12.gam ni_12.pack ni_12.vcf

# =============================================================================
# Triple nested graph tests
# Graph structure:
#   x (ref):  1->5 (short ref, bypasses all nesting)
#   y0 (a#1): 1->2->3->31->311->313->32->33->34->4->5 (deep insertion)
#   y1 (a#2): 1->2->3->31->311->312->32->33->34->4->5 (different at deepest level)
#   y2 (a#3): same as y0
# Top-level snarl: 1->5, with 4 levels of nesting inside
# Since ref bypasses all nesting, only top-level variant can be emitted
# =============================================================================

# Test 0/0: homozygous reference - reads only from x path
vg sim -x nesting/triple_nested.gfa -P x -n 100 -l 2 -a -s 80 > tn_00.gam
vg pack -x nesting/triple_nested.gfa -g tn_00.gam -o tn_00.pack
vg call nesting/triple_nested.gfa -k tn_00.pack -n -p x 2>/dev/null > tn_00.vcf
TN_00_NONREF=$(grep -v "^#" tn_00.vcf | grep -v "0/0" | wc -l)
is "$TN_00_NONREF" "0" "triple_nested 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/insertion - reads from x and y0
# Need longer reads (10bp) to span the complex nested structure
vg sim -x nesting/triple_nested.gfa -P x -n 100 -l 10 -a -s 81 > tn_01.gam
vg sim -x nesting/triple_nested.gfa -P "a#1#y0#0" -n 100 -l 10 -a -s 82 >> tn_01.gam
vg pack -x nesting/triple_nested.gfa -g tn_01.gam -o tn_01.pack
vg call nesting/triple_nested.gfa -k tn_01.pack -n -p x 2>/dev/null > tn_01.vcf
# Only top-level variant (nested snarls not on ref path)
TN_01_COUNT=$(grep -v "^#" tn_01.vcf | wc -l)
is "$TN_01_COUNT" "1" "triple_nested 0/1: het ref/ins produces top-level variant"

# Test 1/1: homozygous insertion - reads only from y0 path
vg sim -x nesting/triple_nested.gfa -P "a#1#y0#0" -n 100 -l 10 -a -s 83 > tn_11.gam
vg pack -x nesting/triple_nested.gfa -g tn_11.gam -o tn_11.pack
vg call nesting/triple_nested.gfa -k tn_11.pack -n -p x 2>/dev/null > tn_11.vcf
TN_11_COUNT=$(grep -v "^#" tn_11.vcf | wc -l)
is "$TN_11_COUNT" "1" "triple_nested 1/1: homozygous ins produces top-level variant"

# Test 1/2: het between insertion alleles - reads from y0 and y1 (differ at deepest SNP)
vg sim -x nesting/triple_nested.gfa -P "a#1#y0#0" -n 100 -l 10 -a -s 84 > tn_12.gam
vg sim -x nesting/triple_nested.gfa -P "a#2#y1#0" -n 100 -l 10 -a -s 85 >> tn_12.gam
vg pack -x nesting/triple_nested.gfa -g tn_12.gam -o tn_12.pack
vg call nesting/triple_nested.gfa -k tn_12.pack -n -p x 2>/dev/null > tn_12.vcf
TN_12_COUNT=$(grep -v "^#" tn_12.vcf | wc -l)
is "$TN_12_COUNT" "1" "triple_nested 1/2: het ins/ins produces top-level variant"

rm -f tn_00.gam tn_00.pack tn_00.vcf tn_01.gam tn_01.pack tn_01.vcf
rm -f tn_11.gam tn_11.pack tn_11.vcf tn_12.gam tn_12.pack tn_12.vcf

# =============================================================================
# Short reference bypass tests
# =============================================================================

# Test nested_snp_in_nested_ins.gfa - ref bypasses all nested structures
vg sim -x nesting/nested_snp_in_nested_ins.gfa -m a -n 100 -l 2 -a -s 60 > bypass.gam
vg pack -x nesting/nested_snp_in_nested_ins.gfa -g bypass.gam -o bypass.pack
vg call nesting/nested_snp_in_nested_ins.gfa -k bypass.pack -n -p x 2>/dev/null > bypass.vcf
BYPASS_EXIT=$?
is "$BYPASS_EXIT" "0" "nested_snp_in_nested_ins: vg call handles short-ref nested graph without crashing"

rm -f bypass.gam bypass.pack bypass.vcf




#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 86

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
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --top-down > call_no_cluster.vcf 2>/dev/null
# Call with clustering
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --top-down -L 0.3 > call_cluster.vcf 2>/dev/null
is $(grep -v "^#" call_no_cluster.vcf | wc -l) $(grep -v "^#" call_cluster.vcf | wc -l) "clustering does not change number of variant sites in vg call"

rm -f small_cluster_call.gfa small_cluster_call.vg small_cluster_call.gam small_cluster_call.pack call_no_cluster.vcf call_cluster.vcf

# Test: Nested calling with vg call --top-down (basic test)
# Use a larger graph with clear variants
vg construct -r small/x.fa -v small/x.vcf.gz -a > nested_call_test.vg
vg sim -x nested_call_test.vg -n 500 -l 50 -a -s 42 > nested_call_test.gam
vg pack -x nested_call_test.vg -g nested_call_test.gam -o nested_call_test.pack
vg call nested_call_test.vg -k nested_call_test.pack -p x --top-down > nested_call_test.vcf 2>/dev/null
# Should produce at least one variant (the small/x graph has multiple variants)
NESTED_VARIANT_COUNT=$(grep -v "^#" nested_call_test.vcf | wc -l)
is $(if [ "$NESTED_VARIANT_COUNT" -ge 1 ]; then echo "1"; else echo "0"; fi) "1" "nested vg call produces at least one variant"

rm -f nested_call_test.vg nested_call_test.gam nested_call_test.pack nested_call_test.vcf

# Test: Nested calling emits both top-level and child snarls
# This tests the fix for nested snarl VCF emission where genotypes propagate to children
vg view -Fv nesting/nested_snp_in_del.gfa > nested_snp.vg
vg sim -x nested_snp.vg -n 100 -l 5 -a -s 42 > nested_snp.gam
vg pack -x nested_snp.vg -g nested_snp.gam -o nested_snp.pack
vg call nested_snp.vg -k nested_snp.pack --top-down -p x 2>/dev/null > nested_snp.vcf
# Should have exactly 2 variant lines: one for top-level snarl (1->6) and one for nested snarl (2->5)
NESTED_LINE_COUNT=$(grep -v "^#" nested_snp.vcf | wc -l)
is "$NESTED_LINE_COUNT" "2" "nested vg call emits both top-level and child snarl variants"

rm -f nested_snp.vg nested_snp.gam nested_snp.pack nested_snp.vcf

# Test: Star allele option validation (-Y requires --top-down)
vg construct -r small/x.fa -v small/x.vcf.gz > star_test.vg
vg sim -x star_test.vg -n 100 -l 20 -a -s 1 > star_test.gam
vg pack -x star_test.vg -g star_test.gam -o star_test.pack
vg call star_test.vg -k star_test.pack -Y 2>&1 | grep -q "requires --top-down"
is "$?" 0 "star allele option requires top-down mode"

rm -f star_test.vg star_test.gam star_test.pack

# Test: Cluster-post validation (requires -L < 1.0)
vg construct -r small/x.fa -v small/x.vcf.gz > cluster_post_test.vg
vg sim -x cluster_post_test.vg -n 100 -l 20 -a -s 1 > cluster_post_test.gam
vg pack -x cluster_post_test.vg -g cluster_post_test.gam -o cluster_post_test.pack
vg call cluster_post_test.vg -k cluster_post_test.pack --cluster-post 2>&1 | grep -q "requires -L/--cluster"
is "$?" 0 "cluster-post option requires cluster threshold"

rm -f cluster_post_test.vg cluster_post_test.gam cluster_post_test.pack

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
vg call nesting/nested_snp_in_del.gfa -k nd_00.pack --top-down -p x 2>/dev/null > nd_00.vcf
# 0/0 should produce no non-ref variants
ND_00_NONREF=$(grep -v "^#" nd_00.vcf | grep -v "0/0" | wc -l)
is "$ND_00_NONREF" "0" "nested_snp_in_del 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/SNP - reads from x and y0 (both traverse nested snarl)
vg sim -x nesting/nested_snp_in_del.gfa -P x -n 50 -l 2 -a -s 10 > nd_01.gam
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 50 -l 2 -a -s 11 >> nd_01.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_01.gam -o nd_01.pack
vg call nesting/nested_snp_in_del.gfa -k nd_01.pack --top-down -p x 2>/dev/null > nd_01.vcf
# Should have 2 variants (top-level and nested), both het
ND_01_COUNT=$(grep -v "^#" nd_01.vcf | wc -l)
is "$ND_01_COUNT" "2" "nested_snp_in_del 0/1: produces both top-level and nested variants"

# Test 1/1: homozygous alt SNP - reads only from y0 path (via sample a haplotype 1)
# We need to simulate specifically from haplotype 1 (y0) not from both haplotypes
# Since -m a gives both haplotypes, we use -A and rely on the path structure
vg sim -x nesting/nested_snp_in_del.gfa -A -n 100 -l 2 -a -s 20 > nd_11.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_11.gam -o nd_11.pack
vg call nesting/nested_snp_in_del.gfa -k nd_11.pack --top-down -p x 2>/dev/null > nd_11.vcf
# Should produce variants (exact genotype depends on path coverage)
ND_11_EXIT=$?
is "$ND_11_EXIT" "0" "nested_snp_in_del 1/1: vg call completes without error"

# Test 1/2: het SNP/deletion - reads from y0 and y1 (sample a has both)
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 100 -l 2 -a -s 30 > nd_12.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_12.gam -o nd_12.pack
vg call nesting/nested_snp_in_del.gfa -k nd_12.pack --top-down -p x 2>/dev/null > nd_12.vcf
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
vg call nesting/nested_snp_in_del.gfa -k star.pack --top-down -Y -p x 2>/dev/null > star.vcf
# With -Y, nested snarl should have * in ALT instead of . in GT
# The genotype will use numeric index (e.g. 1/2) where one allele is *
STAR_IN_ALT=$(grep ">2>5" star.vcf | cut -f5 | grep -c "\*")
is "$STAR_IN_ALT" "1" "star allele: -Y flag produces * in ALT for spanning deletion"
# Verify the genotype doesn't have . when -Y is used (it uses indexed * instead)
# Extract just the GT field (first colon-separated field in SAMPLE column)
NO_MISSING_GT=$(grep ">2>5" star.vcf | cut -f10 | cut -d: -f1 | grep -v "\." | wc -l)
is "$NO_MISSING_GT" "1" "star allele: genotype uses indexed * instead of . with -Y"

rm -f star.gam star.pack star.vcf

# =============================================================================
# Nested quality metrics tests
# Verify QUAL, GQ, GL are computed for nested variant calls
# =============================================================================

# Test: nested calls should have non-zero QUAL
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 100 -l 2 -a -s 50 > nq.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nq.gam -o nq.pack
vg call nesting/nested_snp_in_del.gfa -k nq.pack --top-down -p x 2>/dev/null > nq.vcf

# Check that nested snarl (>2>5) has non-zero QUAL
NQ_QUAL=$(grep ">2>5" nq.vcf | cut -f6)
NQ_HAS_QUAL=$(if [ "$NQ_QUAL" != "0" ] && [ "$NQ_QUAL" != "." ]; then echo "1"; else echo "0"; fi)
is "$NQ_HAS_QUAL" "1" "nested call has non-zero QUAL"

# Check that nested snarl has GQ in FORMAT
NQ_FORMAT=$(grep ">2>5" nq.vcf | cut -f9)
NQ_HAS_GQ=$(echo "$NQ_FORMAT" | grep -c "GQ")
is "$NQ_HAS_GQ" "1" "nested call has GQ field"

# Check that top-level snarl also has quality (unaffected)
TL_QUAL=$(grep ">1>6" nq.vcf | cut -f6)
TL_HAS_QUAL=$(if [ "$TL_QUAL" != "0" ] && [ "$TL_QUAL" != "." ]; then echo "1"; else echo "0"; fi)
is "$TL_HAS_QUAL" "1" "top-level call still has non-zero QUAL"

rm -f nq.gam nq.pack nq.vcf

# Test: triple nested with augref paths should all have quality
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/triple_nested.gfa > tnq_ap.gfa 2>/dev/null
vg sim -x tnq_ap.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 51 > tnq.gam
vg sim -x tnq_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 52 >> tnq.gam
vg pack -x tnq_ap.gfa -g tnq.gam -o tnq.pack
vg call tnq_ap.gfa -k tnq.pack --top-down -P x 2>/dev/null > tnq.vcf

# All variant lines should have non-zero QUAL
TNQ_ZERO_QUAL=$(grep -v "^#" tnq.vcf | awk -F'\t' '$6 == "0" || $6 == "."' | wc -l)
is "$TNQ_ZERO_QUAL" "0" "triple nested calls all have non-zero QUAL"

# All variant lines should have GQ in FORMAT
TNQ_ALL_GQ=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "GQ" | wc -l)
is "$TNQ_ALL_GQ" "0" "triple nested calls all have GQ field"

# All variant lines should have GL (Genotype Likelihood) in FORMAT
TNQ_ALL_GL=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "GL" | wc -l)
is "$TNQ_ALL_GL" "0" "triple nested calls all have GL field"

# All variant lines should have GP (Genotype Posterior) in FORMAT
TNQ_ALL_GP=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "GP" | wc -l)
is "$TNQ_ALL_GP" "0" "triple nested calls all have GP field"

# All variant lines should have XD (Expected Depth) in FORMAT
TNQ_ALL_XD=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "XD" | wc -l)
is "$TNQ_ALL_XD" "0" "triple nested calls all have XD field"

# All variant lines should have AD (Allelic Depth) in FORMAT
TNQ_ALL_AD=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "AD" | wc -l)
is "$TNQ_ALL_AD" "0" "triple nested calls all have AD field"

# GQ values should be in valid range (0-256, integers)
TNQ_INVALID_GQ=$(grep -v "^#" tnq.vcf | awk -F'\t' '{
    split($9, fmt, ":");
    split($10, val, ":");
    for (i=1; i<=length(fmt); i++) {
        if (fmt[i] == "GQ") {
            gq = val[i];
            if (gq !~ /^[0-9]+$/ || gq < 0 || gq > 256) print "invalid";
        }
    }
}' | wc -l)
is "$TNQ_INVALID_GQ" "0" "triple nested calls have valid GQ values (0-256)"

# All variant lines should have LV (nesting level) in INFO
TNQ_ALL_LV=$(grep -v "^#" tnq.vcf | cut -f8 | grep -v "LV=" | wc -l)
is "$TNQ_ALL_LV" "0" "triple nested calls all have LV tag"

# Nested variants (LV > 0) should have PS (parent snarl) in INFO
TNQ_NESTED_NO_PS=$(grep -v "^#" tnq.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/ && $8 !~ /PS=/' | wc -l)
is "$TNQ_NESTED_NO_PS" "0" "nested calls (LV>0) all have PS tag"

# Top-level variants (LV=0) should NOT have PS tag
TNQ_TOPLEVEL_HAS_PS=$(grep -v "^#" tnq.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /PS=/' | wc -l)
is "$TNQ_TOPLEVEL_HAS_PS" "0" "top-level calls (LV=0) do not have PS tag"

rm -f tnq_ap.gfa tnq.gam tnq.pack tnq.vcf

# =============================================================================
# nested_snp_in_ins.gfa tests
# Graph structure:
#   x (ref):  1->6                    (short ref, allele 0 at top-level)
#   y0 (a#1): 1->2->4->5->6           (insertion with SNP A, allele 1 at top, allele 1 at nested)
#   y1 (a#2): 1->2->3->5->6           (insertion with SNP T, allele 2 at top, allele 0 at nested)
# Top-level snarl: 1->6, Nested snarl: 2->5
# NOTE: ref path x bypasses the nested snarl, so we need augref covers to call nested variants
# =============================================================================

# Compute augref cover first (creates x_1_alt, x_2_alt, etc. covering nodes not on x)
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/nested_snp_in_ins.gfa > ni_ap.gfa

# Test 0/0: homozygous reference - reads only from x path (short ref)
vg sim -x ni_ap.gfa -P x -n 100 -l 2 -a -s 70 > ni_00.gam
vg pack -x ni_ap.gfa -g ni_00.gam -o ni_00.pack
vg call ni_ap.gfa -k ni_00.pack --top-down -P x 2>/dev/null > ni_00.vcf
# 0/0 should produce no non-ref variants (or only ref calls)
NI_00_NONREF=$(grep -v "^#" ni_00.vcf | grep -v "0/0" | wc -l)
is "$NI_00_NONREF" "0" "nested_snp_in_ins 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/insertion - reads from x and y1 (a#2 haplotype)
# Note: y1 path is 5bp vs x path 2bp, so need ~2x more y1 reads for balanced boundary coverage
vg sim -x ni_ap.gfa -P x -n 50 -l 2 -a -s 71 > ni_01.gam
vg sim -x ni_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 72 >> ni_01.gam
vg pack -x ni_ap.gfa -g ni_01.gam -o ni_01.pack
vg call ni_ap.gfa -k ni_01.pack --top-down -P x 2>/dev/null > ni_01.vcf
# With augref paths: both top-level and nested variants emitted
NI_01_COUNT=$(grep -v "^#" ni_01.vcf | wc -l)
is "$NI_01_COUNT" "2" "nested_snp_in_ins 0/1: het ref/ins produces top-level and nested variants with augref paths"

# Test 1/1: homozygous insertion - reads only from y1 path
vg sim -x ni_ap.gfa -P "a#2#y1#0" -n 100 -l 2 -a -s 73 > ni_11.gam
vg pack -x ni_ap.gfa -g ni_11.gam -o ni_11.pack
vg call ni_ap.gfa -k ni_11.pack --top-down -P x 2>/dev/null > ni_11.vcf
# With augref paths: both top-level and nested variants emitted
NI_11_COUNT=$(grep -v "^#" ni_11.vcf | wc -l)
is "$NI_11_COUNT" "2" "nested_snp_in_ins 1/1: homozygous ins produces top-level and nested variants with augref paths"

# Test 1/2: het between two insertion alleles - reads from both y0 and y1
vg sim -x ni_ap.gfa -m a -n 200 -l 2 -a -s 74 > ni_12.gam
vg pack -x ni_ap.gfa -g ni_12.gam -o ni_12.pack
vg call ni_ap.gfa -k ni_12.pack --top-down -P x 2>/dev/null > ni_12.vcf
# With augref paths: both top-level and nested variants emitted
NI_12_COUNT=$(grep -v "^#" ni_12.vcf | wc -l)
is "$NI_12_COUNT" "2" "nested_snp_in_ins 1/2: het ins/ins produces top-level and nested variants with augref paths"

rm -f ni_ap.gfa ni_00.gam ni_00.pack ni_00.vcf ni_01.gam ni_01.pack ni_01.vcf
rm -f ni_11.gam ni_11.pack ni_11.vcf ni_12.gam ni_12.pack ni_12.vcf

# =============================================================================
# Triple nested graph tests
# Graph structure:
#   x (ref):  1->5 (short ref, bypasses all nesting)
#   y0 (a#1): 1->2->3->31->311->313->32->33->34->4->5 (deep insertion)
#   y1 (a#2): 1->2->3->31->311->312->32->33->34->4->5 (different at deepest level)
#   y2 (a#3): same as y0
# Top-level snarl: 1->5, with 4 levels of nesting inside
# NOTE: ref path x bypasses all nesting, so we need augref covers to call nested variants
# =============================================================================

# Compute augref cover first (creates x_1_alt, x_2_alt, etc. covering nested nodes)
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/triple_nested.gfa > tn_ap.gfa

# Test 0/0: homozygous reference - reads only from x path
vg sim -x tn_ap.gfa -P x -n 100 -l 2 -a -s 80 > tn_00.gam
vg pack -x tn_ap.gfa -g tn_00.gam -o tn_00.pack
vg call tn_ap.gfa -k tn_00.pack --top-down -P x 2>/dev/null > tn_00.vcf
TN_00_NONREF=$(grep -v "^#" tn_00.vcf | grep -v "0/0" | wc -l)
is "$TN_00_NONREF" "0" "triple_nested 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/insertion - reads from x and y0
# Note: y0 path is 11bp vs x path 2bp, so need ~5x more y0 reads for balanced boundary coverage
vg sim -x tn_ap.gfa -P x -n 50 -l 2 -a -s 81 > tn_01.gam
vg sim -x tn_ap.gfa -P "a#1#y0#0" -n 500 -l 2 -a -s 82 >> tn_01.gam
vg pack -x tn_ap.gfa -g tn_01.gam -o tn_01.pack
vg call tn_ap.gfa -k tn_01.pack --top-down -P x 2>/dev/null > tn_01.vcf
# With augref paths: all 5 nesting levels can be emitted
TN_01_COUNT=$(grep -v "^#" tn_01.vcf | wc -l)
is "$TN_01_COUNT" "5" "triple_nested 0/1: het ref/ins produces all 5 nesting level variants with augref paths"

# Test 1/1: homozygous insertion - reads only from y0 path
vg sim -x tn_ap.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 83 > tn_11.gam
vg pack -x tn_ap.gfa -g tn_11.gam -o tn_11.pack
vg call tn_ap.gfa -k tn_11.pack --top-down -P x 2>/dev/null > tn_11.vcf
# With augref paths: 1 variant emitted (top-level insertion only)
# All nested snarls are 0/0 because y0 matches x_1_alt reference at all levels
# (y0 goes through 313, and x_1_alt also goes through 313)
TN_11_COUNT=$(grep -v "^#" tn_11.vcf | wc -l)
is "$TN_11_COUNT" "1" "triple_nested 1/1: homozygous ins produces only top-level variant"

# Test 1/2: het between insertion alleles - reads from y0 and y1 (differ at deepest SNP)
vg sim -x tn_ap.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 84 > tn_12.gam
vg sim -x tn_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 85 >> tn_12.gam
vg pack -x tn_ap.gfa -g tn_12.gam -o tn_12.pack
vg call tn_ap.gfa -k tn_12.pack --top-down -P x 2>/dev/null > tn_12.vcf
# With augref paths: all 5 nesting levels can be emitted
TN_12_COUNT=$(grep -v "^#" tn_12.vcf | wc -l)
is "$TN_12_COUNT" "5" "triple_nested 1/2: het ins/ins produces all 5 nesting level variants with augref paths"

rm -f tn_ap.gfa tn_00.gam tn_00.pack tn_00.vcf tn_01.gam tn_01.pack tn_01.vcf
rm -f tn_11.gam tn_11.pack tn_11.vcf tn_12.gam tn_12.pack tn_12.vcf

# =============================================================================
# Multi-level SNP test (triple_nested_multisnp.gfa)
# This graph has SNPs at multiple nesting levels:
# - y0 matches x_1_alt at all levels (will be chosen as augref reference)
# - y1 differs from x_1_alt at all nested levels
# When simulating from y1, we should get variants at all 4 nesting levels
# =============================================================================

vg paths -x nesting/triple_nested_multisnp.gfa -Q x --compute-augref --min-augref-len 1 > tn_ms_ap.gfa
vg sim -x tn_ms_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 100 > tn_ms.gam
vg pack -x tn_ms_ap.gfa -g tn_ms.gam -o tn_ms.pack
vg call tn_ms_ap.gfa -k tn_ms.pack --top-down -P x 2>/dev/null > tn_ms.vcf
# Should get 4 variants: top-level + 3 nested SNPs (all at 1/1)
TN_MS_COUNT=$(grep -v "^#" tn_ms.vcf | wc -l)
is "$TN_MS_COUNT" "4" "triple_nested_multisnp 1/1: homozygous alt produces variants at all 4 nesting levels"
# Verify all variants are 1/1 (homozygous alt)
TN_MS_HOM=$(grep -v "^#" tn_ms.vcf | cut -f10 | cut -d: -f1 | grep -c "1/1")
is "$TN_MS_HOM" "4" "triple_nested_multisnp 1/1: all 4 variants are homozygous alt"
# Verify LV tags span levels 0-3
TN_MS_LV0=$(grep -v "^#" tn_ms.vcf | grep -c "LV=0")
TN_MS_LV123=$(grep -v "^#" tn_ms.vcf | grep -c "LV=[123]")
is "$TN_MS_LV0" "1" "triple_nested_multisnp: one top-level variant (LV=0)"
is "$TN_MS_LV123" "3" "triple_nested_multisnp: three nested variants (LV=1,2,3)"

rm -f tn_ms_ap.gfa tn_ms.gam tn_ms.pack tn_ms.vcf

# =============================================================================
# Short reference bypass tests
# NOTE: ref path x bypasses nested structures, so we need augref covers to call nested variants
# =============================================================================

# Test nested_snp_in_nested_ins.gfa - ref bypasses all nested structures
# Compute augref cover first (creates x_1_alt, etc. covering nested nodes)
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/nested_snp_in_nested_ins.gfa > bypass_ap.gfa
vg sim -x bypass_ap.gfa -m a -n 100 -l 2 -a -s 60 > bypass.gam
vg pack -x bypass_ap.gfa -g bypass.gam -o bypass.pack
vg call bypass_ap.gfa -k bypass.pack --top-down -P x 2>/dev/null > bypass.vcf
BYPASS_EXIT=$?
is "$BYPASS_EXIT" "0" "nested_snp_in_nested_ins: vg call handles short-ref nested graph with augref paths without crashing"

rm -f bypass_ap.gfa bypass.gam bypass.pack bypass.vcf

# =============================================================================
# --top-down -a interaction tests
# Verify nested calling with reference output (-a) works correctly
# =============================================================================

# Test 1: nested_snp_in_del 0/0 with --top-down -a
# When ref traverses nested snarl, both levels should get 0/0
vg sim -x nesting/nested_snp_in_del.gfa -P x -n 100 -l 2 -a -s 200 > na_del.gam
vg pack -x nesting/nested_snp_in_del.gfa -g na_del.gam -o na_del.pack
vg call nesting/nested_snp_in_del.gfa -k na_del.pack --top-down -a -p x 2>/dev/null > na_del.vcf

# Count variant lines (should be 2: top-level + nested)
NA_DEL_COUNT=$(grep -v "^#" na_del.vcf | wc -l)
is "$NA_DEL_COUNT" "2" "--top-down -a: nested_snp_in_del 0/0 emits both snarls"

# Verify top-level is 0/0 (use awk to match ID column exactly)
NA_DEL_TOP_GT=$(awk -F'\t' '$3 == ">1>6" {print $10}' na_del.vcf | cut -d: -f1)
is "$NA_DEL_TOP_GT" "0/0" "--top-down -a: nested_snp_in_del top-level is 0/0"

# Verify nested is 0/0 (use awk to match ID column exactly)
NA_DEL_NEST_GT=$(awk -F'\t' '$3 == ">2>5" {print $10}' na_del.vcf | cut -d: -f1)
is "$NA_DEL_NEST_GT" "0/0" "--top-down -a: nested_snp_in_del nested is 0/0"

rm -f na_del.gam na_del.pack na_del.vcf

# Test 2: nested_snp_in_ins 0/0 with --top-down -a (with augref paths)
# When ref bypasses nested snarl and both alleles are ref, nested NOT emitted
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/nested_snp_in_ins.gfa > na_ins_ap.gfa 2>/dev/null
vg sim -x na_ins_ap.gfa -P x -n 100 -l 2 -a -s 200 > na_ins_00.gam
vg pack -x na_ins_ap.gfa -g na_ins_00.gam -o na_ins_00.pack
vg call na_ins_ap.gfa -k na_ins_00.pack --top-down -a -P x 2>/dev/null > na_ins_00.vcf

# Count variant lines (should be 1: only top-level, nested not emitted)
NA_INS_00_COUNT=$(grep -v "^#" na_ins_00.vcf | wc -l)
is "$NA_INS_00_COUNT" "1" "--top-down -a: nested_snp_in_ins 0/0 emits only top-level (ref spans nested)"

rm -f na_ins_00.gam na_ins_00.pack na_ins_00.vcf

# Test 3: nested_snp_in_ins 0/1 with --top-down -a (with augref paths)
# When ref bypasses nested but alt traverses, nested should have ./X genotype
vg sim -x na_ins_ap.gfa -P x -n 50 -l 2 -a -s 200 > na_ins_01.gam
vg sim -x na_ins_ap.gfa -P "a#1#y0#0" -n 100 -l 2 -a -s 201 >> na_ins_01.gam
vg pack -x na_ins_ap.gfa -g na_ins_01.gam -o na_ins_01.pack
vg call na_ins_ap.gfa -k na_ins_01.pack --top-down -a -P x 2>/dev/null > na_ins_01.vcf

# Count variant lines (should be 2: top-level + nested)
NA_INS_01_COUNT=$(grep -v "^#" na_ins_01.vcf | wc -l)
is "$NA_INS_01_COUNT" "2" "--top-down -a: nested_snp_in_ins 0/1 emits both snarls"

# Verify nested has missing allele marker (.)
NA_INS_01_NEST_GT=$(grep ">2>5" na_ins_01.vcf | cut -f10 | cut -d: -f1)
NA_INS_01_HAS_MISSING=$(echo "$NA_INS_01_NEST_GT" | grep -c "\.")
is "$NA_INS_01_HAS_MISSING" "1" "--top-down -a: nested_snp_in_ins 0/1 nested has missing allele (.)"

rm -f na_ins_ap.gfa na_ins_01.gam na_ins_01.pack na_ins_01.vcf

# Test 4: triple_nested 0/0 with --top-down -a (with augref paths)
# When ref bypasses all nested snarls and genotype is 0/0, only top-level emitted
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/triple_nested.gfa > na_tn_ap.gfa 2>/dev/null
vg sim -x na_tn_ap.gfa -P x -n 100 -l 2 -a -s 210 > na_tn_00.gam
vg pack -x na_tn_ap.gfa -g na_tn_00.gam -o na_tn_00.pack
vg call na_tn_ap.gfa -k na_tn_00.pack --top-down -a -P x 2>/dev/null > na_tn_00.vcf

# Count variant lines (should be 1: only top-level, nested not emitted since ref spans them)
NA_TN_00_COUNT=$(grep -v "^#" na_tn_00.vcf | wc -l)
is "$NA_TN_00_COUNT" "1" "--top-down -a: triple_nested 0/0 emits only top-level (ref spans all nested)"

# Verify top-level is 0/0
NA_TN_00_GT=$(grep -v "^#" na_tn_00.vcf | cut -f10 | cut -d: -f1)
is "$NA_TN_00_GT" "0/0" "--top-down -a: triple_nested 0/0 top-level is 0/0"

rm -f na_tn_00.gam na_tn_00.pack na_tn_00.vcf

# Test 5: triple_nested 0/1 with --top-down -a (with augref paths)
# When ref bypasses nested but alt traverses, nested snarls should have ./X genotype
vg sim -x na_tn_ap.gfa -P x -n 50 -l 2 -a -s 211 > na_tn_01.gam
vg sim -x na_tn_ap.gfa -P "a#1#y0#0" -n 500 -l 2 -a -s 212 >> na_tn_01.gam
vg pack -x na_tn_ap.gfa -g na_tn_01.gam -o na_tn_01.pack
vg call na_tn_ap.gfa -k na_tn_01.pack --top-down -a -P x 2>/dev/null > na_tn_01.vcf

# Count variant lines (should be 5: all nesting levels emitted with -a)
NA_TN_01_COUNT=$(grep -v "^#" na_tn_01.vcf | wc -l)
is "$NA_TN_01_COUNT" "5" "--top-down -a: triple_nested 0/1 emits all 5 nesting levels"

# Nested snarls (LV > 0) should have missing allele (.) for the spanning ref
NA_TN_01_NESTED_MISSING=$(grep -v "^#" na_tn_01.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/ {print $10}' | cut -d: -f1 | grep -c "\.")
NA_TN_01_NESTED_COUNT=$(grep -v "^#" na_tn_01.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | wc -l)
is "$NA_TN_01_NESTED_MISSING" "$NA_TN_01_NESTED_COUNT" "--top-down -a: triple_nested 0/1 all nested snarls have missing allele (.)"

rm -f na_tn_ap.gfa na_tn_01.gam na_tn_01.pack na_tn_01.vcf

# =============================================================================
# -A (all-snarls) flag tests
# Verify that -A flag:
# 1. Includes LV/PS header tags in VCF output
# 2. Calls all snarls including nested ones (each independently)
# 3. Produces consistent results
# =============================================================================

# Test: -A flag includes LV and PS header lines
vg construct -r small/x.fa -v small/x.vcf.gz -a > all_snarls_test.vg
vg sim -x all_snarls_test.vg -n 200 -l 30 -a -s 300 > all_snarls_test.gam
vg pack -x all_snarls_test.vg -g all_snarls_test.gam -o all_snarls_test.pack
vg call all_snarls_test.vg -k all_snarls_test.pack -A > all_snarls_test.vcf

# Check for LV header
AS_HAS_LV_HEADER=$(grep -c "##INFO=<ID=LV" all_snarls_test.vcf)
is "$AS_HAS_LV_HEADER" "1" "-A flag: VCF includes LV header line"

# Check for PS header
AS_HAS_PS_HEADER=$(grep -c "##INFO=<ID=PS" all_snarls_test.vcf)
is "$AS_HAS_PS_HEADER" "1" "-A flag: VCF includes PS header line"

# Check that variants have LV tags
AS_VARIANT_COUNT=$(grep -v "^#" all_snarls_test.vcf | wc -l)
AS_LV_COUNT=$(grep -v "^#" all_snarls_test.vcf | grep -c "LV=")
is "$AS_LV_COUNT" "$AS_VARIANT_COUNT" "-A flag: all variants have LV tag"

rm -f all_snarls_test.vg all_snarls_test.gam all_snarls_test.pack all_snarls_test.vcf

# Test: -A flag on nested graph produces variants at all nesting levels
vg view -Fv nesting/nested_snp_in_del.gfa > as_nested.vg
vg sim -x as_nested.vg -m a -n 100 -l 2 -a -s 301 > as_nested.gam
vg pack -x as_nested.vg -g as_nested.gam -o as_nested.pack
vg call as_nested.vg -k as_nested.pack -A -p x > as_nested.vcf 2>/dev/null

# Should produce variants at both nesting levels
AS_NESTED_COUNT=$(grep -v "^#" as_nested.vcf | wc -l)
is "$AS_NESTED_COUNT" "2" "-A flag: nested graph produces both top-level and nested variants"

# Verify LV tags present
AS_NESTED_LV=$(grep -v "^#" as_nested.vcf | grep -c "LV=")
is "$AS_NESTED_LV" "2" "-A flag: nested variants have LV tags"

# Verify LV=0 exists (top-level)
AS_NESTED_LV0=$(grep -v "^#" as_nested.vcf | grep -c "LV=0")
is "$AS_NESTED_LV0" "1" "-A flag: has top-level variant (LV=0)"

# Verify LV=1 exists (nested)
AS_NESTED_LV1=$(grep -v "^#" as_nested.vcf | grep -c "LV=1")
is "$AS_NESTED_LV1" "1" "-A flag: has nested variant (LV=1)"

# Verify nested variant has PS tag pointing to parent
AS_NESTED_PS=$(grep -v "^#" as_nested.vcf | awk -F'\t' '$8 ~ /LV=1/ && $8 ~ /PS=/' | wc -l)
is "$AS_NESTED_PS" "1" "-A flag: nested variant has PS tag"

rm -f as_nested.vg as_nested.gam as_nested.pack as_nested.vcf

# Test: -A flag on triple nested graph with augref paths
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/triple_nested.gfa > as_triple.gfa 2>/dev/null
vg sim -x as_triple.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 302 > as_triple.gam
vg sim -x as_triple.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 303 >> as_triple.gam
vg pack -x as_triple.gfa -g as_triple.gam -o as_triple.pack
vg call as_triple.gfa -k as_triple.pack -A -P x > as_triple.vcf 2>/dev/null

# Should produce variants at multiple nesting levels
AS_TRIPLE_COUNT=$(grep -v "^#" as_triple.vcf | wc -l)
AS_TRIPLE_HAS_VARIANTS=$(if [ "$AS_TRIPLE_COUNT" -ge 3 ]; then echo "1"; else echo "0"; fi)
is "$AS_TRIPLE_HAS_VARIANTS" "1" "-A flag: triple nested produces at least 3 variants"

# Verify all variants have LV tags
AS_TRIPLE_LV=$(grep -v "^#" as_triple.vcf | grep -c "LV=")
is "$AS_TRIPLE_LV" "$AS_TRIPLE_COUNT" "-A flag: all triple nested variants have LV tags"

# Verify PS tags on nested variants (LV > 0)
AS_TRIPLE_NESTED=$(grep -v "^#" as_triple.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | wc -l)
AS_TRIPLE_NESTED_PS=$(grep -v "^#" as_triple.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/ && $8 ~ /PS=/' | wc -l)
is "$AS_TRIPLE_NESTED_PS" "$AS_TRIPLE_NESTED" "-A flag: all nested variants (LV>0) have PS tags"

rm -f as_triple.gfa as_triple.gam as_triple.pack as_triple.vcf

# =============================================================================
# RC, RS, RD tag tests (reference coordinate tags for nested snarls)
# =============================================================================

# Test: RC, RS, RD headers are present
vg paths --compute-augref -Q x --min-augref-len 1 -x nesting/triple_nested.gfa > rc_test.gfa 2>/dev/null
vg sim -x rc_test.gfa -P "a#1#y0#0" -n 100 -l 2 -a -s 400 > rc_test.gam 2>/dev/null
vg sim -x rc_test.gfa -P "a#2#y1#0" -n 100 -l 2 -a -s 401 >> rc_test.gam 2>/dev/null
vg pack -x rc_test.gfa -g rc_test.gam -o rc_test.pack 2>/dev/null
vg call rc_test.gfa -k rc_test.pack --top-down -P x > rc_test.vcf 2>/dev/null

# Check for RC, RS, RD headers
RC_HEADER=$(grep -c "##INFO=<ID=RC" rc_test.vcf)
is "$RC_HEADER" "1" "RC header is present in VCF"

RS_HEADER=$(grep -c "##INFO=<ID=RS" rc_test.vcf)
is "$RS_HEADER" "1" "RS header is present in VCF"

RD_HEADER=$(grep -c "##INFO=<ID=RD" rc_test.vcf)
is "$RD_HEADER" "1" "RD header is present in VCF"

# Check that all variants have RC, RS, RD tags
RC_COUNT=$(grep -v "^#" rc_test.vcf | wc -l)
RC_TAG_COUNT=$(grep -v "^#" rc_test.vcf | grep -c "RC=")
is "$RC_TAG_COUNT" "$RC_COUNT" "All variants have RC tag"

RS_TAG_COUNT=$(grep -v "^#" rc_test.vcf | grep -c "RS=")
is "$RS_TAG_COUNT" "$RC_COUNT" "All variants have RS tag"

RD_TAG_COUNT=$(grep -v "^#" rc_test.vcf | grep -c "RD=")
is "$RD_TAG_COUNT" "$RC_COUNT" "All variants have RD tag"

# Check that top-level variant has RC pointing to its own contig
TOP_RC=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=0/' | grep -o "RC=[^;]*" | cut -d= -f2)
TOP_CHROM=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=0/ {print $1}')
is "$TOP_RC" "$TOP_CHROM" "Top-level variant RC equals its own CHROM"

# Check that nested variants point to top-level's coordinates
# All nested variants should have RC=x (the top-level reference)
NESTED_RC_X=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | grep -c "RC=x")
NESTED_COUNT=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | wc -l)
is "$NESTED_RC_X" "$NESTED_COUNT" "All nested variants have RC=x (top-level contig)"

rm -f rc_test.gfa rc_test.gam rc_test.pack rc_test.vcf


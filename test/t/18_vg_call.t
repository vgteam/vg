#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

# Never index a FORMAT sub-field by position; find it by name in column 9 first. Adding a
# FORMAT field shifts every later one, which broke four assertions here that were not
# testing field order at all -- one of them silently compared BL against a GQ threshold.

plan tests 286

# Toy example of hand-made pileup (and hand inspected truth) to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg view -J -v call/tiny.json > tiny.vg

# With an empty pileup and loci mode we should assert the primery path.
true > empty.gam
vg augment tiny.vg empty.gam -A empty_aug.gam > tiny_aug.vg
vg index tiny_aug.vg -x tiny_aug.xg
vg pack -x tiny_aug.xg -g empty_aug.gam -o tiny_aug.pack
vg call tiny_aug.xg -k tiny_aug.pack > tiny_aug.vcf

is $(grep -v '#' tiny_aug.vcf | wc -l | tr -d ' ') 0 "calling empty gam gives empty VCF"

rm -f tiny.vg tiny_aug.vg tiny_aug.xg empty_aug.gam tiny_aug.pack tiny_aug.vcf empty.gam

# No ref paths test
vg convert --gfa-in graphs/three_samples.gfa > only_haps.vg
vg sim --sample sample2 -n 20 -l 5 --align-out --random-seed 79 --xg-name only_haps.vg > sample2.gam
vg augment only_haps.vg sample2.gam -A sample2.aug.gam > only_haps.aug.vg
vg index only_haps.aug.vg -x only_haps.aug.xg
vg pack -x only_haps.aug.xg -g sample2.aug.gam -o sample2.aug.pack

vg call only_haps.aug.xg -k sample2.aug.pack 2> error.txt
is "$?" 1 "If no path is specified, there must be a reference path to fall back to"
grep "REFERENCE or GENERIC" error.txt
is "$?" 0 "Problem is explained in some detail"
grep "Changing-References" error.txt
is "$?" 0 "Hint towards solution provided"

vg call only_haps.aug.xg -k sample2.aug.pack -p "sample1#1#A#0" > sample2.vcf
is $(fgrep -v "#" sample2.vcf | wc -l | tr -d ' ') 2 "can call against a HAPLOTYPE sense -p path"

vg call only_haps.aug.xg -k sample2.aug.pack -P "sample1#1" > sample2.vcf
is $(fgrep -v "#" sample2.vcf | wc -l | tr -d ' ') 2 "can call against HAPLOTYPE sense -P paths"

vg call only_haps.aug.xg -k sample2.aug.pack -S "sample1" > sample2.vcf
is $(fgrep -v "#" sample2.vcf | wc -l | tr -d ' ') 2 "can call against HAPLOTYPE sense -S paths"

vg call only_haps.aug.xg -k sample2.aug.pack -S "missing" 2> error.txt
is "$?" 1 "-S sample must have usable paths"
grep "REFERENCE or HAPLOTYPE" error.txt
is "$?" 0 "Problem is explained in some detail"
grep "Changing-References" error.txt
is "$?" 0 "Hint towards solution provided"

rm -f only_haps.vg only_haps.xg only_haps.aug.vg only_haps.aug.xg error.txt sample2.gam sample2.vcf sample2.aug.gam sample2.aug.pack

vg construct -r inverting/miniFasta.fa -v inverting/miniFasta_VCFinversion.vcf.gz -S > miniFastaGraph.vg
vg index -x miniFastaGraph.xg -g miniFastaGraph.gcsa miniFastaGraph.vg
vg sim -x miniFastaGraph.xg -n 1000 -l 30 -a -s 1 > miniFasta.gam
vg map -G miniFasta.gam -g miniFastaGraph.gcsa -x miniFastaGraph.xg > miniFastaGraph.gam
vg augment  miniFastaGraph.vg miniFastaGraph.gam -A mappedminitest_aug.gam > mappedminitest_aug.vg
vg index mappedminitest_aug.vg -x mappedminitest_aug.xg
vg pack -x mappedminitest_aug.xg -g mappedminitest_aug.gam -o mappedminitest_aug.pack
vg call  mappedminitest_aug.xg -k mappedminitest_aug.pack > calledminitest.vcf

L_COUNT=$(cat calledminitest.vcf | grep "#" -v | wc -l | tr -d ' ')
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

L_COUNT=$(cat calledminitest.vcf | grep "#" -v | wc -l | tr -d ' ')
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
DIFF_COUNT=$(diff -U1000 <(nl baseline_gts.txt) <(nl gts.txt) | tail -n +4 | grep '^+' | wc -l | tr -d ' ')
LESS_EIGHT=$(if (( $DIFF_COUNT < 8 )); then echo 1; else echo 0; fi)
is "${LESS_EIGHT}" "1" "Fewer than 8 differences between called and true SV genotypes"

# genotype the VCF in haploid mode
vg call HGSVC_alts.xg -k HGSVC_alts.pack -v call/HGSVC_chr22_17200000_17800000.vcf.gz -s HG00514 -d 1 > HGSVC1.vcf
# extract the "true" calls
gzip -dc call/HGSVC_chr22_17200000_17800000.vcf.gz | grep -v '#' | awk '{print $10}' | awk -F ':' '{print $1}' | awk -F '|' '{print $1}'  > baseline_gts1.txt
# extract the called genotypes
grep -v '#' HGSVC1.vcf | sort -k1,1d -k2,2n | awk '{print $10}' | awk -F ':' '{print $1}' | sed 's/\//\|/g' > gts1.txt
DIFF_COUNT=$(diff -U1000 <(nl baseline_gts1.txt) <(nl gts1.txt) | tail -n +4 | grep '^+' | wc -l | tr -d ' ')
LESS_EIGHT=$(if (( $DIFF_COUNT <= 8 )); then echo 1; else echo 0; fi)
is "${LESS_EIGHT}" "1" "Fewer than 8 differences between called haploid and truncated true SV genotypes"

# call all the snarls with -a
vg call HGSVC_alts.xg -k HGSVC_alts.pack -s HG00514 -a > HGSVC2.vcf
REF_COUNT_V=$(grep "0/0" HGSVC.vcf | wc -l | tr -d ' ')
REF_COUNT_A=$(grep "0/0" HGSVC2.vcf | wc -l | tr -d ' ')
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
DIFF_COUNT=$(diff -U1000 <(nl calls-travs.txt) <(nl calls-direct.txt) | tail -n +4 | grep '^+' | wc -l | tr -d ' ')
LESS_THREE=$(if (( $DIFF_COUNT < 3 )); then echo 1; else echo 0; fi)
# because call makes an attempt to call multiple snarls at once when outputting traversals (to make bigger traversals)
# there is some wobble here
is "${LESS_THREE}" "1" "Fewer than 3 differences between allales called via traversals or directly"

rm -f HGSVC_alts.vg HGSVC_alts.xg HGSVC_alts.pack HGSVC.vcf baseline_gts.txt gts.txt HGSVC1.vcf HGSVC2.vcf HGSVC_travs.gaf.gz HGSVC_travs.gbwt HGSVC_travs.vcf HGSVC_direct.vcf baseline_gts1.txt gts1.txt gts-travs.txt gts-direct.txt calls-travs.txt calls-direct.txt

## Read-level genotyping (--read-likelihood)
# Uses the same HGSVC fixture: the GAM is already present, so no new test data.
vg index call/HGSVC_chr22_17119590_17880307.vg -x HGSVC_rl.xg
vg pack -x HGSVC_rl.xg -g call/HGSVC_chr22_17119590_17880307.gam -o HGSVC_rl.pack

# The default path must be untouched by the feature existing.
vg call HGSVC_rl.xg -k HGSVC_rl.pack -t 1 > HGSVC_rl_default.vcf 2>/dev/null
is "$?" "0" "vg call default path still works with read-likelihood support compiled in"

vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 > HGSVC_rl.vcf 2>/dev/null
is "$?" "0" "vg call --read-likelihood runs"

# Address FORMAT fields by NAME, never by position. This block previously hard-coded
# "GT:DP:GL:GQ:GP" and read GL from sub-field 3; adding AD, BL and GQI moved GL to 5 and
# broke three assertions at once, none of which was testing field order.
RL_MISSING=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); need="GT DP GL GQ GP AD BL GQI GQN"; m=split(need,w," "); for(i=1;i<=m;i++){found=0; for(j=1;j<=nk;j++) if(k[j]==w[i]) found=1; if(!found){print; next}}}' | wc -l | tr -d ' ')
is "${RL_MISSING}" "0" "every read-likelihood record carries the full FORMAT field set"

# --min-confidence marks on GQN. A threshold above 1 must catch everything that has a GQN
# at all, since GQN is a fraction -- that is the cheapest way to prove it reads the right
# field and compares in the right direction.
MC_VCF=$(vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --min-confidence 1.1 --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -v "^#")
MC_UNMARKED=$(echo "${MC_VCF}" | awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="GQN") n=j; split($10,f,":"); if (f[n]!="." && $7 !~ /lowconf/) print}' | wc -l | tr -d ' ')
is "${MC_UNMARKED}" "0" "--min-confidence above 1 marks every record carrying a GQN"

# A record with no GQN to measure reports '.', which is not low confidence -- it is no
# measurement. Sweeping those up would filter exactly the sites the model declined to
# judge, which is the opposite of what a confidence threshold should do.
MC_DOTMARKED=$(echo "${MC_VCF}" | awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="GQN") n=j; split($10,f,":"); if (f[n]=="." && $7 ~ /lowconf/) print}' | wc -l | tr -d ' ')
is "${MC_DOTMARKED}" "0" "a record with no GQN is never marked lowconf"

# Off by default, and marks rather than drops.
MC_DEFAULT=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '$7 ~ /lowconf/' | wc -l | tr -d ' ')
is "${MC_DEFAULT}" "0" "lowconf does not fire without --min-confidence"
is "$(echo "${MC_VCF}" | wc -l | tr -d ' ')" "$(grep -vc '^#' HGSVC_rl.vcf)" "--min-confidence marks records rather than dropping them"

# GQN is a fraction of what the site could have achieved, so it is bounded. An
# unbounded value would mean the denominator had gone wrong -- which is exactly the
# failure mode of computing it from the reads' own mismap probabilities rather than
# from the floor, and that version calibrated worse than not normalising at all.
GQN_RANGE=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="GQN") n=j; if(n==0){print; next} split($10,f,":"); if (f[n]==".") next; if (f[n]+0 < 0 || f[n]+0 > 1) print}' | wc -l | tr -d ' ')
is "${GQN_RANGE}" "0" "GQN is always a fraction in [0,1] or missing"

# A confident call must not be near zero on the normalised scale either: if GQN were
# constant, or independent of GQ, it would pass the range check and still be useless.
GQN_FLAT=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="GQN") n=j; split($10,f,":"); if (f[n]!=".") print f[n]}' | sort -u | wc -l | tr -d ' ')
is $(test "${GQN_FLAT}" -gt 1 && echo 1 || echo 0) "1" "GQN varies between records rather than being a constant"

# GL must have exactly one entry per genotype of the emitted alleles, or the
# field is silently mislabelled.
GL_BAD=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); gi=0; for(j=1;j<=nk;j++) if(k[j]=="GL") gi=j; if(gi==0){print; next} n=split($5,alts,","); na=n+1; ng=na*(na+1)/2; split($10,f,":"); ngl=split(f[gi],gl,","); if (ngl != ng) print}' | wc -l | tr -d ' ')
is "${GL_BAD}" "0" "every GL field has the VCF-required number of entries"

# The called genotype must be the one GL says is most likely.
# A tie makes the argmax genuinely ambiguous (a flat likelihood means the reads
# say nothing), so require only that no OTHER genotype is strictly more likely.
GT_MISMATCH=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); gi=0; for(j=1;j<=nk;j++) if(k[j]=="GL") gi=j; if(gi==0){print; next} split($10,f,":"); gt=f[1]; ngl=split(f[gi],gl,","); sub("/","|",gt); split(gt,a,"|"); lo=(a[1]<a[2]?a[1]:a[2]); hi=(a[1]<a[2]?a[2]:a[1]); idx=hi*(hi+1)/2+lo+1; if (idx<1 || idx>ngl) {print; next} for(i=1;i<=ngl;i++) if (gl[i]+0 > gl[idx]+0 + 1e-9) {print; next}}' | wc -l | tr -d ' ')
is "${GT_MISMATCH}" "0" "no genotype is strictly more likely than the one called"

# GQ is GQI scaled by the explained-read fraction, so it can only ever be lower.
# A sign error or an unclamped share above 1 would show up here and nowhere else.
GQ_ABOVE_GQI=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); q=0; qi=0; for(j=1;j<=nk;j++){if(k[j]=="GQ") q=j; if(k[j]=="GQI") qi=j} if(q==0||qi==0){print; next} split($10,f,":"); if (f[q]+0 > f[qi]+0) print}' | wc -l | tr -d ' ')
is "${GQ_ABOVE_GQI}" "0" "GQ never exceeds GQI"

# ...and with the discount off the two must be identical, which is what makes
# --no-share-quality a genuine restoration of the previous behaviour.
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --no-share-quality --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 > HGSVC_rl_nosq.vcf 2>/dev/null
GQ_NE_GQI=$(grep -v "^#" HGSVC_rl_nosq.vcf | awk -F'\t' '{nk=split($9,k,":"); q=0; qi=0; for(j=1;j<=nk;j++){if(k[j]=="GQ") q=j; if(k[j]=="GQI") qi=j} split($10,f,":"); if (f[q] != f[qi]) print}' | wc -l | tr -d ' ')
is "${GQ_NE_GQI}" "0" "--no-share-quality makes GQ equal GQI"

# --ploidy-bed: ploidy per region, which -d and -R cannot express. Tested here because it needs a
# contig with enough sites to have an inside and an outside, and because the read-likelihood path
# is the one that cares.
#
# Mind the two coordinate systems, which this got wrong first: the BED is 0-based half-open and
# VCF POS is 1-based, so [0,E) covers 0-based 0..E-1 -- VCF POS 1..E. A record at POS == E is
# therefore *inside* the window.
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -v "^#" | awk -F'\t' '{split($10,f,":"); print $1"\t"$2"\t"f[1]}' > dip_gts.txt
PB_CONTIG=$(awk 'NR==1{print $1}' dip_gts.txt)
PB_FIRST=$(awk 'NR==1{print $2}' dip_gts.txt)
PB_LAST=$(awk 'END{print $2}' dip_gts.txt)
PB_CUT=$(( (PB_FIRST + PB_LAST) / 2 ))
printf '%s\t0\t%s\t1\n' "${PB_CONTIG}" "${PB_CUT}" > window_ploidy.bed
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --ploidy-bed window_ploidy.bed --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -v "^#" | awk -F'\t' '{split($10,f,":"); print $2"\t"f[1]}' > win_gts.txt
# Both sides must be non-empty, or the two assertions after them hold vacuously.
is $(awk -v c="${PB_CUT}" '$1 <= c' win_gts.txt | wc -l | awk '{print ($1>0)?1:0}') "1" "the haploid window contains at least one call"
is $(awk -v c="${PB_CUT}" '$1 > c' win_gts.txt | wc -l | awk '{print ($1>0)?1:0}') "1" "at least one call falls outside the haploid window"
IN_WINDOW_DIP=$(awk -v c="${PB_CUT}" '$1 <= c && $2 ~ /[\/|]/' win_gts.txt | wc -l | tr -d ' ')
is "${IN_WINDOW_DIP}" "0" "--ploidy-bed makes every call inside a haploid window haploid"
OUT_WINDOW_HAP=$(awk -v c="${PB_CUT}" '$1 > c && $2 !~ /[\/|]/' win_gts.txt | wc -l | tr -d ' ')
is "${OUT_WINDOW_HAP}" "0" "--ploidy-bed leaves calls outside the window at the -d ploidy"
rm -f dip_gts.txt win_gts.txt window_ploidy.bed

# --read-likelihood without reads must fail loudly rather than genotype with none.
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood -t 1 > /dev/null 2> rl_err.txt
is "$?" "1" "--read-likelihood without --gam/--gaf is an error"
is $(grep -c "requires reads" rl_err.txt) "1" "the error explains that reads are required"

rm -f HGSVC_rl.xg HGSVC_rl.pack HGSVC_rl.vcf HGSVC_rl_nosq.vcf HGSVC_rl_default.vcf rl_err.txt

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

vg index -x c.xg -g c.gcsa msgas/c1.vg
# True alignment has 3 variants:
# TCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTC
# TCCCTCCTCAAAGGCTTCTCACTACTCCA-ATCAAAGCTACCCAGGCCATTTTAAGTTTC
#            *       *
cat cyclic/cycle.fa | sed s/TCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTC/TCCCTCCTCAAAGGCTTCTCACTACTCCAATCAAAGCTACCCAGGCCATTTTAAGTTTC/ >m.fa
vg construct -r m.fa >m.vg
vg index -x m.xg m.vg
vg sim -n 200 -s 23823 -l 50 -x m.xg -a >m.sim
vg map -x c.xg -g c.gcsa -G m.sim >m.gam
vg augment msgas/c1.vg m.gam -A m.aug.gam >c.aug.vg
vg index -x c.aug.xg c.aug.vg
vg pack -x c.aug.xg -g m.aug.gam -o m.aug.pack
vg call c.aug.xg -k m.aug.pack -p s1 >m.vcf
is $(cat m.vcf | grep -v "^#" | grep -v "0/0" | wc -l | tr -d ' ') 3 "vg call finds true homozygous variants in a cyclic graph"
rm -f c.xg c.gcsa c.gcsa.lcp m.fa m.vg m.xg m.sim m.gam m.aug.gam c.aug.vg c.aug.xg m.aug.pack m.vcf

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
is "$(grep -v 0/0 callg.vcf | grep -v lowad | wc -l | tr -d ' ')" "$(grep -v 0/0 call.vcf | grep -v lowad | wc -l | tr -d ' ')" "vg call finds same variants when using gbwt to enumerate traversals"
# try with gbz
vg call x.gbz -k x.pack -z > callz.vcf
cat callg.vcf | grep -v lowad | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6}' > callg.6
cat callz.vcf | grep -v lowad | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6}' > callz.6
diff callg.6 callz.6
is $? 0 "call produces same output with gbwt and gbz"

# Read-level genotyping needs no pack file when alleles come from a haplotype
# index: GBWTTraversalFinder enumerates from recorded haplotypes rather than from
# support, so nothing in that path consults the pack.
vg call x.vg --read-likelihood --gam sim.gam -g x.gbwt > callrl_nopack.vcf 2>/dev/null
is "$?" "0" "--read-likelihood works with -g and no pack file"
is $(grep -c "^#CHROM" callrl_nopack.vcf) "1" "pack-free read-likelihood output is a valid VCF"
is $(if [ $(grep -v "^#" callrl_nopack.vcf | wc -l) -gt 0 ]; then echo 1; else echo 0; fi) "1" "pack-free read-likelihood emits variants"

# Same via GBZ.
vg call x.gbz --read-likelihood --gam sim.gam -z > callrl_nopack_z.vcf 2>/dev/null
is "$?" "0" "--read-likelihood works with -z and no pack file"

# Panel enumeration is the default under --read-likelihood on a GBZ that carries haplotypes,
# so -z is redundant there and must be exactly redundant: the same command with and without it
# has to produce the same calls, or the flag and the default disagree about what they select.
vg call x.gbz --read-likelihood --gam sim.gam > rl_autoz_full.vcf 2>/dev/null
is "$?" "0" "--read-likelihood defaults to panel enumeration, so it needs neither -z nor a pack"
grep -v "^#" rl_autoz_full.vcf > rl_autoz.vcf
grep -v "^#" callrl_nopack_z.vcf > rl_explicit_z.vcf
is $(diff rl_explicit_z.vcf rl_autoz.vcf | wc -l | tr -d ' ') "0" \
   "the default enumeration mode is identical to an explicit -z"

# --enumerate-support opts back out. It has to actually change the enumeration, not just be
# accepted: asserted by the pack requirement returning, since support enumeration is the only
# mode that consults one.
vg call x.gbz --read-likelihood --gam sim.gam --enumerate-support >/dev/null 2>es_nopack.txt
is "$?" "1" "--enumerate-support returns the pack file requirement"
vg call x.gbz --read-likelihood --gam sim.gam --enumerate-support -k x.pack > rl_support_full.vcf 2>/dev/null
is "$?" "0" "--enumerate-support works with a pack file"
grep -v "^#" rl_support_full.vcf > rl_support.vcf
is $(if [ $(diff rl_autoz.vcf rl_support.vcf | wc -l | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "--enumerate-support enumerates different alleles from the panel default"

# Two ways of asking for enumeration at once is a contradiction, not a preference order: a
# silent winner would mean the run did not do what the command line says.
vg call x.gbz --read-likelihood --gam sim.gam --enumerate-support -z >/dev/null 2>es_z.txt
is "$?" "1" "--enumerate-support with -z is refused"
vg call x.gbz --read-likelihood --gam sim.gam --enumerate-support -g x.gbwt >/dev/null 2>es_g.txt
is "$?" "1" "--enumerate-support with -g is refused"

# A GBZ always contains a GBWT, but it may hold nothing but reference paths. Enumerating from
# that would offer the reference allele and nothing else -- near-zero alt recall, silently.
# The default has to notice and decline rather than take the empty panel at its word.
vg construct -r small/x.fa -v small/x.vcf.gz > nopanel.vg 2>/dev/null
vg gbwt -E -o nopanel.gbwt -x nopanel.vg --gbz-format -g nopanel.gbz 2>/dev/null
vg call nopanel.gbz --read-likelihood --gam sim.gam >/dev/null 2>nopanel_err.txt
is "$?" "1" "the panel-enumeration default declines on a GBZ with no haplotypes"
is $(grep -c "too few to enumerate alleles from" nopanel_err.txt) "1" \
   "declining on an empty panel says so"
vg pack -x nopanel.vg -o nopanel.pack -g sim.gam 2>/dev/null
vg call nopanel.gbz --read-likelihood --gam sim.gam -k nopanel.pack >/dev/null 2>/dev/null
is "$?" "0" "the empty-panel fallback runs once given the pack file it asked for"

# The support caller keeps its own default. Haplotype enumeration measures worse for it on
# structural variants, so -z must stay opt-in outside --read-likelihood: if the two modes ever
# agreed here, the default had leaked across callers.
vg call x.gbz -k x.pack 2>/dev/null | grep -v "^#" > poisson_default.vcf
vg call x.gbz -k x.pack -z 2>/dev/null | grep -v "^#" > poisson_z.vcf
is $(if [ $(diff poisson_default.vcf poisson_z.vcf | wc -l | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "the support caller does not default to panel enumeration"

# Thread count must not change the output. This is asserted for its own sake and because the
# linkage pass collects sites from parallel threads and re-decides them afterwards: an ordering
# that depended on which thread finished first would make results irreproducible in a way no
# accuracy metric would reveal.
vg call x.gbz --read-likelihood --gam sim.gam -z -t 1 2>/dev/null | grep -v "^#" > rl_t1.vcf
vg call x.gbz --read-likelihood --gam sim.gam -z -t 4 2>/dev/null | grep -v "^#" > rl_t4.vcf
is $(diff rl_t1.vcf rl_t4.vcf | wc -l | tr -d ' ') "0" "--read-likelihood output does not depend on thread count"

vg call x.gbz --read-likelihood --gam sim.gam -z --linkage-weight 1 -t 1 2>/dev/null | grep -v "^#" > rl_link_t1.vcf
vg call x.gbz --read-likelihood --gam sim.gam -z --linkage-weight 1 -t 4 2>/dev/null | grep -v "^#" > rl_link_t4.vcf
is $(diff rl_link_t1.vcf rl_link_t4.vcf | wc -l | tr -d ' ') "0" "--linkage-weight output does not depend on thread count"

# This block used to assert "--linkage-weight 0 is inert" against a run with no flag. That
# reference *was* the default, and the default is now 2, so the same comparison would pit
# linkage-on against linkage-off and pass only if the layer did nothing. The property is not
# expressible against an independent reference any more -- weight 0 is the per-site caller, so
# there is nothing else to compare it with -- and the assertion worth having is the converse: that
# the default is not silently doing nothing.
vg call x.gbz --read-likelihood --gam sim.gam -z --linkage-weight 0 2>/dev/null | grep -v "^#" > rl_link_off.vcf
grep -v "^#" callrl_nopack_z.vcf > rl_link_default.vcf
is $(if [ $(diff rl_link_default.vcf rl_link_off.vcf | wc -l | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" "the default --linkage-weight changes calls that --linkage-weight 0 leaves alone"

# GQ = GQI * share, so GQ can never exceed GQI -- and that has to hold on records the linkage
# pass rewrote too, or one file carries two definitions of the field. It did not: the rewritten
# GQ was the phred complement of the posterior with no discount, and on real data 4.63% of
# records came out with GQ > GQI. Asserted over the linkage-on run because that is the only
# path that can break it.
GQ_OVER_GQI=$(grep -v "^#" callrl_nopack_z.vcf | awk -F'\t' '{
    split($9, k, ":"); split($10, v, ":");
    gq = ""; gqi = "";
    for (i = 1; i <= length(k); i++) { if (k[i] == "GQ") gq = v[i]; if (k[i] == "GQI") gqi = v[i]; }
    if (gq != "" && gqi != "" && gq + 0 > gqi + 0 + 0.001) n++;
} END { print n + 0 }')
is "${GQ_OVER_GQI}" "0" "GQ never exceeds GQI, including on linkage-changed records"

# Weight 0 never constructs the collector, so it is a different code path from a small weight and
# earns its own determinism check.
vg call x.gbz --read-likelihood --gam sim.gam -z --linkage-weight 0 -t 4 2>/dev/null | grep -v "^#" > rl_link_off_t4.vcf
is $(diff rl_link_off.vcf rl_link_off_t4.vcf | wc -l | tr -d ' ') "0" "--linkage-weight 0 output does not depend on thread count"

# Without haplotype enumeration the default must decline in silence, or every support-enumeration
# run would fail on a flag the user never typed.
vg call x.vg --read-likelihood --gam sim.gam -k x.pack 2>/dev/null | grep -v "^#" > rl_nolink_pack.vcf
vg call x.vg --read-likelihood --gam sim.gam -k x.pack --linkage-weight 0 2>/dev/null | grep -v "^#" > rl_nolink_pack0.vcf
is $(diff rl_nolink_pack.vcf rl_nolink_pack0.vcf | wc -l | tr -d ' ') "0" "the default --linkage-weight declines without haplotype enumeration"

# An explicit request that cannot be honoured is an error, though. Silence is right for a default
# and wrong for something the user typed.
vg call x.vg --read-likelihood --gam sim.gam -k x.pack --linkage-weight 1 2>rl_link_err.txt >/dev/null
is "$?" "1" "explicit --linkage-weight without haplotype enumeration is refused"

vg call x.gbz --gam sim.gam -z -k x.pack --linkage-weight 1 2>rl_link_err2.txt >/dev/null
is "$?" "1" "explicit --linkage-weight without --read-likelihood is refused"

# --phased emits the linkage layer's Viterbi path as phased genotypes. The property that makes it
# safe to turn on is that it re-orders a genotype without re-deciding it, so the phased GT must be
# a permutation of the unphased one at every record -- if it is not, phasing has silently changed
# a call, which no amount of downstream checking would catch.
vg call x.gbz --read-likelihood --gam sim.gam --phased 2>/dev/null > rl_phased.vcf
is "$?" "0" "--phased produces output"
is $(grep -c "FORMAT=<ID=PS" rl_phased.vcf) "1" "--phased declares FORMAT/PS in the header"
is $(grep -v "^#" rl_phased.vcf | cut -f10 | cut -d: -f1 | grep -c "/") "0" \
   "every phased record uses | rather than /"
vg call x.gbz --read-likelihood --gam sim.gam 2>/dev/null | grep -v "^#" | cut -f1,2,10 | cut -d: -f1 > rl_unphased_gt.txt
grep -v "^#" rl_phased.vcf | cut -f1,2,10 | cut -d: -f1 > rl_phased_gt.txt
is $(paste rl_unphased_gt.txt rl_phased_gt.txt | awk '{
       split($3,a,"/"); split($6,b,"|");
       if ((a[1]!=b[1] || a[2]!=b[2]) && (a[1]!=b[2] || a[2]!=b[1])) bad++
     } END { print bad+0 }') "0" \
   "the phased genotype is a permutation of the unphased one at every record"

# The mosaic is the same phasing, run-length encoded. Its value rests entirely on being *shorter*
# than the site list -- one line per maximal run rather than per site -- so a mosaic with as many
# segments as sites would mean the encoding had bought nothing and the format is the wrong one.
vg call x.gbz --read-likelihood --gam sim.gam --mosaic-out rl_mosaic.tsv 2>/dev/null >/dev/null
is "$?" "0" "--mosaic-out produces output"
is $(grep -c "^#mosaic-version" rl_mosaic.tsv) "1" "the mosaic declares its version"
is $(if [ $(grep -vc "^#" rl_mosaic.tsv) -lt $(grep -vc "^#" rl_phased.vcf) ]; then echo 1; else echo 0; fi) "1" \
   "the mosaic has fewer segments than there are sites"
# Two strands, so every contig must appear on both -- a mosaic naming only one strand would
# describe a haploid genome.
is $(grep -v "^#" rl_mosaic.tsv | cut -f3 | sort -u | tr -d '\n') "01" "the mosaic covers both strands"
# Segments must name a haplotype the panel actually has, or the file cannot be read back.
is $(grep -v "^#" rl_mosaic.tsv | awk -F'\t' '$9 == "?" {n++} END {print n+0}') "0" \
   "every mosaic segment names a known haplotype"

# --- the mosaic is self-describing ------------------------------------------
# Three properties make the file readable without out-of-band knowledge: it names its reference, it
# names its haplotypes portably, and it hands the consumer a GBWT position. One test each.
is $(awk -F'\t' '$1=="#mosaic-version" {print $2}' rl_mosaic.tsv) "2" "the mosaic declares version 2"

# The set of header keys is the format's contract with a parser, and doc/read-likelihood-genotyping.md
# enumerates it. Adding a key without documenting it is exactly the drift that shipped a header with
# five undocumented #note lines while the doc's example showed none, so pin the set: if this fails,
# update the doc's key list in the same commit.
is "$(grep "^#" rl_mosaic.tsv | cut -f1 | sort -u | paste -sd, -)" \
   "#H,#decoding,#graph,#haplotype,#mosaic-version,#note,#reference,#sample" \
   "the mosaic header uses exactly the documented set of keys"

# (1) ref_start/ref_end need a stated coordinate system: a graph can carry several references, so a
# contig name alone does not fix one.
is $(grep -c "^#reference" rl_mosaic.tsv) "1" "the mosaic names the reference its coordinates use"

# (2) hap_index is a position in this run's GBWT metadata and means nothing outside it, so the file
# also carries the portable sample#phase name and a table binding the two. Every index in the body
# must resolve in that table.
is $(grep -c "^#haplotype" rl_mosaic.tsv | awk '{print ($1>0)?1:0}') "1" \
   "the mosaic carries a haplotype table"
is $(awk -F'\t' '$1=="#haplotype" {name[$2]=$3; next}
     /^H\t/ && $8 != "*" { if (!($8 in name) || name[$8] != $9) bad++ }
     END {print bad+0}' rl_mosaic.tsv) "0" \
   "every segment's hap_index resolves to its named haplotype in the header table"

# (3) The point of the whole format: a consumer holding a segment must be able to find that
# haplotype in the GBWT without a locate query or an r-index, so the file hands them the GBWT
# position outright. Both position columns must be present and last.
is $(awk -F'\t' '$1=="#H" {print NF"/"$(NF-1)"/"$NF}' rl_mosaic.tsv) "12/gbwt_node/gbwt_offset" \
   "the header ends with the GBWT position columns"
is $(awk -F'\t' '/^H\t/ && NF != 12 {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "every segment row carries all 12 columns"

# gbwt_node is an oriented GBWT node, so it encodes start_node as id*2 + is_reverse. If that
# identity fails the position points somewhere other than the segment it is attached to, which is
# the one way this column can be silently useless.
is $(awk -F'\t' '/^H\t/ && $11 != "." && int($11/2) != $6 {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "each GBWT position sits on that segment's own start node"
is $(awk -F'\t' '/^H\t/ && $12 != "." && $12 !~ /^[0-9]+$/ {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "every resolved GBWT offset is a non-negative integer"
is $(awk -F'\t' '/^H\t/ {if (($11==".") != ($12==".")) n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "gbwt_node and gbwt_offset are either both resolved or both absent"
# An unresolvable position is allowed -- it means no panel haplotype of that name visits the node --
# but it must stay the exception, or the column is not doing its job.
is $(awk -F'\t' '/^H\t/ {t++; if ($11 != ".") p++} END {print (p > t/2) ? 1 : 0}' rl_mosaic.tsv) "1" \
   "most segments resolve to a GBWT position"

# A segment carrying one GBWT position is a claim that the position walks the whole segment, which
# fails if the segment crosses a fragment boundary. The caller splits on those boundaries, so
# splitting can only ever add segments, never lose sites: the site total is the invariant to check.
vg call x.gbz --read-likelihood --gam sim.gam --mosaic-out rl_mosaic2.tsv 2>/dev/null >/dev/null
is $(awk -F'\t' '/^H\t/ {n += $10} END {print n+0}' rl_mosaic.tsv) \
   $(awk -F'\t' '/^H\t/ {n += $10} END {print n+0}' rl_mosaic2.tsv) \
   "fragment splitting is deterministic and conserves the site total"

# Haploid chains. chrX outside the pseudoautosomal regions and all of chrY are haploid in a male
# sample, and they used to be dropped from the linkage pass by a guard that only accepted a
# two-allele genotype -- silently costing them both the transition model and any mosaic. The
# mosaic must name one strand, not two: a second would assert a homozygote where there is one copy.
vg call x.gbz --read-likelihood --gam sim.gam -d 1 --phased --mosaic-out rl_hap_mosaic.tsv 2>/dev/null > rl_hap.vcf
is "$?" "0" "--read-likelihood -d 1 runs with phasing and a mosaic"
is $(grep -v "^#" rl_hap.vcf | cut -f10 | cut -d: -f1 | grep -cE "[|/]") "0" \
   "a haploid genotype is a single allele, not a pair"
is $(grep -v "^#" rl_hap_mosaic.tsv | cut -f3 | sort -u | tr -d '\n') "0" \
   "a haploid mosaic names one strand"
is $(if [ $(grep -vc "^#" rl_hap_mosaic.tsv) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "a haploid chain still produces mosaic segments"
is $(grep -v "^#" rl_hap_mosaic.tsv | awk -F'\t' '$9 == "?" {n++} END {print n+0}') "0" \
   "every haploid mosaic segment names a known haplotype"

# A haploid record's GT is a bare allele, and apply_linkage_change used to build the genotype it
# expected as "i/j" regardless. The guard therefore rejected every haploid change: the linkage
# layer did the work, said so in the progress line, and dropped all of it -- so chrY and
# non-pseudoautosomal chrX got no linkage correction, and the mosaic (built from the *post*-linkage
# genotypes) described genotypes the VCF did not contain.
#
# Stated as an implication so it can never pass vacuously: if the layer reports changes, the output
# must differ from the same run with the layer off.
vg call x.gbz --read-likelihood --gam sim.gam -d 1 --linkage-weight 0 2>/dev/null | grep -v "^#" > rl_hap_lw0.vcf
vg call x.gbz --read-likelihood --gam sim.gam -d 1 --linkage-weight 8 --progress 2>rl_hap_lw8.err | grep -v "^#" > rl_hap_lw8.vcf
HAP_CHANGED=$(grep -o '[0-9]* genotypes changed' rl_hap_lw8.err | awk '{print $1}')
HAP_CHANGED=${HAP_CHANGED:-0}
is $(if [ "${HAP_CHANGED}" -eq 0 ] || ! cmp -s rl_hap_lw0.vcf rl_hap_lw8.vcf; then echo 1; else echo 0; fi) "1" \
   "haploid linkage changes reach the VCF rather than being dropped by the genotype guard"
rm -f rl_hap_lw0.vcf rl_hap_lw8.vcf rl_hap_lw8.err

# Phasing is the linkage layer's path, so asking for it with the layer switched off cannot be
# quietly satisfied by writing unphased output.
vg call x.gbz --read-likelihood --gam sim.gam --phased --linkage-weight 0 >/dev/null 2>rl_ph_err.txt
is "$?" "1" "--phased with --linkage-weight 0 is refused"

# PS is per chain, so within one contig every phased record must share a phase set -- that is the
# claim the switch-error benchmark tests, and it should be visible in the output rather than
# implied.
is $(grep -v "^#" rl_phased.vcf | grep -o "PS" | head -1 | wc -l | tr -d ' ') "1" \
   "phased records carry a PS field"
is $(grep -v "^#" rl_phased.vcf | awk -F'\t' '{n=split($9,k,":"); for(i=1;i<=n;i++) if(k[i]=="PS"){split($10,v,":"); print v[i]}}' | sort -u | wc -l | tr -d ' ') "1" \
   "one phase set per contig"

# Every read-likelihood option is meaningless without --read-likelihood, and silently ignoring one
# means the run did not do what the command line says. --linkage-weight was refused while
# --depth-term and the rest were accepted and dropped; this pins the consistent behaviour.
vg call x.vg -k x.pack --depth-term 0.5 2>/dev/null >/dev/null
is "$?" "1" "--depth-term without --read-likelihood is refused"

vg call x.vg -k x.pack --flat-mixture 2>/dev/null >/dev/null
is "$?" "1" "--flat-mixture without --read-likelihood is refused"

# The value-taking form too, since the check reads argv and has to handle --opt=value.
vg call x.vg -k x.pack --mismap-min=0.4 2>/dev/null >/dev/null
is "$?" "1" "--mismap-min=VALUE without --read-likelihood is refused"

# And it must not fire when the mode *is* on: this run differs only in --read-likelihood.
vg call x.vg --read-likelihood --gam sim.gam -k x.pack --depth-term 0.5 2>/dev/null | grep -v "^#" > rl_dt.vcf
is $(if [ $(wc -l < rl_dt.vcf | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" "read-likelihood options are accepted with --read-likelihood"

rm -f rl_dt.vcf rl_t1.vcf rl_t4.vcf rl_link_t1.vcf rl_link_t4.vcf rl_link_off.vcf rl_link_off_t4.vcf rl_link_default.vcf rl_nolink_pack.vcf rl_nolink_pack0.vcf rl_link_err.txt rl_link_err2.txt

# The real assertion: since nothing on the GBWT path consults support, supplying a
# pack file must make no difference whatsoever to the calls.
vg call x.vg -k x.pack --read-likelihood --gam sim.gam -g x.gbwt > callrl_withpack.vcf 2>/dev/null
diff <(grep -v "^#" callrl_withpack.vcf) <(grep -v "^#" callrl_nopack.vcf) > /dev/null
is "$?" "0" "pack-free read-likelihood calls are identical to those made with a pack file"

# But the flow traversal finder is driven entirely by node/edge weights, so
# dropping the pack file there has to be refused rather than silently genotyped
# against zero support.
vg call x.vg --read-likelihood --gam sim.gam > /dev/null 2> nopack_err.txt
is $(grep -c "requires haplotype-based allele enumeration" nopack_err.txt) "1" "--read-likelihood without -k and without -g/-z is refused"

# Indexed GAM read source: reads fetched per site from a .gai instead of all held
# in memory. The calls must be identical -- the backend is a pure substitution.
vg gamsort -i sim.sorted.gam.gai sim.gam > sim.sorted.gam 2>/dev/null
is "$?" "0" "vg gamsort produces a sorted GAM and .gai index"

vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam -t 1 2>/dev/null > rl_inmem.vcf
vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --gam-index sim.sorted.gam.gai -t 1 2>/dev/null > rl_indexed.vcf
is "$?" "0" "--gam-index runs"
diff <(grep -v "^#" rl_inmem.vcf) <(grep -v "^#" rl_indexed.vcf) > /dev/null
is "$?" "0" "indexed GAM read source produces identical calls to the in-memory one"

# --gam-index without --gam is an error rather than a silently ignored flag.
vg call x.vg -k x.pack --gam-index sim.sorted.gam.gai -t 1 >/dev/null 2>gi_err.txt
is $(grep -c "requires --gam" gi_err.txt) "1" "--gam-index without --gam is refused"

# GAF-Base read source. The flag validation and the missing-binary path are checked
# unconditionally, since none of them ever runs gbz-base; the rest needs it installed.
vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --gbz-base x.gbz -t 1 >/dev/null 2>gb_err.txt
is $(grep -c "requires --gaf-base" gb_err.txt) "1" "--gbz-base without --gaf-base is refused"

vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --gaf-base sim.gaf.db -t 1 >/dev/null 2>gb_excl.txt
is $(grep -c "mutually exclusive" gb_excl.txt) "1" "--gaf-base and --gam together are refused"

vg call x.vg -k x.pack --gaf-base sim.gaf.db -t 1 >/dev/null 2>gb_norl.txt
is $(grep -c "only applies to --read-likelihood" gb_norl.txt) "1" "--gaf-base without --read-likelihood is refused"

# A missing gbz-base is the user's setup, not a vg bug: it must be an error with a fix
# in it, not a crash telling them to file an issue.
vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db --gaf-base-binary /nonexistent/gbz-base -t 1 >/dev/null 2>gb_nobin.txt
is $(grep -c "could not execute" gb_nobin.txt) "1" "a missing gbz-base binary is reported with a remedy"
is $(grep -c "VG has crashed" gb_nobin.txt) "0" "a missing gbz-base binary is not reported as a vg crash"

# The equivalence test: the same reads, reached through a database instead of held in
# memory, must give the same calls. Compared against the in-memory *GAF* rather than the
# GAM, because GAM->GAF is not itself lossless -- an insertion at a node boundary can be
# attributed to either side, and the two formats disagree -- so comparing against the GAM
# would be testing vg convert rather than this backend.
GAFBASE_NOTE=""
if ! command -v gbz-base >/dev/null 2>&1 || ! command -v gaf-base >/dev/null 2>&1; then
    GAFBASE_NOTE=" (skipped: gbz-base/gaf-base not on PATH)"
fi

if [ -z "$GAFBASE_NOTE" ]; then
    vg convert -G sim.sorted.gam x.gbz 2>/dev/null | grep -v "^@" > sim.gaf
    gaf-base construct sim.gaf -r x.gbz -o sim.gaf.db --overwrite >/dev/null 2>&1
    gbz-base construct x.gbz -o x.gbz.db >/dev/null 2>&1

    vg call x.vg -k x.pack --read-likelihood --gaf-reads sim.gaf -t 1 2>/dev/null > rl_gafmem.vcf
    vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db --gbz-base x.gbz.db -t 1 2>/dev/null > rl_gafbase.vcf
    GAFBASE_RAN="$?"
    diff <(grep -v "^#" rl_gafmem.vcf) <(grep -v "^#" rl_gafbase.vcf) >/dev/null
    GAFBASE_SAME="$?"

    # Several threads, each with its own cache and output file.
    vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db --gbz-base x.gbz.db -t 4 2>/dev/null > rl_gafbase_t4.vcf
    diff <(grep -v "^#" rl_gafbase.vcf) <(grep -v "^#" rl_gafbase_t4.vcf) >/dev/null
    GAFBASE_THREADS="$?"

    # The window is a performance knob and must not change what is called.
    vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db --gbz-base x.gbz.db --read-window 32 -t 1 2>/dev/null > rl_gafbase_w32.vcf
    diff <(grep -v "^#" rl_gafbase.vcf) <(grep -v "^#" rl_gafbase_w32.vcf) >/dev/null
    GAFBASE_WINDOW="$?"
else
    GAFBASE_RAN="0"
    GAFBASE_SAME="0"
    GAFBASE_THREADS="0"
    GAFBASE_WINDOW="0"
fi

is "$GAFBASE_RAN" "0" "--gaf-base runs$GAFBASE_NOTE"
is "$GAFBASE_SAME" "0" "GAF-Base read source produces identical calls to the in-memory one$GAFBASE_NOTE"
is "$GAFBASE_THREADS" "0" "GAF-Base calls do not depend on the thread count$GAFBASE_NOTE"
is "$GAFBASE_WINDOW" "0" "GAF-Base calls do not depend on the read window size$GAFBASE_NOTE"


# --- nested calling: symbolic alleles and ploidy propagation ------------------
# The design's whole case in one graph, and it fails on a build without --nested.
#
# A top-level snarl (1,6) carries a deletion edge 1->6, and inside it a chain of two nested snarls
# each carrying a single-base SNP. The sample is heterozygous for both SNPs and carries no deletion,
# so its two haplotypes cross the chain and differ only inside it.
#
# By default the two called traversals differ at two separated bases, so flattening cannot reduce
# them and the caller emits one 22 bp substitution -- the shape that buries 55,222 SNVs genome-wide,
# where 90.6% of same-length structural false positives differ at ten bases or fewer. With --nested
# the traversals are symbolically equal, the top-level record collapses to reference, and each
# nested snarl is called at the propagated ploidy and emits its SNP.
rm -f nest.gfa nest.gbz nest.gam nest_default.vcf nest_nested.vcf
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTA\n'
  printf 'S\t2\tGGCATTAGCCTTAGACCGAT\n'
  printf 'S\t3\tA\nS\t4\tT\n'
  printf 'S\t5\tTTGACCAGTTCAGGACTTAC\n'
  printf 'S\t7\tC\nS\t8\tG\n'
  printf 'S\t9\tAACCGGTTACGTTGCAATCG\n'
  printf 'S\t6\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\n'
  printf 'L\t5\t+\t7\t+\t0M\nL\t5\t+\t8\t+\t0M\n'
  printf 'L\t7\t+\t9\t+\t0M\nL\t8\t+\t9\t+\t0M\nL\t9\t+\t6\t+\t0M\n'
  printf 'L\t1\t+\t6\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+,5+,7+,9+,6+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t122\t>1>2>3>5>7>9>6\n'
  printf 'W\tp1\t1\tchr1\t0\t122\t>1>2>4>5>8>9>6\n'
  printf 'W\tp2\t0\tchr1\t0\t60\t>1>6\n'
  printf 'W\tp2\t1\tchr1\t0\t122\t>1>2>4>5>7>9>6\n'
  printf 'W\tp3\t0\tchr1\t0\t122\t>1>2>3>5>8>9>6\n'
  printf 'W\tp3\t1\tchr1\t0\t60\t>1>6\n'
} > nest.gfa
vg gbwt -G nest.gfa --gbz-format -g nest.gbz --set-reference GRCh 2>/dev/null
is "$?" 0 "nested-calling test graph builds"
vg sim -x nest.gbz -n 300 -l 40 -a -s 17 --path "p1#0#chr1#0" > nest.gam 2>/dev/null
vg sim -x nest.gbz -n 300 -l 40 -a -s 23 --path "p1#1#chr1#0" >> nest.gam 2>/dev/null
vg call nest.gbz --read-likelihood --gam nest.gam -t 1 -s samp 2>/dev/null > nest_default.vcf
vg call nest.gbz --read-likelihood --gam nest.gam -t 1 -s samp --nested 2>/dev/null > nest_nested.vcf

is $(grep -vc "^#" nest_default.vcf) "1" \
   "by default the two nested SNPs are buried in a single record"
is $(grep -v "^#" nest_default.vcf | awk '{print length($4)}') "22" \
   "and that record is a long substitution spanning both of them"
is $(grep -vc "^#" nest_nested.vcf) "2" \
   "--nested emits one record per nested SNP instead"
is $(grep -v "^#" nest_nested.vcf | awk 'length($4) == 1 && length($5) == 1' | wc -l | tr -d ' ') "2" \
   "--nested emits them as single-base records, not as one compensating substitution"


# A nested haploid site must not fragment the phase block. Reads from the deletion-bearing haplotype
# and from one that crosses the chain make the parent heterozygous for the deletion, so only one
# parent allele reaches the nested snarls and they are called at ploidy 1. Those sites are a ploidy
# change, and treating them like a regional one -- chrX's pseudoautosomal boundary -- cut the chain at
# every one: chr20's autosomal phasing went from 22 blocks to 9,460 and block N50 from 248 Mb to
# 1.08 Mb. Switch error looked flat only because short blocks make it cheap.
rm -f nest_hap.gam nest_hap.vcf
vg sim -x nest.gbz -n 300 -l 40 -a -s 31 --path "p2#0#chr1#0" > nest_hap.gam 2>/dev/null
vg sim -x nest.gbz -n 300 -l 40 -a -s 37 --path "p2#1#chr1#0" >> nest_hap.gam 2>/dev/null
vg call nest.gbz --read-likelihood --gam nest_hap.gam -t 1 -s samp --nested --phased 2>/dev/null > nest_hap.vcf
is "$?" 0 "--nested --phased runs with a heterozygous deletion over nested sites"
# A nested site reached by only one parent allele is called haploid -- a single allele, not a pair.
# Only one of the two nested snarls emits here: on the surviving haplotype the other carries the
# reference allele, so it is hom-ref and correctly writes nothing.
#
# Bare, without a strand, and that is the point of this fixture rather than an accident. The parent
# is called hom for its deletion, which deletes the chain on both haplotypes, so the site is
# `nested_unreachable` and there is no strand to name. A coherent nested site is written as "a|." or
# ".|a" instead; this asserts the caller does not invent one where the parent says there is none.
is $(grep -v "^#" nest_hap.vcf | awk -F'\t' '{split($10,a,":"); print a[1]}' | grep -cvE "[|/]" ) \
   "1" "a nested site reached by one parent allele is called haploid"
is $(grep -v "^#" nest_hap.vcf | awk -F'\t' '$7 ~ /nested_unreachable/ {n++} END {print n+0}') "1" \
   "a nested site its parent's final genotype deletes on both haplotypes is flagged"
is $(grep -v "^#" nest_hap.vcf | awk -F'\t' '{split($10,a,":"); if (a[1] ~ /\|\.$|^\.\|/) n++} END {print n+0}') \
   "0" "no strand is claimed for a nested site that has none"
# And the parent still reports its deletion: symbolic collapsing must not swallow a real event.
is $(grep -v "^#" nest_hap.vcf | awk -F'\t' 'length($4) > 50 {n++} END {print n+0}') "1" \
   "the parent deletion is still emitted alongside the nested call"
# One phase set across every phased record: the nested haploid sites join the parent's block rather
# than starting their own.
is $(grep -v "^#" nest_hap.vcf | awk -F'\t' '{
       n=split($9,f,":"); for (i=1;i<=n;i++) if (f[i]=="PS") k=i;
       if (k) { split($10,a,":"); if (a[k] != ".") print a[k] }
     }' | sort -u | wc -l | tr -d ' ') "1" \
   "nested haploid sites share one phase set with their parent, not one each"

# Strand assignment has to reach the output, because it is what step three groups on: nested sites
# hanging off opposite parent strands are not on the same haplotype, and chaining them together
# would link sequences that never co-occur. The mosaic is where the per-strand assignment is visible.
rm -f nest_hap.mosaic.tsv
vg call nest.gbz --read-likelihood --gam nest_hap.gam -t 1 -s samp --nested --phased \
    --mosaic-out nest_hap.mosaic.tsv >/dev/null 2>/dev/null
is $(awk -F'\t' '/^H\t/ {print $3}' nest_hap.mosaic.tsv | sort -u | tr -d '\n') "01" \
   "the mosaic carries both strands when nested sites are assigned to a parent strand"
# And every record still reaches the mosaic: a nested site dropped here would break the invariant
# that the mosaic accounts for the whole call set.
is $(awk -F'\t' '/^H\t/ && $3=="0" {n += $10} END {print n+0}' nest_hap.mosaic.tsv) \
   $(grep -vc "^#" nest_hap.vcf) \
   "the mosaic still accounts for every emitted record with nested sites present"
# Segment rows must be in reference order within a strand. A segment is a maximal run of one
# haplotype, so a consumer walks it from start_node to end_node -- and a run assembled out of order
# spans everything between two unrelated loci. Nested sites are what break this: placing one needs
# its parent already phased, so they are phased after every chain and appended. On chr20 that put
# five segments across tens of megabases, one claiming 284 sites between 451 kb and 65.5 Mb, while
# the site totals the check above asserts still added up exactly.
is $(awk -F'\t' '/^H\t/ {if ($3 in prev && $4 < prev[$3]) bad++; prev[$3]=$4} END {print bad+0}' \
     nest_hap.mosaic.tsv) "0" \
   "mosaic segments are in reference order within each strand"

rm -f nest.gfa nest.gbz nest.gam nest_default.vcf nest_nested.vcf nest_hap.gam nest_hap.vcf nest_hap.mosaic.tsv x.vg x.gbz x.gbwt sim.gam x.pack call.vcf callg.vcf callz.vcf callg.6 callz.6 callrl_nopack.vcf callrl_nopack_z.vcf callrl_withpack.vcf nopack_err.txt sim.sorted.gam sim.sorted.gam.gai rl_inmem.vcf rl_indexed.vcf gi_err.txt gb_err.txt gb_excl.txt gb_norl.txt gb_nobin.txt sim.gaf sim.gaf.db x.gbz.db rl_gafmem.vcf rl_gafbase.vcf rl_gafbase_t4.vcf rl_gafbase_w32.vcf rl_autoz_full.vcf rl_autoz.vcf rl_explicit_z.vcf rl_support_full.vcf rl_support.vcf es_nopack.txt es_z.txt es_g.txt nopanel.vg nopanel.gbwt nopanel.gbz nopanel.pack nopanel_err.txt poisson_default.vcf poisson_z.vcf rl_phased.vcf rl_unphased_gt.txt rl_phased_gt.txt rl_ph_err.txt rl_mosaic.tsv rl_mosaic2.tsv rl_hap.vcf rl_hap_mosaic.tsv


# subpath test
sed -e 's/x/x[100]/g' small/x.fa > x_sub1.fa
sed -e 's/x/x[10000]/g' small/x.fa > x_sub2.fa
gzip -dc small/x.vcf.gz | sed -e 's/x/x[100]/g' | bgzip > x_sub1.vcf.gz && tabix -fp vcf x_sub1.vcf.gz
gzip -dc small/x.vcf.gz | sed -e 's/x/x[10000]/g' | bgzip > x_sub2.vcf.gz && tabix -fp vcf x_sub2.vcf.gz
vg construct -r x_sub1.fa -v x_sub1.vcf.gz -r x_sub2.fa -v x_sub2.vcf.gz > x_subs.vg
vg sim -x x_subs.vg -n 1000 -a -s 23 > sim.gam
vg pack -x x_subs.vg -o x_subs.pack -g sim.gam
vg call x_subs.vg -k x_subs.pack > x_subs.vcf
is $(grep "^##contig=<ID=x,length=11001>" x_subs.vcf | wc -l | tr -d ' ') 1 "vg call makes currect base path header with subpath input"
is $(grep "^##contig" x_subs.vcf | wc -l | tr -d ' ') 1 "vg call makes only currect base path header with subpath input"
is "$(grep -v "^#" x_subs.vcf | wc -l | tr -d ' ')" "$(grep "^x" x_subs.vcf | grep -v "\[" | wc -l | tr -d ' ')" "vg call only reports base paths with subpath input"
vg call x_subs.vg -k x_subs.pack -p x -l 50000 > x_subs_override.vcf
is $(grep "^##contig=<ID=x,length=50000>" x_subs_override.vcf | wc -l | tr -d ' ') 1 "vg call makes currect base path header with subpath input and override"
is $(grep "^##contig" x_subs_override.vcf | wc -l | tr -d ' ') 1 "vg call makes only currect base path header with subpath input and override"
grep -v "##contig" x_subs.vcf > x_subs_nocontig.vcf
grep -v "##contig" x_subs_override.vcf > x_subs_override_nocontig.vcf
diff x_subs_nocontig.vcf x_subs_override_nocontig.vcf
is $? 0 "overriding contig length does not change calls"

rm -f x_sub1.fa x_sub1.fa.fai x_sub2.fa x_sub2.fa.fai x_sub1.vcf.gz x_sub1.vcf.gz.tbi  x_sub2.vcf.gz x_sub2.vcf.gz.tbi sim.gam x_subs.vcf x_subs_override.vcf x_subs_nocontig.vcf x_subs_override_nocontig.vcf

# Test: -L/--cluster merges near-identical called alt alleles after genotyping, so a 1/2 call of two
# effectively-identical alleles becomes 1/1.  Two alt paths differing by one 1bp node; reads are
# simulated from each separately so both alleles get support and the site really is called 1/2.
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - > small_cluster_call.gfa
printf "L\t1\t+\t9\t+\t0M\n" >> small_cluster_call.gfa
printf "P\tyA\t1+,2+,4+,6+,7+,9+\t*\n" >> small_cluster_call.gfa
printf "P\tyB\t1+,2+,4+,6+,8+,9+\t*\n" >> small_cluster_call.gfa
vg view -Fv small_cluster_call.gfa > small_cluster_call.vg
vg sim -x small_cluster_call.vg -P yA -n 40 -l 20 -a -s 1 > call_a.gam
vg sim -x small_cluster_call.vg -P yB -n 40 -l 20 -a -s 2 > call_b.gam
cat call_a.gam call_b.gam > small_cluster_call.gam
vg pack -x small_cluster_call.vg -g small_cluster_call.gam -o small_cluster_call.pack

call_site() { awk -F'\t' -v c="$2" '$2 == 9 {print $c}' "$1"; }

vg call small_cluster_call.vg -k small_cluster_call.pack -p x > call_no_cluster.vcf 2>/dev/null
is "$(call_site call_no_cluster.vcf 5)"  "ATTTGA,ATTTGG" "vg call emits two near-identical alts without -L"
is "$(call_site call_no_cluster.vcf 10 | cut -f1 -d:)" "1/2" "vg call genotypes them 1/2 without -L"

# -L 1.0 is the default, so this pins that the merge is off there and that the MAT header line the
# feature adds does not appear.  (The diff alone would be a tautology: there is no "was -L given"
# flag, so -L 1.0 and no -L are the same internal state.)
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 1.0 > call_cluster_off.vcf 2>/dev/null
diff call_no_cluster.vcf call_cluster_off.vcf
is "$?" 0 "-L 1.0 is byte-identical to no -L"
is "$(grep -c 'ID=MAT' call_cluster_off.vcf)" "0" "-L 1.0 does not declare MAT in the header"

vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 --cluster-min-len 0 > call_cluster.vcf 2>/dev/null
is "$(call_site call_cluster.vcf 5)" "ATTTGA" "-L merges the near-identical alt alleles"
is "$(call_site call_cluster.vcf 10 | cut -f1 -d:)" "1/1" "-L turns the 1/2 call into 1/1"
is "$(call_site call_cluster.vcf 10 | cut -f3 -d:)" "0,47" "-L sums AD onto the surviving allele"
is "$(call_site call_cluster.vcf 10 | cut -f8 -d:)" "47" "-L recomputes MAD from the folded AD"
# the whole point of merging after update_vcf_info: the genotyper saw every candidate, so the
# site-level evidence is identical to the unmerged run
is "$(call_site call_cluster.vcf 6)" "$(call_site call_no_cluster.vcf 6)" "-L leaves QUAL unchanged"
is "$(echo "$(call_site call_cluster.vcf 8)" | grep -o 'DP=[0-9]*')" "$(echo "$(call_site call_no_cluster.vcf 8)" | grep -o 'DP=[0-9]*')" "-L leaves DP unchanged"
# merging after flatten_common_allele_ends keeps the surviving allele anchored where it was
is "$(call_site call_cluster.vcf 4)" "$(call_site call_no_cluster.vcf 4)" "-L leaves REF unchanged"
is "$(awk -F'\t' '$3==">1>9"{print $2}' call_cluster.vcf)" "$(awk -F'\t' '$3==">1>9"{print $2}' call_no_cluster.vcf)" "-L leaves POS unchanged"
# assert the VALUES, not just the count: vg emits GL i-major, and a spec-ordered fold would
# still produce three fields while putting the wrong number in the middle
is "$(call_site call_cluster.vcf 10 | cut -f4 -d:)" "-111.563773,-26.511924,-4.458384" "-L folds GL by max over the merged genotype classes"
is "$(grep -v '^#' call_cluster.vcf | grep -c '')" "$(grep -v '^#' call_no_cluster.vcf | grep -c '')" "-L never deletes a variant record"
is "$(call_site call_cluster.vcf 8 | grep -o 'MAT=[^;]*')" "MAT=2>1:0.714" "-L records what it merged in MAT"
is "$(grep -c 'ID=MAT' call_no_cluster.vcf)" "0" "the MAT header is absent when -L is not used"
is "$(grep -c 'ID=MAT' call_cluster.vcf)" "1" "and present when it is"

# the threshold is the same weighted-Jaccard as vg deconstruct -L: these two alts score 5/7
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.7143 --cluster-min-len 0 > call_cluster_hi.vcf 2>/dev/null
is "$(call_site call_cluster_hi.vcf 5)" "ATTTGA,ATTTGG" "-L just above the Jaccard does not merge"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.7142 --cluster-min-len 0 > call_cluster_lo.vcf 2>/dev/null
is "$(call_site call_cluster_lo.vcf 5)" "ATTTGA" "-L just below the Jaccard merges"
# and vg deconstruct flips at the same place on the same graph.  Assert the merge DECISION -- how
# many alts survive -- not which one survives: the two tools pick a cluster's survivor by different
# rules (vg call by allele depth, vg deconstruct by get_traversal_order's ranking), so pinning the
# sequence here would make an unrelated path rename in the fixture fail this test with a message
# blaming the threshold.
n_alts() { awk -F'\t' '$2==9{print split($5,a,",")}' "$1"; }
vg deconstruct small_cluster_call.gfa -p x -L 0.7143 --cluster-min-len 0 2>/dev/null > decon_hi.vcf
vg deconstruct small_cluster_call.gfa -p x -L 0.7142 --cluster-min-len 0 2>/dev/null > decon_lo.vcf
is "$(n_alts decon_hi.vcf)" "$(n_alts call_cluster_hi.vcf)" "vg call and vg deconstruct agree above the threshold"
is "$(n_alts decon_lo.vcf)" "$(n_alts call_cluster_lo.vcf)" "vg call and vg deconstruct agree below the threshold"

# --cluster-min-len gates per site, exactly as in vg deconstruct.  This site's core length is 6 bp.
for N in 0 6 7 50 ; do
    vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 --cluster-min-len $N > call_minlen_$N.vcf 2>/dev/null
done
is "$(call_site call_minlen_0.vcf 5)"  "ATTTGA"        "--cluster-min-len 0 merges everywhere"
is "$(call_site call_minlen_6.vcf 5)"  "ATTTGA"        "--cluster-min-len at the site length still merges"
is "$(call_site call_minlen_7.vcf 5)"  "ATTTGA,ATTTGG" "--cluster-min-len above the site length gates merging"
is "$(call_site call_minlen_50.vcf 5)" "ATTTGA,ATTTGG" "--cluster-min-len 50 restricts merging to SVs"
# 50 is the DEFAULT, so -L on its own must behave like --cluster-min-len 50 and leave this 6bp site
# alone.  Every merge assertion above has to pass --cluster-min-len 0 for exactly this reason.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 > call_minlen_default.vcf 2>/dev/null
diff call_minlen_50.vcf call_minlen_default.vcf
is "$?" 0 "-L alone gates like --cluster-min-len 50, the default"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --cluster-min-len 50 > call_minlen_only.vcf 2>/dev/null
diff call_no_cluster.vcf call_minlen_only.vcf
is "$?" 0 "--cluster-min-len without -L is a no-op"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --cluster-min-len 50 2>&1 >/dev/null | grep -q "no effect without -L"
is "$?" 0 "--cluster-min-len without -L warns"

# The gate must read the alleles the record emits, not the traversal finder's candidate list.
# Adding a long branch with no reads -- absent from the genotype and from every AT entry -- must not
# change whether merging happens.
cp small_cluster_call.gfa minlen_ghost.gfa
printf "S\t9999\t%s\n" "$(printf 'ACGT%.0s' $(seq 1 50))" >> minlen_ghost.gfa
printf "L\t1\t+\t9999\t+\t0M\n" >> minlen_ghost.gfa
printf "L\t9999\t+\t9\t+\t0M\n" >> minlen_ghost.gfa
vg view -Fv minlen_ghost.gfa > minlen_ghost.vg
vg pack -x minlen_ghost.vg -g small_cluster_call.gam -o minlen_ghost.pack
vg call minlen_ghost.vg -k minlen_ghost.pack -p x -L 0.6 --cluster-min-len 50 > call_minlen_ghost.vcf 2>/dev/null
is "$(call_site call_minlen_ghost.vcf 5)" "ATTTGA,ATTTGG" "--cluster-min-len ignores uncalled candidate traversals"

# When two alleles merge, the surviving sequence must be the better-supported one: the finder ranks
# by length-weighted average flow on large snarls, which can put a short lightly-supported allele
# first, and merging into it would emit the minority sequence as a hom-alt carrying the pooled depth.
vg view -Fv nesting/merge_support.gfa > merge_support.vg
vg sim -x merge_support.vg -P hapA -n 40  -l 120 -a -s 5 > merge_a.gam
vg sim -x merge_support.vg -P hapB -n 120 -l 120 -a -s 6 > merge_b.gam
cat merge_a.gam merge_b.gam > merge_support.gam
vg pack -x merge_support.vg -g merge_support.gam -o merge_support.pack
vg call merge_support.vg -k merge_support.pack -p x > merge_plain.vcf 2>/dev/null
vg call merge_support.vg -k merge_support.pack -p x -L 0.4 --cluster-min-len 0 > merge_L.vcf 2>/dev/null
is "$(grep -v '^#' merge_plain.vcf | cut -f10 | cut -f3 -d:)" "0,20,62" "the unmerged call has a clear minority and majority allele"
is "$(grep -v '^#' merge_L.vcf | cut -f5)" "$(grep -v '^#' merge_plain.vcf | cut -f5 | cut -f2 -d,)" "the merge keeps the better-supported allele's sequence"
is "$(grep -v '^#' merge_L.vcf | cut -f10 | cut -f3 -d:)" "0,82" "the merge folds both alleles' depth onto the survivor"

rm -f minlen_ghost.gfa minlen_ghost.vg minlen_ghost.pack call_minlen_ghost.vcf
rm -f merge_support.vg merge_support.gam merge_support.pack merge_a.gam merge_b.gam merge_plain.vcf merge_L.vcf

# -L must be rejected rather than silently ignored: -v and -G/-T never reach
# VCFOutputCaller::emit_variant at all, while -B reaches it but reports QUAL/XADL/lowxadl for a
# het the merge would erase.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 -v small/x.vcf.gz >/dev/null 2>&1
is "$?" 1 "-L is rejected when genotyping a VCF"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 -G >/dev/null 2>&1
is "$?" 1 "-L is rejected with GAF output"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 -B >/dev/null 2>&1
is "$?" 1 "-L is rejected with the ratio caller"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 --legacy >/dev/null 2>&1
is "$?" 1 "-L is rejected with the legacy caller"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 -d 1 2>&1 >/dev/null | grep -q "no effect at ploidy 1"
is "$?" 0 "-L warns at ploidy 1"
# -R sets ploidy per contig and the contig list is not known until the graph is loaded, so a -R rule
# does not trigger the warning -- and must not fail the run either
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 0.6 -R 'nosuchcontig:1' >/dev/null 2>&1
is "$?" 0 "-L accepts a -R rule that matches no called contig"
# unlike vg deconstruct, which clamps, we reject: -L 5 is a plausible typo for -L 0.5
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -L 5 >/dev/null 2>&1
is "$?" 1 "-L out of range is an error"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --cluster-min-len=-1 >/dev/null 2>&1
is "$?" 1 "negative --cluster-min-len is an error"
# --cluster-post shipped through v1.76 as an accepted no-op.  It stays accepted for one release so
# pipelines carrying it do not die on an unrecognized option, but it now says so.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --cluster-post >/dev/null 2>&1
is "$?" 0 "the deprecated --cluster-post option is still accepted"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --cluster-post 2>&1 >/dev/null | grep -q "deprecated and ignored"
is "$?" 0 "--cluster-post warns that it is ignored"

# In a nested run a merged parent deliberately disagrees with its own child records: the parent
# gives the collapsed view of a large variant, the children give the precise one.  MAT records it.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --top-down -L 0.6 --cluster-min-len 0 > call_cluster_td.vcf 2>/dev/null
is "$(call_site call_cluster_td.vcf 5)" "ATTTGA" "-L merges under --top-down"
# -Y is the one nested combination refused: it writes "*" in a child record to mean "an upstream
# deletion covers this site", and the merge can absorb the parent allele that deletion came from,
# leaving the "*" referring to nothing in the file.  Malformed, not merely lossy.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --top-down -Y -L 0.6 >/dev/null 2>&1
is "$?" 1 "-L is rejected with -Y/--star-allele"

# The gref workflow, which is what nested -L is for: the reference bypasses a 121bp insertion, so
# the insertion's interior is only reachable as a gref fragment contig.  The insertion is an SV by
# the --cluster-min-len 50 default, so it merges without the gate having to be turned off.
A20=$(printf 'A%.0s' $(seq 20)); C60=$(printf 'C%.0s' $(seq 60))
A60=$(printf 'A%.0s' $(seq 60)); T20=$(printf 'T%.0s' $(seq 20))
{ printf "H\tVN:Z:1.0\n"
  printf "S\t1\t$A20\nS\t2\t$C60\nS\t3\tG\nS\t4\tT\nS\t5\t$A60\nS\t6\t$T20\n"
  printf "L\t1\t+\t6\t+\t0M\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n"
  printf "L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\n"
  printf "P\tx\t1+,6+\t*,*\n"
  printf "P\ts#1#hapA\t1+,2+,3+,5+,6+\t*,*,*,*,*\n"
  printf "P\ts#2#hapB\t1+,2+,4+,5+,6+\t*,*,*,*,*\n"; } > nlg.gfa
vg paths --compute-gref -Q x --min-gref-len 1 -x nlg.gfa > nlg.pg 2>/dev/null
vg sim -x nlg.pg -P "s#1#hapA#0" -n 200 -l 30 -a -s 3 > nlg_a.gam 2>/dev/null
vg sim -x nlg.pg -P "s#2#hapB#0" -n 200 -l 30 -a -s 4 > nlg_b.gam 2>/dev/null
cat nlg_a.gam nlg_b.gam > nlg.gam
vg pack -x nlg.pg -g nlg.gam -o nlg.pack 2>/dev/null
nlg_field() { awk -F'\t' -v c="$2" -v f="$3" '$1==c && $1!~/^#/{print $f}' "$1"; }
vg call nlg.pg -k nlg.pack -P gref_x -A > nlg_plain.vcf 2>/dev/null
vg call nlg.pg -k nlg.pack -P gref_x -A -L 0.9 > nlg_L.vcf 2>/dev/null
is "$(nlg_field nlg_plain.vcf gref_x 5 | awk '{print split($0,a,",")}')" "2" "the unmerged gref parent has both insertion alleles"
is "$(nlg_field nlg_L.vcf gref_x 5 | awk '{print split($0,a,",")}')" "1" "-L merges them at the default --cluster-min-len, the insertion being an SV"
is "$(nlg_field nlg_L.vcf gref_x 8 | grep -c 'MAT=')" "1" "and says so in MAT"
# the merge is parent-only: the nested record on the gref fragment keeps the precise variant
is "$(nlg_field nlg_L.vcf gref_x_1_alt 1-11)" "$(nlg_field nlg_plain.vcf gref_x_1_alt 1-11)" "the child record on the gref fragment is untouched by the merge"
rm -f nlg.gfa nlg.pg nlg_a.gam nlg_b.gam nlg.gam nlg.pack nlg_plain.vcf nlg_L.vcf

# --bottom-up had zero coverage anywhere in the suite.  It aborted on any nested graph whose
# reference traversal crosses a child snarl, because a child-snarl Visit carries no node.
vg view -Fv nesting/nested_snp_in_del.gfa > bu.vg
vg sim -x bu.vg -n 30 -l 10 -a -s 1 > bu.gam
vg pack -x bu.vg -g bu.gam -o bu.pack
vg call bu.vg -k bu.pack -p x --bottom-up > bu.vcf 2>/dev/null
is "$?" 0 "--bottom-up does not abort on a nested graph"
is "$(grep -vc '^#' bu.vcf)" "1" "--bottom-up emits a record"
# NestedFlowCaller's Snarl-carrying Visits break the GAF emitters: -T aborted in to_mapping and -G
# emitted a header with no records.  Rejected rather than left silently inert.
vg call bu.vg -k bu.pack -p x --bottom-up -T >/dev/null 2>&1
is "$?" 1 "--bottom-up is rejected with -T"
vg call bu.vg -k bu.pack -p x --bottom-up -G >/dev/null 2>&1
is "$?" 1 "--bottom-up is rejected with -G"
rm -f bu.vg bu.gam bu.pack bu.vcf

# --cluster-min-len gates on CORE LENGTH -- the longest allele once the prefix and suffix shared by
# every allele are stripped -- not on the raw snarl interior.  Same rule, same helper, as
# vg deconstruct: without it a 1bp SNP gets gated on the size of whatever snarl contains it.
core_pack() {
    vg view -Fv nesting/$1.gfa > core_$1.vg
    vg sim -x core_$1.vg -P hapA -n 40 -l 60 -a -s 5 > core_a.gam
    vg sim -x core_$1.vg -P hapB -n 40 -l 60 -a -s 6 > core_b.gam
    cat core_a.gam core_b.gam > core_$1.gam
    vg pack -x core_$1.vg -g core_$1.gam -o core_$1.pack
}
core_nalt() { vg call core_$1.vg -k core_$1.pack -p x -L 0.6 --cluster-min-len $2 2>/dev/null | awk -F'\t' '$1!~/^#/{printf "%d", split($5,a,",")}'; }
for g in core_snp_in_flanks core_sv60 core_del61 core_ins49 core_ins50 ; do core_pack $g ; done

is "$(core_nalt core_snp_in_flanks 50)" "2" "a 1bp SNP in a large snarl is not merged at --cluster-min-len 50"
is "$(core_nalt core_snp_in_flanks 1)"  "1" "the same SNP is merged at --cluster-min-len 1"
is "$(core_nalt core_sv60 50)"          "1" "a 60bp SV is still merged at --cluster-min-len 50"
is "$(core_nalt core_del61 50)"         "1" "a 61bp deletion is still merged (the reference allele is measured)"
is "$(core_nalt core_ins49 50)"         "2" "a 49bp insertion is not merged at 50 (the anchor base is not counted)"
is "$(core_nalt core_ins50 50)"         "1" "a 50bp insertion is merged at 50"

# Pure deletions: see the matching block in 26_deconstruct.t.  Both tools must make the SAME merge
# decision, so the ladders are compared directly.  Compare COUNTS, never sequences -- vg call keeps
# the highest-depth allele and vg deconstruct keeps whichever get_traversal_order ranks first, so the
# surviving sequence legitimately differs.
for g in del59_vs_del60 del60_vs_del59ins1 del60_vs_snp inv60_vs_del60 inv60_in_2kb unrelated10_in_2kb ; do core_pack $g ; done
cladder() { for L in 0.99 0.983334 0.983333 0.6 0.1 ; do
    vg call core_$1.vg -k core_$1.pack -p x -L $L 2>/dev/null |
      awk -F'\t' '$1!~/^#/{n=split($5,a,","); if($5=="."||$5=="")n=0; printf "%d",n}' ; done ; }
dladder2() { for L in 0.99 0.983334 0.983333 0.6 0.1 ; do
    vg deconstruct nesting/$1.gfa -p x -L $L 2>/dev/null |
      awk -F'\t' '$1!~/^#/{n=split($5,a,","); if($5=="."||$5=="")n=0; printf "%d",n}' ; done ; }

is "$(cladder del59_vs_del60)"     "22111" "a 59bp and a 60bp deletion merge, flipping at 59/60"
is "$(cladder del60_vs_del59ins1)" "22111" "a deletion merges with a deletion carrying a novel base"
is "$(cladder del60_vs_snp)"       "22222" "a 60bp deletion never merges with a 60bp substitution"
is "$(cladder inv60_vs_del60)"     "22222" "a 60bp inversion never merges with a 60bp deletion"
is "$(cladder inv60_in_2kb)"       "22222" "a 2kb snarl does not make an inversion mergeable"
is "$(cladder unrelated10_in_2kb)" "22222" "a 2kb snarl does not make two unrelated 10bp alleles mergeable"
is "$(cladder del59_vs_del60)"     "$(dladder2 del59_vs_del60)"     "vg call and vg deconstruct agree on deletions"
is "$(cladder inv60_in_2kb)"       "$(dladder2 inv60_in_2kb)"       "vg call and vg deconstruct agree on non-merges"

rm -f core_*.vg core_*.gam core_*.pack

# Crash regressions.  None of these modes had any coverage, which is how all three survived.
# -R/--ploidy-regex bypassed the {1,2} ploidy check that -d gets, and reached the caller as an
# unsupported ploidy.  The check runs where the rule is applied, so a rule that matches none of the
# called contigs is still accepted.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -R 'x:3' >/dev/null 2>&1
is "$?" 1 "-R rejects a ploidy the callers do not implement"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -R 'x:2' >/dev/null 2>&1
is "$?" 0 "-R still accepts ploidy 2"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -R 'chrY:3' >/dev/null 2>&1
is "$?" 0 "-R accepts an unsupported ploidy on a contig that is not being called"

# --ploidy-bed: ploidy per region, which -d and -R cannot express. The BED is validated up front
# rather than left to fail inside a caller, where the message says nothing about the file.
printf 'x\t0\t100\t3\n' > bad_ploidy.bed
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed bad_ploidy.bed >/dev/null 2>&1
is "$?" 1 "--ploidy-bed rejects a ploidy the callers do not implement"
printf 'x\t0\t100\t2\nx\t50\t150\t1\n' > overlap_ploidy.bed
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed overlap_ploidy.bed >/dev/null 2>&1
is "$?" 1 "--ploidy-bed rejects overlapping intervals rather than picking one"
printf 'x\t100\t50\t2\n' > reversed_ploidy.bed
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed reversed_ploidy.bed >/dev/null 2>&1
is "$?" 1 "--ploidy-bed rejects a reversed interval"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed nosuchfile.bed >/dev/null 2>&1
is "$?" 1 "--ploidy-bed rejects a missing file"

# A BED covering nothing on the called contig must leave the run exactly as it was: this is what
# makes the option safe to pass unconditionally in a pipeline.
printf 'chrNotHere\t0\t100\t1\n' > elsewhere_ploidy.bed
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed elsewhere_ploidy.bed 2>/dev/null > ploidy_bed_elsewhere.vcf
vg call small_cluster_call.vg -k small_cluster_call.pack -p x 2>/dev/null > ploidy_bed_none.vcf
is "$(grep -v '^#' ploidy_bed_elsewhere.vcf | md5sum | cut -f1 -d' ')" "$(grep -v '^#' ploidy_bed_none.vcf | md5sum | cut -f1 -d' ')" "a --ploidy-bed matching no called contig changes nothing"


rm -f bad_ploidy.bed overlap_ploidy.bed reversed_ploidy.bed elsewhere_ploidy.bed
rm -f ploidy_bed_elsewhere.vcf ploidy_bed_none.vcf
# -G built a re-indexed traversal vector but kept the original reference index, then indexed the
# small vector with it.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -G > call_gaf.gaf 2>/dev/null
is "$?" 0 "-G does not crash"
is "$(grep -c . call_gaf.gaf)" "11" "-G emits a GAF record per called traversal"
rm -f call_gaf.gaf

rm -f small_cluster_call.gfa small_cluster_call.vg small_cluster_call.gam small_cluster_call.pack call_a.gam call_b.gam
rm -f call_no_cluster.vcf call_cluster.vcf call_cluster_off.vcf call_cluster_hi.vcf call_cluster_lo.vcf call_cluster_td.vcf
rm -f call_minlen_0.vcf call_minlen_6.vcf call_minlen_7.vcf call_minlen_50.vcf call_minlen_only.vcf call_minlen_default.vcf decon_hi.vcf decon_lo.vcf

# Test: Nested calling with vg call --top-down (basic test)
# Use a larger graph with clear variants
vg construct -r small/x.fa -v small/x.vcf.gz -a > nested_call_test.vg
vg sim -x nested_call_test.vg -n 500 -l 50 -a -s 42 > nested_call_test.gam
vg pack -x nested_call_test.vg -g nested_call_test.gam -o nested_call_test.pack
vg call nested_call_test.vg -k nested_call_test.pack -p x --top-down > nested_call_test.vcf 2>/dev/null
# Should produce at least one variant (the small/x graph has multiple variants)
NESTED_VARIANT_COUNT=$(grep -v "^#" nested_call_test.vcf | wc -l | tr -d ' ')
is $(if [ "$NESTED_VARIANT_COUNT" -ge 1 ]; then echo "1"; else echo "0"; fi) "1" "nested vg call produces at least one variant"

rm -f nested_call_test.vg nested_call_test.gam nested_call_test.pack nested_call_test.vcf

# Test: Nested calling emits both top-level and child snarls
# This tests the fix for nested snarl VCF emission where genotypes propagate to children
vg view -Fv nesting/nested_snp_in_del.gfa > nested_snp.vg
vg sim -x nested_snp.vg -n 100 -l 5 -a -s 42 > nested_snp.gam
vg pack -x nested_snp.vg -g nested_snp.gam -o nested_snp.pack
vg call nested_snp.vg -k nested_snp.pack --top-down -p x 2>/dev/null > nested_snp.vcf
# Should have exactly 2 variant lines: one for top-level snarl (1->6) and one for nested snarl (2->5)
NESTED_LINE_COUNT=$(grep -v "^#" nested_snp.vcf | wc -l | tr -d ' ')
is "$NESTED_LINE_COUNT" "2" "nested vg call emits both top-level and child snarl variants"

rm -f nested_snp.vg nested_snp.gam nested_snp.pack nested_snp.vcf

## Read-level genotyping over a nested site with a star allele
# nested_snp_in_del.gfa: x = 1>2>3>5>6, y0 = 1>2>4>5>6, y1 = 1>6 (deletes the
# region holding the nested SNP). Reads from x and y1 only, so the top-level site
# is a het deletion and the nested site is traversed by one haplotype.
vg view -Fv nesting/nested_snp_in_del.gfa > ns.vg 2>/dev/null
vg sim -x ns.vg -P x -n 150 -l 4 -a -s 7 > ns_het.gam 2>/dev/null
vg sim -x ns.vg -P 'a#2#y1#0' -n 150 -l 4 -a -s 8 >> ns_het.gam 2>/dev/null
vg pack -x ns.vg -g ns_het.gam -o ns_het.pack 2>/dev/null
vg call ns.vg -k ns_het.pack --top-down -Y -p x --read-likelihood --gam ns_het.gam 2>/dev/null > ns_rl.vcf

is $(grep -v "^#" ns_rl.vcf | wc -l | tr -d ' ') "2" "nested read-likelihood emits both the parent and the child site"

# Regression: reads traversing the deletion edge touch only the snarl's boundary
# nodes. A node-based informativeness test discarded them, which left the parent
# with reference-supporting reads only, called it hom-ref, and dropped the record
# entirely. The deletion must be called, and its DP must include those reads.
is $(grep -v "^#" ns_rl.vcf | awk '$4=="CATG" && $5=="C"' | wc -l | tr -d ' ') "1" "the parent deletion is called from boundary-to-boundary reads"
is $(grep -v "^#" ns_rl.vcf | awk -F'\t' '$4=="CATG"{nk=split($9,k,":"); di=0; for(j=1;j<=nk;j++) if(k[j]=="DP") di=j; split($10,f,":"); print f[di]}') "300" "deletion-spanning reads are counted, not discarded"

# Regression: half these reads are reverse strand. Failing to flip them meant they
# anchored on nothing and scored against the wrong allele, which turned the nested
# site into a spurious het SNP instead of a star allele.
is $(grep -v "^#" ns_rl.vcf | awk '$5=="*"' | wc -l | tr -d ' ') "1" "the nested site is a star allele, not a strand-artefact het SNP"

# Independent nested calling (-A): every snarl genotyped on its own reads, with no
# parent restriction and no phase propagation.
vg call ns.vg -k ns_het.pack -A -p x --read-likelihood --gam ns_het.gam 2>/dev/null > ns_rl_a.vcf
is "$?" "0" "--read-likelihood works with -A independent nested calling"
is $(grep -v "^#" ns_rl_a.vcf | awk '$4=="CATG" && $5=="C"' | wc -l | tr -d ' ') "1" "-A independent calling finds the deletion"
is $(grep -c "PS=" ns_rl_a.vcf | tr -d ' ') "0" "-A emits no PS tags, since it does not propagate phase"

# Effective ploidy at a nested site only one parent haplotype traverses.
# Reads: 100 from x (ref SNP), 30 from y0 (alt SNP), 100 from y1 (deletion). The
# child site is reached by one haplotype only, so it must be genotyped haploid.
# Genotyping it as diploid lets a spurious heterozygote absorb the 30 minority
# reads for free, which shows up as a badly depressed GQ (38 rather than 256).
vg sim -x ns.vg -P x -n 100 -l 4 -a -s 11 > ns_mix.gam 2>/dev/null
vg sim -x ns.vg -P 'a#1#y0#0' -n 30 -l 4 -a -s 12 >> ns_mix.gam 2>/dev/null
vg sim -x ns.vg -P 'a#2#y1#0' -n 100 -l 4 -a -s 13 >> ns_mix.gam 2>/dev/null
vg pack -x ns.vg -g ns_mix.gam -o ns_mix.pack 2>/dev/null
# --mismap-max is pinned rather than defaulted, because every read here comes
# straight from `vg sim` and so carries MAPQ 0. That makes e_r equal to the cap for
# all 230 reads, and the GQ assertion below becomes a measurement of the cap rather
# than of the ploidy handling it is named for: the same correct 1/0 call scores GQ
# 161, 83 or 24 at caps 0.5, 0.7 and 0.9. Pinning keeps this a ploidy test. The cap
# itself is swept against real data, where reads have real mapping qualities.
vg call ns.vg -k ns_mix.pack --top-down -Y -p x --read-likelihood --gam ns_mix.gam --mismap-max 0.5 2>/dev/null > ns_mix.vcf

is $(grep -v "^#" ns_mix.vcf | awk '$5=="*"' | wc -l | tr -d ' ') "1" "mixed-read nested site still yields a star allele"
STAR_GQ=$(grep -v "^#" ns_mix.vcf | awk -F'\t' '$5=="*"{nk=split($9,k,":"); qi=0; for(j=1;j<=nk;j++) if(k[j]=="GQ") qi=j; split($10,f,":"); print f[qi]}')
is $(if [ "${STAR_GQ}" -gt 100 ]; then echo 1; else echo 0; fi) "1" "a singly-traversed nested site is genotyped at its own ploidy, not diluted by a spurious het"

rm -f ns.vg ns_het.gam ns_het.pack ns_rl.vcf ns_rl_a.vcf ns_mix.gam ns_mix.pack ns_mix.vcf

# Test: Star allele option validation (-Y requires --top-down)
vg construct -r small/x.fa -v small/x.vcf.gz > star_test.vg
vg sim -x star_test.vg -n 100 -l 20 -a -s 1 > star_test.gam
vg pack -x star_test.vg -g star_test.gam -o star_test.pack
vg call star_test.vg -k star_test.pack -Y 2>&1 | grep -q "requires --top-down"
is "$?" 0 "star allele option requires top-down mode"

rm -f star_test.vg star_test.gam star_test.pack

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
ND_00_NONREF=$(grep -v "^#" nd_00.vcf | grep -v "0/0" | wc -l | tr -d ' ')
is "$ND_00_NONREF" "0" "nested_snp_in_del 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/SNP - reads from x and y0 (both traverse nested snarl)
vg sim -x nesting/nested_snp_in_del.gfa -P x -n 50 -l 2 -a -s 10 > nd_01.gam
vg sim -x nesting/nested_snp_in_del.gfa -m a -n 50 -l 2 -a -s 11 >> nd_01.gam
vg pack -x nesting/nested_snp_in_del.gfa -g nd_01.gam -o nd_01.pack
vg call nesting/nested_snp_in_del.gfa -k nd_01.pack --top-down -p x 2>/dev/null > nd_01.vcf
# Should have 2 variants (top-level and nested), both het
ND_01_COUNT=$(grep -v "^#" nd_01.vcf | wc -l | tr -d ' ')
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
ND_12_COUNT=$(grep -v "^#" nd_12.vcf | wc -l | tr -d ' ')
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
NO_MISSING_GT=$(grep ">2>5" star.vcf | cut -f10 | cut -d: -f1 | grep -v "\." | wc -l | tr -d ' ')
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

# Test: triple nested with gref paths should all have quality
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/triple_nested.gfa > tnq_ap.gfa 2>/dev/null
vg sim -x tnq_ap.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 51 > tnq.gam
vg sim -x tnq_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 52 >> tnq.gam
vg pack -x tnq_ap.gfa -g tnq.gam -o tnq.pack
vg call tnq_ap.gfa -k tnq.pack --top-down -P gref_x 2>/dev/null > tnq.vcf

# All variant lines should have non-zero QUAL
TNQ_ZERO_QUAL=$(grep -v "^#" tnq.vcf | awk -F'\t' '$6 == "0" || $6 == "."' | wc -l | tr -d ' ')
is "$TNQ_ZERO_QUAL" "0" "triple nested calls all have non-zero QUAL"

# All variant lines should have GQ in FORMAT
TNQ_ALL_GQ=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "GQ" | wc -l | tr -d ' ')
is "$TNQ_ALL_GQ" "0" "triple nested calls all have GQ field"

# All variant lines should have GL (Genotype Likelihood) in FORMAT
TNQ_ALL_GL=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "GL" | wc -l | tr -d ' ')
is "$TNQ_ALL_GL" "0" "triple nested calls all have GL field"

# All variant lines should have GP (Genotype Posterior) in FORMAT
TNQ_ALL_GP=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "GP" | wc -l | tr -d ' ')
is "$TNQ_ALL_GP" "0" "triple nested calls all have GP field"

# All variant lines should have XD (Expected Depth) in FORMAT
TNQ_ALL_XD=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "XD" | wc -l | tr -d ' ')
is "$TNQ_ALL_XD" "0" "triple nested calls all have XD field"

# All variant lines should have AD (Allelic Depth) in FORMAT
TNQ_ALL_AD=$(grep -v "^#" tnq.vcf | cut -f9 | grep -v "AD" | wc -l | tr -d ' ')
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
}' | wc -l | tr -d ' ')
is "$TNQ_INVALID_GQ" "0" "triple nested calls have valid GQ values (0-256)"

# All variant lines should have LV (nesting level) in INFO
TNQ_ALL_LV=$(grep -v "^#" tnq.vcf | cut -f8 | grep -v "LV=" | wc -l | tr -d ' ')
is "$TNQ_ALL_LV" "0" "triple nested calls all have LV tag"

# Nested variants (LV > 0) should have PS (parent snarl) in INFO
TNQ_NESTED_NO_PS=$(grep -v "^#" tnq.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/ && $8 !~ /PS=/' | wc -l | tr -d ' ')
is "$TNQ_NESTED_NO_PS" "0" "nested calls (LV>0) all have PS tag"

# A record with no ancestor anywhere is LV=0 with no contig hop, and has no PS tag.
TNQ_TOPLEVEL_HAS_PS=$(grep -v "^#" tnq.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/ && $8 ~ /PS=/' | wc -l | tr -d ' ')
is "$TNQ_TOPLEVEL_HAS_PS" "0" "calls with no ancestor at all (LV=0, CH=0) do not have PS tag"

# gref_x_2_alt is in the cover but carries no call, so vg call does not declare it either.
is $(grep -c "^##contig" tnq.vcf) 2 "vg call does not declare reference contigs with no records"
is $(grep -v "^#" tnq.vcf | cut -f1 | sort -u | wc -l | tr -d ' ') 2 "every contig with a call is declared"

# But a record that is top-level on its OWN contig (LV=0) while being nested in the snarl
# tree (CH>0) MUST keep PS.  It is the only in-VCF link back to the enclosing base-contig
# site, and vcfbub's rescue of the children of popped bubbles is keyed on it.
TNQ_PERCONTIG_TOP_HAS_PS=$(grep -v "^#" tnq.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 !~ /CH=0/ && $8 !~ /PS=/' | wc -l | tr -d ' ')
is "$TNQ_PERCONTIG_TOP_HAS_PS" "0" "per-contig top-level calls nested in the tree still carry PS"

rm -f tnq_ap.gfa tnq.gam tnq.pack tnq.vcf

# =============================================================================
# nested_snp_in_ins.gfa tests
# Graph structure:
#   x (ref):  1->6                    (short ref, allele 0 at top-level)
#   y0 (a#1): 1->2->4->5->6           (insertion with SNP A, allele 1 at top, allele 1 at nested)
#   y1 (a#2): 1->2->3->5->6           (insertion with SNP T, allele 2 at top, allele 0 at nested)
# Top-level snarl: 1->6, Nested snarl: 2->5
# NOTE: ref path x bypasses the nested snarl, so we need gref covers to call nested variants
# =============================================================================

# Compute gref cover first (creates gref_x plus gref_x_1_alt, gref_x_2_alt, ... covering nodes not on x)
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/nested_snp_in_ins.gfa > ni_ap.gfa

# Test 0/0: homozygous reference - reads only from x path (short ref)
vg sim -x ni_ap.gfa -P x -n 100 -l 2 -a -s 70 > ni_00.gam
vg pack -x ni_ap.gfa -g ni_00.gam -o ni_00.pack
vg call ni_ap.gfa -k ni_00.pack --top-down -P gref_x 2>/dev/null > ni_00.vcf
# 0/0 should produce no non-ref variants (or only ref calls)
NI_00_NONREF=$(grep -v "^#" ni_00.vcf | grep -v "0/0" | wc -l | tr -d ' ')
is "$NI_00_NONREF" "0" "nested_snp_in_ins 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/insertion - reads from x and y1 (a#2 haplotype)
# Note: y1 path is 5bp vs x path 2bp, so need ~2x more y1 reads for balanced boundary coverage
vg sim -x ni_ap.gfa -P x -n 50 -l 2 -a -s 71 > ni_01.gam
vg sim -x ni_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 72 >> ni_01.gam
vg pack -x ni_ap.gfa -g ni_01.gam -o ni_01.pack
vg call ni_ap.gfa -k ni_01.pack --top-down -P gref_x 2>/dev/null > ni_01.vcf
# With gref paths: both top-level and nested variants emitted
NI_01_COUNT=$(grep -v "^#" ni_01.vcf | wc -l | tr -d ' ')
is "$NI_01_COUNT" "2" "nested_snp_in_ins 0/1: het ref/ins produces top-level and nested variants with gref paths"

# Test 1/1: homozygous insertion - reads only from y1 path
vg sim -x ni_ap.gfa -P "a#2#y1#0" -n 100 -l 2 -a -s 73 > ni_11.gam
vg pack -x ni_ap.gfa -g ni_11.gam -o ni_11.pack
vg call ni_ap.gfa -k ni_11.pack --top-down -P gref_x 2>/dev/null > ni_11.vcf
# With gref paths: both top-level and nested variants emitted
NI_11_COUNT=$(grep -v "^#" ni_11.vcf | wc -l | tr -d ' ')
is "$NI_11_COUNT" "2" "nested_snp_in_ins 1/1: homozygous ins produces top-level and nested variants with gref paths"

# Test 1/2: het between two insertion alleles - reads from both y0 and y1
vg sim -x ni_ap.gfa -m a -n 200 -l 2 -a -s 74 > ni_12.gam
vg pack -x ni_ap.gfa -g ni_12.gam -o ni_12.pack
vg call ni_ap.gfa -k ni_12.pack --top-down -P gref_x 2>/dev/null > ni_12.vcf
# With gref paths: both top-level and nested variants emitted
NI_12_COUNT=$(grep -v "^#" ni_12.vcf | wc -l | tr -d ' ')
is "$NI_12_COUNT" "2" "nested_snp_in_ins 1/2: het ins/ins produces top-level and nested variants with gref paths"

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
# NOTE: ref path x bypasses all nesting, so we need gref covers to call nested variants
# =============================================================================

# Compute gref cover first (creates gref_x plus gref_x_1_alt, gref_x_2_alt, ... covering nested nodes)
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/triple_nested.gfa > tn_ap.gfa

# Test 0/0: homozygous reference - reads only from x path
vg sim -x tn_ap.gfa -P x -n 100 -l 2 -a -s 80 > tn_00.gam
vg pack -x tn_ap.gfa -g tn_00.gam -o tn_00.pack
vg call tn_ap.gfa -k tn_00.pack --top-down -P gref_x 2>/dev/null > tn_00.vcf
TN_00_NONREF=$(grep -v "^#" tn_00.vcf | grep -v "0/0" | wc -l | tr -d ' ')
is "$TN_00_NONREF" "0" "triple_nested 0/0: homozygous ref produces no non-ref variants"

# Test 0/1: het ref/insertion - reads from x and y0
# Note: y0 path is 11bp vs x path 2bp, so need ~5x more y0 reads for balanced boundary coverage
vg sim -x tn_ap.gfa -P x -n 50 -l 2 -a -s 81 > tn_01.gam
vg sim -x tn_ap.gfa -P "a#1#y0#0" -n 500 -l 2 -a -s 82 >> tn_01.gam
vg pack -x tn_ap.gfa -g tn_01.gam -o tn_01.pack
vg call tn_ap.gfa -k tn_01.pack --top-down -P gref_x 2>/dev/null > tn_01.vcf
# With gref paths: all 5 nesting levels can be emitted
TN_01_COUNT=$(grep -v "^#" tn_01.vcf | wc -l | tr -d ' ')
is "$TN_01_COUNT" "5" "triple_nested 0/1: het ref/ins produces all 5 nesting level variants with gref paths"

# Test 1/1: homozygous insertion - reads only from y0 path
vg sim -x tn_ap.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 83 > tn_11.gam
vg pack -x tn_ap.gfa -g tn_11.gam -o tn_11.pack
vg call tn_ap.gfa -k tn_11.pack --top-down -P gref_x 2>/dev/null > tn_11.vcf
# With gref paths: 1 variant emitted (top-level insertion only)
# All nested snarls are 0/0 because y0 matches gref_x_1_alt reference at all levels
# (y0 goes through 313, and gref_x_1_alt also goes through 313)
TN_11_COUNT=$(grep -v "^#" tn_11.vcf | wc -l | tr -d ' ')
is "$TN_11_COUNT" "1" "triple_nested 1/1: homozygous ins produces only top-level variant"

# Test 1/2: het between insertion alleles - reads from y0 and y1 (differ at deepest SNP)
vg sim -x tn_ap.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 84 > tn_12.gam
vg sim -x tn_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 85 >> tn_12.gam
vg pack -x tn_ap.gfa -g tn_12.gam -o tn_12.pack
vg call tn_ap.gfa -k tn_12.pack --top-down -P gref_x 2>/dev/null > tn_12.vcf
# With gref paths: all 5 nesting levels can be emitted
TN_12_COUNT=$(grep -v "^#" tn_12.vcf | wc -l | tr -d ' ')
is "$TN_12_COUNT" "5" "triple_nested 1/2: het ins/ins produces all 5 nesting level variants with gref paths"

rm -f tn_ap.gfa tn_00.gam tn_00.pack tn_00.vcf tn_01.gam tn_01.pack tn_01.vcf
rm -f tn_11.gam tn_11.pack tn_11.vcf tn_12.gam tn_12.pack tn_12.vcf

# =============================================================================
# Multi-level SNP test (triple_nested_multisnp.gfa)
# This graph has SNPs at multiple nesting levels:
# - y0 matches gref_x_1_alt at all levels (will be chosen as gref reference)
# - y1 differs from gref_x_1_alt at all nested levels
# When simulating from y1, we should get variants at all 4 nesting levels
# =============================================================================

vg paths -x nesting/triple_nested_multisnp.gfa -Q x --compute-gref --min-gref-len 1 > tn_ms_ap.gfa
vg sim -x tn_ms_ap.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 100 > tn_ms.gam
vg pack -x tn_ms_ap.gfa -g tn_ms.gam -o tn_ms.pack
vg call tn_ms_ap.gfa -k tn_ms.pack --top-down -P gref_x 2>/dev/null > tn_ms.vcf
# Should get 4 variants: top-level + 3 nested SNPs (all at 1/1)
TN_MS_COUNT=$(grep -v "^#" tn_ms.vcf | wc -l | tr -d ' ')
is "$TN_MS_COUNT" "4" "triple_nested_multisnp 1/1: homozygous alt produces variants at all 4 nesting levels"
# Verify all variants are 1/1 (homozygous alt)
TN_MS_HOM=$(grep -v "^#" tn_ms.vcf | cut -f10 | cut -d: -f1 | grep -c "1/1")
is "$TN_MS_HOM" "4" "triple_nested_multisnp 1/1: all 4 variants are homozygous alt"
# Verify LV tags span levels 0-3
TN_MS_TOP=$(grep -v "^#" tn_ms.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/' | wc -l | tr -d ' ')
TN_MS_NESTED=$(grep -v "^#" tn_ms.vcf | awk -F'\t' '$8 !~ /LV=0/ || $8 !~ /CH=0/' | wc -l | tr -d ' ')
is "$TN_MS_TOP" "1" "triple_nested_multisnp: one variant with no ancestor at all"
is "$TN_MS_NESTED" "3" "triple_nested_multisnp: three nested variants"

rm -f tn_ms_ap.gfa tn_ms.gam tn_ms.pack tn_ms.vcf

# =============================================================================
# Short reference bypass tests
# NOTE: ref path x bypasses nested structures, so we need gref covers to call nested variants
# =============================================================================

# Test nested_snp_in_nested_ins.gfa - ref bypasses all nested structures
# Compute gref cover first (creates gref_x plus gref_x_1_alt, ... covering nested nodes)
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/nested_snp_in_nested_ins.gfa > bypass_ap.gfa
vg sim -x bypass_ap.gfa -m a -n 100 -l 2 -a -s 60 > bypass.gam
vg pack -x bypass_ap.gfa -g bypass.gam -o bypass.pack
vg call bypass_ap.gfa -k bypass.pack --top-down -P gref_x 2>/dev/null > bypass.vcf
BYPASS_EXIT=$?
is "$BYPASS_EXIT" "0" "nested_snp_in_nested_ins: vg call handles short-ref nested graph with gref paths without crashing"

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
NA_DEL_COUNT=$(grep -v "^#" na_del.vcf | wc -l | tr -d ' ')
is "$NA_DEL_COUNT" "2" "--top-down -a: nested_snp_in_del 0/0 emits both snarls"

# Verify top-level is 0/0 (use awk to match ID column exactly)
NA_DEL_TOP_GT=$(awk -F'\t' '$3 == ">1>6" {print $10}' na_del.vcf | cut -d: -f1)
is "$NA_DEL_TOP_GT" "0/0" "--top-down -a: nested_snp_in_del top-level is 0/0"

# Verify nested is 0/0 (use awk to match ID column exactly)
NA_DEL_NEST_GT=$(awk -F'\t' '$3 == ">2>5" {print $10}' na_del.vcf | cut -d: -f1)
is "$NA_DEL_NEST_GT" "0/0" "--top-down -a: nested_snp_in_del nested is 0/0"

rm -f na_del.gam na_del.pack na_del.vcf

# Test 2: nested_snp_in_ins 0/0 with --top-down -a (with gref paths)
# When ref bypasses nested snarl and both alleles are ref, nested NOT emitted
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/nested_snp_in_ins.gfa > na_ins_ap.gfa 2>/dev/null
vg sim -x na_ins_ap.gfa -P x -n 100 -l 2 -a -s 200 > na_ins_00.gam
vg pack -x na_ins_ap.gfa -g na_ins_00.gam -o na_ins_00.pack
vg call na_ins_ap.gfa -k na_ins_00.pack --top-down -a -P gref_x 2>/dev/null > na_ins_00.vcf

# Count variant lines (should be 1: only top-level, nested not emitted)
NA_INS_00_COUNT=$(grep -v "^#" na_ins_00.vcf | wc -l | tr -d ' ')
is "$NA_INS_00_COUNT" "1" "--top-down -a: nested_snp_in_ins 0/0 emits only top-level (ref spans nested)"

rm -f na_ins_00.gam na_ins_00.pack na_ins_00.vcf

# Test 3: nested_snp_in_ins 0/1 with --top-down -a (with gref paths)
# When ref bypasses nested but alt traverses, nested should have ./X genotype
vg sim -x na_ins_ap.gfa -P x -n 50 -l 2 -a -s 200 > na_ins_01.gam
vg sim -x na_ins_ap.gfa -P "a#1#y0#0" -n 100 -l 2 -a -s 201 >> na_ins_01.gam
vg pack -x na_ins_ap.gfa -g na_ins_01.gam -o na_ins_01.pack
vg call na_ins_ap.gfa -k na_ins_01.pack --top-down -a -P gref_x 2>/dev/null > na_ins_01.vcf

# Count variant lines (should be 2: top-level + nested)
NA_INS_01_COUNT=$(grep -v "^#" na_ins_01.vcf | wc -l | tr -d ' ')
is "$NA_INS_01_COUNT" "2" "--top-down -a: nested_snp_in_ins 0/1 emits both snarls"

# Verify nested has missing allele marker (.)
NA_INS_01_NEST_GT=$(grep ">2>5" na_ins_01.vcf | cut -f10 | cut -d: -f1)
NA_INS_01_HAS_MISSING=$(echo "$NA_INS_01_NEST_GT" | grep -c "\.")
is "$NA_INS_01_HAS_MISSING" "1" "--top-down -a: nested_snp_in_ins 0/1 nested has missing allele (.)"

rm -f na_ins_ap.gfa na_ins_01.gam na_ins_01.pack na_ins_01.vcf

# Test 4: triple_nested 0/0 with --top-down -a (with gref paths)
# When ref bypasses all nested snarls and genotype is 0/0, only top-level emitted
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/triple_nested.gfa > na_tn_ap.gfa 2>/dev/null
vg sim -x na_tn_ap.gfa -P x -n 100 -l 2 -a -s 210 > na_tn_00.gam
vg pack -x na_tn_ap.gfa -g na_tn_00.gam -o na_tn_00.pack
vg call na_tn_ap.gfa -k na_tn_00.pack --top-down -a -P gref_x 2>/dev/null > na_tn_00.vcf

# Count variant lines (should be 1: only top-level, nested not emitted since ref spans them)
NA_TN_00_COUNT=$(grep -v "^#" na_tn_00.vcf | wc -l | tr -d ' ')
is "$NA_TN_00_COUNT" "1" "--top-down -a: triple_nested 0/0 emits only top-level (ref spans all nested)"

# Verify top-level is 0/0
NA_TN_00_GT=$(grep -v "^#" na_tn_00.vcf | cut -f10 | cut -d: -f1)
is "$NA_TN_00_GT" "0/0" "--top-down -a: triple_nested 0/0 top-level is 0/0"

rm -f na_tn_00.gam na_tn_00.pack na_tn_00.vcf

# Test 5: triple_nested 0/1 with --top-down -a (with gref paths)
# When ref bypasses nested but alt traverses, nested snarls should have ./X genotype
vg sim -x na_tn_ap.gfa -P x -n 50 -l 2 -a -s 211 > na_tn_01.gam
vg sim -x na_tn_ap.gfa -P "a#1#y0#0" -n 500 -l 2 -a -s 212 >> na_tn_01.gam
vg pack -x na_tn_ap.gfa -g na_tn_01.gam -o na_tn_01.pack
vg call na_tn_ap.gfa -k na_tn_01.pack --top-down -a -P gref_x 2>/dev/null > na_tn_01.vcf

# Count variant lines (should be 5: all nesting levels emitted with -a)
NA_TN_01_COUNT=$(grep -v "^#" na_tn_01.vcf | wc -l | tr -d ' ')
is "$NA_TN_01_COUNT" "5" "--top-down -a: triple_nested 0/1 emits all 5 nesting levels"

# Nested snarls (LV > 0) should have missing allele (.) for the spanning ref
NA_TN_01_NESTED_MISSING=$(grep -v "^#" na_tn_01.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/ {print $10}' | cut -d: -f1 | grep -c "\.")
NA_TN_01_NESTED_COUNT=$(grep -v "^#" na_tn_01.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | wc -l | tr -d ' ')
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
AS_VARIANT_COUNT=$(grep -v "^#" all_snarls_test.vcf | wc -l | tr -d ' ')
AS_LV_COUNT=$(grep -v "^#" all_snarls_test.vcf | grep -c "LV=")
is "$AS_LV_COUNT" "$AS_VARIANT_COUNT" "-A flag: all variants have LV tag"

rm -f all_snarls_test.vg all_snarls_test.gam all_snarls_test.pack all_snarls_test.vcf

# Test: -A flag on nested graph produces variants at all nesting levels
vg view -Fv nesting/nested_snp_in_del.gfa > as_nested.vg
vg sim -x as_nested.vg -m a -n 100 -l 2 -a -s 301 > as_nested.gam
vg pack -x as_nested.vg -g as_nested.gam -o as_nested.pack
vg call as_nested.vg -k as_nested.pack -A -p x > as_nested.vcf 2>/dev/null

# Should produce variants at both nesting levels
AS_NESTED_COUNT=$(grep -v "^#" as_nested.vcf | wc -l | tr -d ' ')
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
AS_NESTED_PS=$(grep -v "^#" as_nested.vcf | awk -F'\t' '$8 ~ /LV=1/ && $8 ~ /PS=/' | wc -l | tr -d ' ')
is "$AS_NESTED_PS" "1" "-A flag: nested variant has PS tag"

rm -f as_nested.vg as_nested.gam as_nested.pack as_nested.vcf

# Test: -A flag on triple nested graph with gref paths
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/triple_nested.gfa > as_triple.gfa 2>/dev/null
vg sim -x as_triple.gfa -P "a#1#y0#0" -n 200 -l 2 -a -s 302 > as_triple.gam
vg sim -x as_triple.gfa -P "a#2#y1#0" -n 200 -l 2 -a -s 303 >> as_triple.gam
vg pack -x as_triple.gfa -g as_triple.gam -o as_triple.pack
vg call as_triple.gfa -k as_triple.pack -A -P gref_x > as_triple.vcf 2>/dev/null

# Should produce variants at multiple nesting levels
AS_TRIPLE_COUNT=$(grep -v "^#" as_triple.vcf | wc -l | tr -d ' ')
AS_TRIPLE_HAS_VARIANTS=$(if [ "$AS_TRIPLE_COUNT" -ge 3 ]; then echo "1"; else echo "0"; fi)
is "$AS_TRIPLE_HAS_VARIANTS" "1" "-A flag: triple nested produces at least 3 variants"

# Verify all variants have LV tags
AS_TRIPLE_LV=$(grep -v "^#" as_triple.vcf | grep -c "LV=")
is "$AS_TRIPLE_LV" "$AS_TRIPLE_COUNT" "-A flag: all triple nested variants have LV tags"

# Verify PS tags on nested variants (LV > 0)
AS_TRIPLE_NESTED=$(grep -v "^#" as_triple.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | wc -l | tr -d ' ')
AS_TRIPLE_NESTED_PS=$(grep -v "^#" as_triple.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/ && $8 ~ /PS=/' | wc -l | tr -d ' ')
is "$AS_TRIPLE_NESTED_PS" "$AS_TRIPLE_NESTED" "-A flag: all nested variants (LV>0) have PS tags"

rm -f as_triple.gfa as_triple.gam as_triple.pack as_triple.vcf

# =============================================================================
# RC, RS, RD tag tests (reference coordinate tags for nested snarls)
# =============================================================================

# Test: RC, RS, RD headers are present
vg paths --compute-gref -Q x --min-gref-len 1 -x nesting/triple_nested.gfa > rc_test.gfa 2>/dev/null
vg sim -x rc_test.gfa -P "a#1#y0#0" -n 100 -l 2 -a -s 400 > rc_test.gam 2>/dev/null
vg sim -x rc_test.gfa -P "a#2#y1#0" -n 100 -l 2 -a -s 401 >> rc_test.gam 2>/dev/null
vg pack -x rc_test.gfa -g rc_test.gam -o rc_test.pack 2>/dev/null
vg call rc_test.gfa -k rc_test.pack --top-down -P gref_x > rc_test.vcf 2>/dev/null

# Check for RC, RS, RD headers
RC_HEADER=$(grep -c "##INFO=<ID=RC" rc_test.vcf)
is "$RC_HEADER" "1" "RC header is present in VCF"

RS_HEADER=$(grep -c "##INFO=<ID=RS" rc_test.vcf)
is "$RS_HEADER" "1" "RS header is present in VCF"

RD_HEADER=$(grep -c "##INFO=<ID=RD" rc_test.vcf)
is "$RD_HEADER" "1" "RD header is present in VCF"

# Check that all variants have RC, RS, RD tags
RC_COUNT=$(grep -v "^#" rc_test.vcf | wc -l | tr -d ' ')
RC_TAG_COUNT=$(grep -v "^#" rc_test.vcf | grep -c "RC=")
is "$RC_TAG_COUNT" "$RC_COUNT" "All variants have RC tag"

RS_TAG_COUNT=$(grep -v "^#" rc_test.vcf | grep -c "RS=")
is "$RS_TAG_COUNT" "$RC_COUNT" "All variants have RS tag"

RD_TAG_COUNT=$(grep -v "^#" rc_test.vcf | grep -c "RD=")
is "$RD_TAG_COUNT" "$RC_COUNT" "All variants have RD tag"

# Check that top-level variant has RC pointing to its own contig
TOP_RC=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/' | grep -o "RC=[^;]*" | cut -d= -f2)
TOP_CHROM=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/ {print $1}')
is "$TOP_RC" "$TOP_CHROM" "Top-level variant RC equals its own CHROM"

# Check that nested variants point to top-level's coordinates
# All nested variants should have RC=gref_x (the top-level reference)
NESTED_RC_X=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | grep -c "RC=gref_x")
NESTED_COUNT=$(grep -v "^#" rc_test.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | wc -l | tr -d ' ')
is "$NESTED_RC_X" "$NESTED_COUNT" "All nested variants have RC=x (top-level contig)"

rm -f rc_test.gfa rc_test.gam rc_test.pack rc_test.vcf


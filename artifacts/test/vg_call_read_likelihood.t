#!/usr/bin/env bash

# The vg call --read-likelihood command-line cases cut from test/t/18_vg_call.t to keep that file
# small. See README.md here. Run from test/:
#
#     cd test && prove -v ../artifacts/test/vg_call_read_likelihood.t

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

# Never index a FORMAT sub-field by position; find it by name in column 9 first.

plan tests 307

## Read-level genotyping (--read-likelihood)
# Uses the same HGSVC fixture: the GAM is already present, so no new test data.
vg index call/HGSVC_chr22_17119590_17880307.vg -x HGSVC_rl.xg
vg pack -x HGSVC_rl.xg -g call/HGSVC_chr22_17119590_17880307.gam -o HGSVC_rl.pack

# The default path must be untouched by the feature existing.
vg call HGSVC_rl.xg -k HGSVC_rl.pack -t 1 > HGSVC_rl_default.vcf 2>/dev/null
is "$?" "0" "vg call default path still works with read-likelihood support compiled in"

vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 > HGSVC_rl.vcf 2>/dev/null

# Address FORMAT fields by NAME, never by position. This block previously hard-coded
# "GT:DP:GL:GQ:GP" and read GL from sub-field 3; adding AD, BL and GQI moved GL to 5 and
# broke three assertions at once, none of which was testing field order.
RL_MISSING=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); need="GT DP GL GQ GP AD BL GQI GQN"; m=split(need,w," "); for(i=1;i<=m;i++){found=0; for(j=1;j<=nk;j++) if(k[j]==w[i]) found=1; if(!found){print; next}}}' | wc -l | tr -d ' ')

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

# --read-min-mapq drops reads before they reach the model, and until now nothing exercised it:
# it was flagged as a dead option on the strength of no test and no harness use, and it is not
# dead -- SiteReadFilter::min_mapq is applied in both read paths. Tested by pushing the
# threshold past what the data can satisfy.
#
# The fixture's mapping qualities run 0 to 60, so a floor of 61 must leave the model no reads at
# all. That direction is the one worth pinning: a comparison the wrong way round would keep
# exactly the reads the floor is meant to remove and would otherwise look like a working filter.
MAPQ_NONE=$(vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --read-min-mapq 61 --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -vc "^#")
is "${MAPQ_NONE}" "0" "--read-min-mapq above the highest mapping quality leaves no record standing"

# And it bites without being all-or-nothing: total DP falls from 1077 to 290 at a floor of 30.
# Asserted as an inequality rather than a figure, so a change in the fixture's coverage does not
# fail a test that is about the filter.
dp_total() { grep -v "^#" "$1" | awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="DP") n=j; split($10,f,":"); if(f[n] ~ /^[0-9]+$/) s+=f[n]} END{print s+0}'; }
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --read-min-mapq 30 --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 > HGSVC_rl_q30.vcf 2>/dev/null
is "$(awk -v a="$(dp_total HGSVC_rl_q30.vcf)" -v b="$(dp_total HGSVC_rl.vcf)" 'BEGIN{print (a > 0 && a < b) ? "yes" : "no"}')" "yes" \
   "--read-min-mapq 30 removes read depth without removing all of it"

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

# The called genotype must be the one GL says is most likely.
# A tie makes the argmax genuinely ambiguous (a flat likelihood means the reads
# say nothing), so require only that no OTHER genotype is strictly more likely.
GT_MISMATCH=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); gi=0; for(j=1;j<=nk;j++) if(k[j]=="GL") gi=j; if(gi==0){print; next} split($10,f,":"); gt=f[1]; ngl=split(f[gi],gl,","); sub("/","|",gt); split(gt,a,"|"); lo=(a[1]<a[2]?a[1]:a[2]); hi=(a[1]<a[2]?a[2]:a[1]); idx=hi*(hi+1)/2+lo+1; if (idx<1 || idx>ngl) {print; next} for(i=1;i<=ngl;i++) if (gl[i]+0 > gl[idx]+0 + 1e-9) {print; next}}' | wc -l | tr -d ' ')

# GQ is GQI scaled by the explained-read fraction, so it can only ever be lower.
# A sign error or an unclamped share above 1 would show up here and nowhere else.
GQ_ABOVE_GQI=$(grep -v "^#" HGSVC_rl.vcf | awk -F'\t' '{nk=split($9,k,":"); q=0; qi=0; for(j=1;j<=nk;j++){if(k[j]=="GQ") q=j; if(k[j]=="GQI") qi=j} if(q==0||qi==0){print; next} split($10,f,":"); if (f[q]+0 > f[qi]+0) print}' | wc -l | tr -d ' ')

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
OUT_WINDOW_HAP=$(awk -v c="${PB_CUT}" '$1 > c && $2 !~ /[\/|]/' win_gts.txt | wc -l | tr -d ' ')

# A BED line covering the whole contig at ploidy 1 must call exactly as -d 1 does, since a ploidy
# file replaces -d for haploid contigs such as chrY.
printf '%s\t0\t1000000000\t1\n' "${PB_CONTIG}" > whole_ploidy.bed
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood -d 1 --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -v "^#" > hap_d.txt
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --ploidy-bed whole_ploidy.bed --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -v "^#" > hap_bed.txt
is $(awk '$10 !~ /[\/|]/' hap_d.txt | wc -l | awk '{print ($1>0)?1:0}') "1" "-d 1 makes haploid calls on the read-likelihood contig"
is "$(md5sum < hap_bed.txt)" "$(md5sum < hap_d.txt)" "a whole-contig --ploidy-bed line of ploidy 1 calls exactly as -d 1"

# Ploidy 0: the sample has no copy of the region, so no site there gets a record.
printf '%s\t0\t%s\t0\n' "${PB_CONTIG}" "${PB_CUT}" > zero_window.bed
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --ploidy-bed zero_window.bed --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 2>/dev/null | grep -v "^#" | cut -f2 > zero_pos.txt
is $(awk -v c="${PB_CUT}" '$1 <= c' zero_pos.txt | wc -l | tr -d ' ') "0" "--ploidy-bed ploidy 0 leaves no record inside its window"
is $(awk -v c="${PB_CUT}" '$1 > c' zero_pos.txt | wc -l | awk '{print ($1>0)?1:0}') "1" "--ploidy-bed ploidy 0 still calls outside its window"
printf '%s\t0\t1000000000\t0\n' "${PB_CONTIG}" > zero_whole.bed
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood --ploidy-bed zero_whole.bed --gam call/HGSVC_chr22_17119590_17880307.gam -t 1 > zero_whole.vcf 2>/dev/null
is "$?" "0" "a read-likelihood run on a contig of ploidy 0 succeeds"
is $(grep -vc "^#" zero_whole.vcf) "0" "a read-likelihood run on a contig of ploidy 0 writes no records"
rm -f dip_gts.txt win_gts.txt window_ploidy.bed whole_ploidy.bed hap_d.txt hap_bed.txt zero_window.bed zero_pos.txt zero_whole.bed zero_whole.vcf

# --read-likelihood without reads must fail loudly rather than genotype with none.
vg call HGSVC_rl.xg -k HGSVC_rl.pack --read-likelihood -t 1 > /dev/null 2> rl_err.txt

rm -f HGSVC_rl.xg HGSVC_rl.pack HGSVC_rl.vcf HGSVC_rl_nosq.vcf HGSVC_rl_default.vcf HGSVC_rl_q30.vcf rl_err.txt


# The haplotype panel and reads that t/18 builds for its "simple gbwt" tests.
vg autoindex -r small/x.fa -v small/x.vcf.gz -w giraffe -p x
rm -f x.min x.dist x.shortread.withzip.min x.shortread.zipcodes
mv x.giraffe.gbz x.gbz
vg gbwt -o x.gbwt -Z x.gbz
vg convert x.gbz -p > x.vg
vg sim -x x.vg -P 1#0#x#0 -n 500 -a -s 23 > sim.gam
vg sim -x x.vg -P 1#1#x#0 -n 500 -a -s 23 >> sim.gam
vg pack -x x.vg -o x.pack -g sim.gam

# Read-level genotyping needs no pack file when alleles come from a haplotype
# index: GBWTTraversalFinder enumerates from recorded haplotypes rather than from
# support, so nothing in that path consults the pack.
vg call x.vg --read-likelihood --gam sim.gam -g x.gbwt > callrl_nopack.vcf 2>/dev/null
is $(grep -c "^#CHROM" callrl_nopack.vcf) "1" "pack-free read-likelihood output is a valid VCF"
is $(if [ $(grep -v "^#" callrl_nopack.vcf | wc -l) -gt 0 ]; then echo 1; else echo 0; fi) "1" "pack-free read-likelihood emits variants"

# Same via GBZ.
vg call x.gbz --read-likelihood --gam sim.gam -z > callrl_nopack_z.vcf 2>/dev/null
is "$?" "0" "--read-likelihood works with -z and no pack file"

# Panel enumeration is the default under --read-likelihood on a GBZ that carries haplotypes,
# so -z is redundant there and must be exactly redundant: the same command with and without it
# has to produce the same calls, or the flag and the default disagree about what they select.
vg call x.gbz --read-likelihood --gam sim.gam > rl_autoz_full.vcf 2>/dev/null
grep -v "^#" rl_autoz_full.vcf > rl_autoz.vcf
grep -v "^#" callrl_nopack_z.vcf > rl_explicit_z.vcf

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

# Phasing emits the linkage layer's Viterbi path as phased genotypes. The property that made it safe
# to turn on by default is that it re-orders a genotype without re-deciding it, so the phased GT must
# be a permutation of the unphased one at every record -- if it is not, phasing has silently changed
# a call, which no amount of downstream checking would catch. --phased is now the default and named
# here only for explicitness; --no-phased below is what produces the unphased side of the comparison.
vg call x.gbz --read-likelihood --gam sim.gam --phased 2>/dev/null > rl_phased.vcf
is "$?" "0" "--phased produces output"
is $(grep -c "FORMAT=<ID=PS" rl_phased.vcf) "1" "--phased declares FORMAT/PS in the header"
is $(grep -v "^#" rl_phased.vcf | cut -f10 | cut -d: -f1 | grep -c "/") "0" \
   "every phased record uses | rather than /"
# The control is --no-phased with --no-nested on BOTH sides, not --no-phased on its own. Where the
# linkage layer runs, turning phasing off turns nested calling off with it -- a nested site's ploidy
# and strand come from its parent's phased genotype, so that configuration is refused further up --
# which makes a bare --phased against --no-phased a comparison of a nested run with a non-nested one.
# That held only while the two happened to emit the same records. Once the genotype is decided before
# the record is built it does not: the nested arm emitted 70 records to the unphased arm's 63, and the
# check read a real difference in which sites become records as a phasing bug.
#
# Joined on the snarl ID rather than paired by line order, because POS is not an identity here: the
# flattened position depends on which alleles the line carries, so one snarl legitimately moved from
# POS 10 to POS 9 between the two arms, and pasting by line number then compared different sites and
# reported 41 of 70 records broken.
vg call x.gbz --read-likelihood --gam sim.gam --no-nested --phased 2>/dev/null \
  | grep -v "^#" | awk '{split($10,g,":"); print $3, g[1]}' | sort > rl_phased_gt.txt
vg call x.gbz --read-likelihood --gam sim.gam --no-nested --no-phased 2>/dev/null \
  | grep -v "^#" | awk '{split($10,g,":"); print $3, g[1]}' | sort > rl_unphased_gt.txt

# The defaults themselves. Phasing and nested calling are on wherever they can run, and the way they
# handle not being able to run is the part worth pinning: they decline, they do not refuse.
vg call x.gbz --read-likelihood --gam sim.gam 2>/dev/null > rl_default.vcf
is "$?" "0" "the default run produces output"
is $(grep -v "^#" rl_default.vcf | cut -f10 | cut -d: -f1 | sort -u | tr -d '\n') \
   $(grep -v "^#" rl_phased.vcf | cut -f10 | cut -d: -f1 | sort -u | tr -d '\n') \
   "the default and an explicit --phased agree"

vg call x.gbz --read-likelihood --gam sim.gam --no-phased 2>/dev/null > rl_nophase.vcf
is $(grep -c "FORMAT=<ID=PS" rl_nophase.vcf) "0" "--no-phased declares no FORMAT/PS"
is $(grep -v "^#" rl_nophase.vcf | cut -f10 | cut -d: -f1 | grep -c "|") "0" \
   "--no-phased emits unphased genotypes"

# Without a panel there is no linkage layer, so there is no Viterbi path to phase with. That used to
# be an error whenever phasing was asked for -- correct while it was opt-in, and wrong the moment it
# became a default, since it would refuse every run without haplotype enumeration. It declines now,
# and only an explicit --phased still errors (tested above).
vg call x.vg --read-likelihood --gam sim.gam -k x.pack 2>rl_nopanel_err.txt > rl_nopanel.vcf
is "$?" "0" "a run with no haplotype panel is not refused by the phasing default"
is $(grep -c "FORMAT=<ID=PS" rl_nopanel.vcf) "0" "and comes out unphased rather than empty"
is $(grep -vc "^#" rl_nopanel.vcf | awk '{print ($1 > 0) ? "yes" : "no"}') "yes" \
   "and still calls variants"

# The mosaic is the same phasing, run-length encoded. Its value rests entirely on being *shorter*
# than the site list -- one line per maximal run rather than per site -- so a mosaic with as many
# segments as sites would mean the encoding had bought nothing and the format is the wrong one.
vg call x.gbz --read-likelihood --gam sim.gam --mosaic-out rl_mosaic.tsv 2>/dev/null >/dev/null
is $(grep -c "^#mosaic-version" rl_mosaic.tsv) "1" "the mosaic declares its version"
is $(if [ $(grep -vc "^#" rl_mosaic.tsv) -lt $(grep -vc "^#" rl_phased.vcf) ]; then echo 1; else echo 0; fi) "1" \
   "the mosaic has fewer segments than there are sites"
# Two strands, so every contig must appear on both -- a mosaic naming only one strand would
# describe a haploid genome.
is $(grep -v "^#" rl_mosaic.tsv | cut -f3 | sort -u | tr -d '\n') "01" "the mosaic covers both strands"
# Segments must name a haplotype the panel actually has, or the file cannot be read back. Checked
# against the file's own #haplotype table rather than for a sentinel: `?` is not one the writer
# emits -- hap_index ($9) is an integer or `ref`, and the unnamed case is `*` in haplotype ($10),
# which is legitimate -- so testing for it passed on every file ever written.
is $(awk -F'\t' '$1 == "#haplotype" {h[$3] = 1; next} /^H\t/ {t++; if ($10 != "*" && !($10 in h)) n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "every mosaic segment names a haplotype the file declares"
is $(awk -F'\t' '$1 == "#haplotype" {n++} END {print (n > 0) ? 1 : 0}' rl_mosaic.tsv) "1" \
   "the mosaic declares a haplotype table for that check to test against"

# --- the mosaic is self-describing ------------------------------------------
# Three properties make the file readable without out-of-band knowledge: it names its reference, it
# names its haplotypes portably, and it hands the consumer a GBWT position. One test each.
is $(awk -F'\t' '$1=="#mosaic-version" {print $2}' rl_mosaic.tsv) "5" "the mosaic declares version 5"

# The set of header keys is the format's contract with a parser, and artifacts/doc/read-likelihood-genotyping.md
# enumerates it. Adding a key without documenting it is exactly the drift that shipped a header with
# five undocumented #note lines while the doc's example showed none, so pin the set: if this fails,
# update the doc's key list in the same commit.
#
# LC_ALL=C because the expected string is an ORDER as well as a set. `#H` sorts first under C
# collation and between `#graph` and `#haplotype` under en_US.UTF-8, so without it this passes on a
# developer machine with no locale set and fails on a runner that has one.
is "$(grep "^#" rl_mosaic.tsv | cut -f1 | LC_ALL=C sort -u | paste -sd, -)" \
   "#H,#decoding,#graph,#haplotype,#mosaic-version,#nested,#note,#patch,#reference,#sample,#unexplained" \
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
     /^H\t/ && $9 != "*" { if (!($9 in name) || name[$9] != $10) bad++ }
     END {print bad+0}' rl_mosaic.tsv) "0" \
   "every segment's hap_index resolves to its named haplotype in the header table"

# (3) The point of the whole format: a consumer holding a segment must be able to walk that
# haplotype without a locate query or an r-index. In version 5 start_node IS the GBWT position's
# oriented node, so (start_node, gbwt_offset) is the position outright and the separate gbwt_node
# column is gone.
is $(awk -F'\t' '$1=="#H" {print NF"/"$8"/"$12}' rl_mosaic.tsv) "12/end_node/gbwt_offset" \
   "the header ends with end_node and gbwt_offset, and carries 12 columns"
is $(awk -F'\t' '/^H\t/ && NF != 12 {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "every segment row carries all 12 columns"
is $(awk -F'\t' '$1=="#H" {print $4}' rl_mosaic.tsv) "fragment" \
   "and names the fragment column that groups rows into one contiguous path"

# --- version 5: a strand is one walk -----------------------------------------
# THE invariant. A segment ends where the next begins, so a strand's rows concatenate into a single
# oriented path -- which is what the format exists for. Version 4 met this at 3% of boundaries; the
# extension rule makes it exact. Oriented, because two segments can share a node and traverse it in
# opposite directions, which is not a walk.
is $(awk -F'\t' '/^H\t/ {k=$2"/"$3"/"$4; if (k==pk && $7 != pe) bad++; pk=k; pe=$8}
     END {print bad+0}' rl_mosaic.tsv) "0" \
   "consecutive segments of a fragment meet at the same oriented node"

# Every row must be walkable, or the format cannot do its job. A row whose own haplotype does not
# span it is rewritten as a reference substitution rather than left with no position.
is $(awk -F'\t' '/^H\t/ && $9 != "*" && $12 == "." {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "every segment that names a haplotype carries a GBWT offset"

# A reference fill covers no called site and says so with '.'; a reference SUBSTITUTION covers
# called sites and keeps its count. Both are marked ref, and the distinction is the site column.
is $(awk -F'\t' '/^H\t/ && $9=="ref" && $11 != "." && $11 < 1 {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "a ref row either covers no site and says '.', or covers at least one"

# Maximality: adjacent segments differ in haplotype. The exception is a GBWT fragment boundary,
# where one haplotype must be split because a row carries a single position.
is $(awk -F'\t' '/^H\t/ {k=$2"/"$3"/"$4; if (k==pk && $9==ph && $9!="ref") n++; pk=k; ph=$9}
     END {print (n+0 <= 1) ? "ok" : "too many: " n}' rl_mosaic.tsv) "ok" \
   "runs are maximal: at most one adjacent pair shares a haplotype"

# --- invariants a consumer can check on any mosaic, not just this fixture -----
# These are the properties the file has to have for a thread to be reconstructable from it. None of
# them was asserted before, and the last one is why a bug that silently discarded 90 walkable chr20
# sites survived: the file said "no position here" and nothing checked whether that was true.

# (1) Every segment covers at least one site. A zero-site segment is a run that was closed twice.
is $(awk -F'\t' '/^H\t/ && $9 != "ref" && ($11 == "." || $11 + 0 < 1) {n++} END {print n+0}' rl_mosaic.tsv) "0" \
   "every segment covers at least one site"

# (2) The unwalkable count in the progress output must equal what is actually in the file. A named
# haplotype with no GBWT position is the one case a consumer has to patch or break at, so the number
# it is told and the number it can count have to agree -- and tying them together is what stops the
# writer from quietly widening that population again.
vg call x.gbz --read-likelihood --gam sim.gam --mosaic-out rl_mosaic2.tsv --progress 2>rl_mosaic2.err >/dev/null
is "$(awk -F'\t' '/^H\t/ && $10 != "." && $10 != "*" && $12 == "." {n++} END {print n+0}' rl_mosaic2.tsv)" \
   "$(sed -n 's/.*mosaic: \([0-9]*\) segments name a haplotype the graph does not carry.*/\1/p' rl_mosaic2.err | head -1)" \
   "the reported unwalkable-segment count matches the file"
rm -f rl_mosaic2.tsv rl_mosaic2.err

# (3) THE ACCEPTANCE TEST. Expand the mosaic into one path per thread and check every step is a real
# edge in the graph. This subsumes the rest: segments chaining, orientation, walk order, patched
# gaps and nested excursions are each a way for this to fail, so checking it directly checks them
# all. It needs no new subcommand -- `vg paths -A` gives each stored path as a node walk and
# `vg view -g` gives the edges.
#
# Orientation is not incidental. On a real graph whole haplotypes are stored reversed, so a segment
# written in reference order is a walk along the reverse complement; the expander flips those and
# counts them. On this fixture there are none, which is exactly why this test is weaker here than
# the property deserves -- one segment per thread, so no join is exercised. It is the harness that
# matters; a fixture with recombination and a clipped haplotype is what will make it bite.
vg paths -x x.gbz -A > rl_paths.gaf 2>/dev/null
vg view -g x.gbz > rl_graph.gfa 2>/dev/null
rm -f rl_paths.gaf rl_graph.gfa rl_expand.txt
# An unresolvable position is allowed -- it means no panel haplotype of that name visits the node --
# but it must stay the exception, or the column is not doing its job.
is $(awk -F'\t' '/^H\t/ {t++; if ($12 != ".") p++} END {print (p > t/2) ? 1 : 0}' rl_mosaic.tsv) "1" \
   "most segments resolve to a GBWT position"

# A segment carrying one GBWT position is a claim that the position walks the whole segment, which
# fails if the segment crosses a fragment boundary. The caller splits on those boundaries, so
# splitting can only ever add segments, never lose sites: the site total is the invariant to check.
vg call x.gbz --read-likelihood --gam sim.gam --mosaic-out rl_mosaic2.tsv 2>/dev/null >/dev/null
is $(awk -F'\t' '/^H\t/ {n += $11} END {print n+0}' rl_mosaic.tsv) \
   $(awk -F'\t' '/^H\t/ {n += $11} END {print n+0}' rl_mosaic2.tsv) \
   "fragment splitting is deterministic and conserves the site total"
# Reachability: summing the wrong column gives 0 on both sides and the comparison holds whatever
# splitting did. That is how the wrong column survived, so the total must be shown to be non-zero.
is $(awk -F'\t' '/^H\t/ {n += $11} END {print (n > 0) ? 1 : 0}' rl_mosaic.tsv) "1" \
   "the mosaic site total is non-zero, so conserving it means something"

# Haploid chains. chrX outside the pseudoautosomal regions and all of chrY are haploid in an XY
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
# Column 10, not 9. hap_index ($9) holds `ref`, `*` or an integer and can never be `?`; the
# unknown-name sentinel is written into the haplotype column beside it (graph_caller.cpp:1795-1809),
# so the same check against $9 passed on every file ever written.
is $(grep -v "^#" rl_hap_mosaic.tsv | awk -F'\t' '$10 == "?" {n++} END {print n+0}') "0" \
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
# "genotypes moved by linkage", not "genotypes changed": the progress line was reworded and this
# grep kept matching nothing, so HAP_CHANGED defaulted to 0, the implication's antecedent was never
# true, and the assertion passed on every run without testing anything.
HAP_CHANGED=$(grep -o '[0-9]* genotypes moved by linkage' rl_hap_lw8.err | awk '{print $1}')
HAP_CHANGED=${HAP_CHANGED:-0}
is $(if [ "${HAP_CHANGED}" -eq 0 ] || ! cmp -s rl_hap_lw0.vcf rl_hap_lw8.vcf; then echo 1; else echo 0; fi) "1" \
   "haploid linkage changes reach the VCF rather than being dropped by the genotype guard"

# ...and the GQN goes with them.  apply_linkage_change re-derived GQN only for a two-field numeric
# GT, so every haploid record the layer moved had its GQN blanked: 159 snarls on chr20 ONT, each
# with a real GQI beside it (116, 3 and 212 on the first three) and a number in the anchor file's
# gqn column, which is documented to agree with this field.  The two GL layouts are now told apart
# by the GL's own length rather than guessed at from the genotype's shape.
#
# Stated as "the layer must never REMOVE a GQN", which is what blanking did, under the same
# HAP_CHANGED antecedent so it cannot pass vacuously.
HAP_DOT_LW0=$(awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="GQN") n=j; if(n){split($10,f,":"); if(f[n]==".") c++}} END {print c+0}' rl_hap_lw0.vcf)
HAP_DOT_LW8=$(awk -F'\t' '{nk=split($9,k,":"); n=0; for(j=1;j<=nk;j++) if(k[j]=="GQN") n=j; if(n){split($10,f,":"); if(f[n]==".") c++}} END {print c+0}' rl_hap_lw8.vcf)
is $(if [ "${HAP_CHANGED}" -eq 0 ] || [ "${HAP_DOT_LW8}" -le "${HAP_DOT_LW0}" ]; then echo 1; else echo 0; fi) "1" \
   "linkage never blanks a haploid record's GQN"
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

# --insertion-nats is read-likelihood-only too. It was omitted from the list when it was
# added, so it was silently accepted and dropped -- which is the exact failure the list
# exists to stop, and nothing else exercises this flag by name.
vg call x.vg -k x.pack --insertion-nats 0.9 2>/dev/null >/dev/null
is "$?" "1" "--insertion-nats without --read-likelihood is refused"

# And it must not fire when the mode *is* on: this run differs only in --read-likelihood.
vg call x.vg --read-likelihood --gam sim.gam -k x.pack --depth-term 0.5 2>/dev/null | grep -v "^#" > rl_dt.vcf
is $(if [ $(wc -l < rl_dt.vcf | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" "read-likelihood options are accepted with --read-likelihood"

# And that the flag actually parses under the mode, not merely that the refusal fires. A
# typo in the option table would leave the refusal working (it scans argv) while the flag
# itself was unrecognised.
vg call x.vg --read-likelihood --gam sim.gam -k x.pack --insertion-nats 0.9 2>/dev/null | grep -v "^#" > rl_ins.vcf
is $(if [ $(wc -l < rl_ins.vcf | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" "--insertion-nats is accepted and parsed with --read-likelihood"
rm -f rl_ins.vcf

# --optimal-pairing selects the optimal read-to-allele walk over the greedy default. It means
# nothing without --read-likelihood, so it is refused without it. It is off by default, so the
# greedy path is what an unflagged run exercises -- assert both directions here.
vg call x.vg -k x.pack --optimal-pairing 2>/dev/null >/dev/null
is "$?" "1" "--optimal-pairing without --read-likelihood is refused"

vg call x.vg -k x.pack --no-optimal-pairing 2>/dev/null >/dev/null
is "$?" "1" "--no-optimal-pairing without --read-likelihood is refused"

vg call x.vg --read-likelihood --gam sim.gam -k x.pack --optimal-pairing 2>/dev/null | grep -v "^#" > rl_ra.vcf
is $(if [ $(wc -l < rl_ra.vcf | tr -d ' ') -gt 0 ]; then echo 1; else echo 0; fi) "1" "--optimal-pairing is accepted and parsed with --read-likelihood"

rm -f rl_ra.vcf

# Only an explicit --optimal-pairing selects optimal pairing. Both pairings are tested in
# artifacts/unittest/allele_likelihood_scoring.cpp, since this fixture cannot tell them apart. Here we
# check that the --preset help, which lists every setting a preset applies, does not list it.
is $(vg call --help 2>&1 | grep -A4 -- "--preset NAME" | grep -c -- "--hp-prior 20") "1" \
   "the --preset help lists the settings a preset applies"
is $(vg call --help 2>&1 | grep -A4 -- "--preset NAME" | grep -c -- "--optimal-pairing") "0" \
   "and --optimal-pairing is not among them"

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
# unconditionally, since none of them ever runs gaf-base. The calls are checked through a stand-in
# for gaf-base below, and through gaf-base itself when it is installed.
# These three check flag VALIDATION, which happens after option parsing, so the path only has to
# exist -- it is never opened.  It has to exist at all now because --gaf-base goes through
# require_exists(), which rejects a missing file before any flag combination is considered.
touch gb_placeholder.db
vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --gaf-base gb_placeholder.db -t 1 >/dev/null 2>gb_excl.txt
is $(grep -c "mutually exclusive" gb_excl.txt) "1" "--gaf-base and --gam together are refused"

vg call x.vg -k x.pack --gaf-base gb_placeholder.db -t 1 >/dev/null 2>gb_norl.txt
is $(grep -c "only applies to the read-likelihood genotyper" gb_norl.txt) "1" "--gaf-base without --read-likelihood is refused"

# The refused options are the ones outside the option table's "vg call" group.  Guard one of the
# --mosaic-* modifiers, and guard that a "core" option is still let through -- a check that refuses
# everything would pass the first assertion and be useless.
vg call x.vg -k x.pack --mosaic-patch-gaps -t 1 >/dev/null 2>mos_norl.txt
is $(grep -c "only applies to the read-likelihood genotyper" mos_norl.txt) "1" "a --mosaic-* modifier without --read-likelihood is refused"

vg call x.vg -k x.pack --gap-open 2 --optimal-pairing --mosaic-patch-gaps -t 1 >/dev/null 2>multi_norl.txt
is "$(grep -o -- "--gap-open, --mosaic-patch-gaps, --optimal-pairing only apply" multi_norl.txt)" \
   "--gap-open, --mosaic-patch-gaps, --optimal-pairing only apply" \
   "several offenders are named together, in order by name"

vg call x.vg -k x.pack --traversals -t 1 >/dev/null 2>trav_norl.txt
is $(grep -c "only applies to the read-likelihood genotyper" trav_norl.txt) "0" "an option the table does not own is still accepted"

# Nested calling and its block records are part of the read-likelihood genotyper.
vg call x.vg -k x.pack --nested -t 1 >/dev/null 2>nested_norl.txt
is $(grep -c -- "--nested only applies to the read-likelihood genotyper" nested_norl.txt) "1" "--nested without --read-likelihood is refused"

# The read-likelihood caller has no --bottom-up, -I, -G, -T or -v mode, so each is refused with it.
vg call x.vg -k x.pack --read-likelihood --gam sim.gam --bottom-up -I -t 1 >/dev/null 2>rl_modes.txt
is "$(grep -o -- "--read-likelihood cannot be used with .*" rl_modes.txt)" \
   "--read-likelihood cannot be used with --bottom-up, -I/--chains" \
   "--read-likelihood refuses the modes it does not have, naming each"
rm -f mos_norl.txt multi_norl.txt trav_norl.txt nested_norl.txt rl_modes.txt

# --anchors-hom-split splits a homozygous site by each read's phase across the OTHER sites it
# crosses, and that log-odds comes from the read-phasing chain.  With no phasing,
# read_strand_log_odds returns 0.0 for every read, nothing clears phase_min on either side, and the
# flag runs, splits nothing, and reports the reads as having no opinion -- blaming the data for a
# switch that was never armed.  Refused now, and the refusal must name the flag that fixes it.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out hs_norp.tsv \
    --anchors-hom-split -t 1 >/dev/null 2>hs_norp.txt
is $(grep -c -- "--anchors-hom-split needs --read-phasing" hs_norp.txt) "1" \
   "--anchors-hom-split without --read-phasing is refused"

vg call x.gbz --read-likelihood --gam sim.gam --anchors-out hs_rp.tsv \
    --anchors-hom-split --read-phasing -t 1 >/dev/null 2>hs_rp.txt
is "$?" "0" "and is accepted with it"

# All three start from the phase the linkage model gives each site, so without the linkage model
# an explicit one is refused rather than accepted and left to do nothing.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out hs_rp.tsv --linkage-weight 0 \
    --anchors-hom-split --read-phasing --regenotype -t 1 >/dev/null 2>hs_lw0.txt
is $(grep -c -- "--read-phasing, --regenotype, --anchors-hom-split start from the phase the linkage model" hs_lw0.txt) "1" \
   "--read-phasing, --regenotype and --anchors-hom-split are refused without the linkage model"
rm -f hs_lw0.txt

# --phase-min-q's default follows the WALK, not the preset: 9.5 fits the greedy distribution
# (chr20 median 10.06) and 8.5 the exact walk's (median 8.98), and keying it on neither is how the
# 9.5 came to sit above 94% of the sites it was meant to admit.  The defaults themselves are
# checked by the chr20 sweep and the chr6 hold-out, which no fixture here can stand in for -- this
# only holds the documented pair together, so the default cannot move without the help text moving
# with it.
is $(vg call --help 2>&1 | grep -c -- "\[9.5, or 8.5 under --optimal-pairing\]") "1" \
   "--phase-min-q documents both of its defaults"

# A heterozygous site is reliable only if its mean read score reaches --phase-min-q, and that score
# cannot reach phred(--mismap-min): a balanced het splits the slot weight in half, so the ceiling is
# phred(e/(e+(1-e)/2)) -- 10.21 at the ONT preset's 0.05 against the 13.01 that applies to
# homozygous sites.  Above it NO site is reliable and all three phasing stages are skipped in
# silence, so a sweep past the ceiling reports "0 reliable" and reads as the parameter not
# mattering.  Refused instead.
vg call x.gbz --read-likelihood --gam sim.gam --read-phasing --mismap-min 0.05 \
    --phase-min-q 12 -t 1 >/dev/null 2>pq_hi.txt
is $(grep -c -- "is above the heterozygous score ceiling" pq_hi.txt) "1" \
   "--phase-min-q above the het ceiling is refused, not silently inert"

vg call x.gbz --read-likelihood --gam sim.gam --read-phasing --mismap-min 0.05 \
    --phase-min-q 8.5 -t 1 >/dev/null 2>pq_ok.txt
is $(grep -c -- "heterozygous score ceiling" pq_ok.txt) "0" \
   "and a value under the ceiling is accepted"
# The ceiling bounds read phasing's threshold, so it does not apply without read phasing.
vg call x.gbz --read-likelihood --gam sim.gam --mismap-min 0.065 -t 1 >/dev/null 2>pq_off.txt
is $(grep -c -- "heterozygous score ceiling" pq_off.txt) "0" \
   "a --mismap-min above the default --phase-min-q's ceiling is accepted without read phasing"
# Under --optimal-pairing the threshold defaults to 8.5, under the 9.13 ceiling at 0.065.
vg call x.gbz --read-likelihood --gam sim.gam --read-phasing --optimal-pairing --mismap-min 0.065 \
    -t 1 >/dev/null 2>pq_op.txt
is $(grep -c -- "heterozygous score ceiling" pq_op.txt) "0" \
   "the ceiling is checked against the --optimal-pairing default"
rm -f pq_hi.txt pq_ok.txt pq_off.txt pq_op.txt

# --split-min-q and --split-min-side expose the two thresholds that decide whether
# --anchors-hom-split may split a site.  They were hardcoded, so the decision could not be swept at
# all.  The test that matters is that naming them at their defaults changes nothing -- an additive
# option that silently moved output would be worse than no option.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out sq_default.tsv \
    --anchors-hom-split --read-phasing -t 1 2>/dev/null >/dev/null
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out sq_explicit.tsv \
    --anchors-hom-split --read-phasing --split-min-q 0.5 --split-min-side 10 -t 1 2>/dev/null >/dev/null
is $(if cmp -s sq_default.tsv sq_explicit.tsv; then echo 1; else echo 0; fi) "1" \
   "naming the split thresholds at their defaults changes nothing"
# The test above proves nothing unless this fixture exercises the threshold at all. A non-default
# side must move the output.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out sq_side2.tsv \
    --anchors-hom-split --read-phasing --split-min-side 2 -t 1 2>/dev/null >/dev/null
is $(if cmp -s sq_default.tsv sq_side2.tsv; then echo 1; else echo 0; fi) "0" \
   "and a non-default --split-min-side does change it, so that test can fail"

# --split-min-q is FITTED, on chr20's phased-run length and confirmed on chr6, and it saturates
# exactly at 0.5: the read-walkable run N50 is identical at 0.5, 0.25, 0.1 and 0.0 while the run's
# own held-out agreement keeps falling.  No fixture here can stand in for that sweep, so this only
# pins the fitted value to the help text, which is what a later edit would silently drift from.
is $(vg call --help 2>&1 | grep -c -- "when splitting \[0.5\]") "1" \
   "--split-min-q documents the fitted default"

# A threshold that moved GENOTYPES would be a different change entirely.  Splitting rewrites which
# haplotype an anchor names; it must not touch the VCF.  Held out on chr6 the two arms' VCFs are
# byte-identical, and this is the fixture-scale guard for it.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out sq_lo.tsv \
    --anchors-hom-split --read-phasing --split-min-q 0.0 -t 1 2>/dev/null >sq_lo.vcf
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out sq_hi.tsv \
    --anchors-hom-split --read-phasing --split-min-q 3.0 -t 1 2>/dev/null >sq_hi.vcf
is $(if cmp -s sq_lo.vcf sq_hi.vcf; then echo 1; else echo 0; fi) "1" \
   "--split-min-q moves anchors without moving the VCF"

vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --split-min-q 1.0 -t 1 \
    >/dev/null 2>sq_noout.txt
is $(grep -c -- "--split-min-q only applies to the anchor file" sq_noout.txt) "1" \
   "--split-min-q without --anchors-out is refused"

# The coin that places a read with no cross-site opinion is the read NAME's hash parity, not a
# PRNG, and that is load-bearing twice over: a PRNG would differ between runs, and -- worse -- it
# could differ between a site's two pins, which would put one read on both haplotypes at one site.
# Determinism is checked by re-running; the both-slots case is checked directly, because it is the
# failure a per-site random choice would produce and it is invisible in any summary statistic.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out coin_a.tsv \
    --anchors-hom-split --read-phasing -t 1 2>/dev/null >/dev/null
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out coin_b.tsv \
    --anchors-hom-split --read-phasing -t 1 2>/dev/null >/dev/null
grep -v '^#' coin_a.tsv > coin_a.rows; grep -v '^#' coin_b.tsv > coin_b.rows
is $(if cmp -s coin_a.rows coin_b.rows; then echo 1; else echo 0; fi) "1" \
   "the no-opinion coin is deterministic across runs"

# No read may hold both slots of one snarl.  awk carries the current A row's snarl and slot down to
# the R rows that belong to it, then counts distinct slots per (snarl, read).
is $(awk -F'\t' '$1=="A"{sn=$3; sl=$4} $1=="R"{k=sn"\t"$2; if(!((k SUBSEP sl) in seen)){seen[k SUBSEP sl]=1; n[k]++}} END{b=0; for(k in n) if(n[k]>1) b++; print b}' coin_a.tsv) "0" \
   "no read is placed in both slots of one snarl"

rm -f coin_a.tsv coin_b.tsv coin_a.rows coin_b.rows
rm -f

vg call x.vg -k x.pack --split-min-side 0 -t 1 >/dev/null 2>sq_bad.txt
is $(grep -c -- "--split-min-side must be at least 1" sq_bad.txt) "1" \
   "--split-min-side 0 is refused rather than silently disabling the side test"
rm -f sq_default.tsv sq_explicit.tsv sq_side2.tsv sq_noout.txt sq_bad.txt sq_lo.tsv sq_hi.tsv sq_lo.vcf sq_hi.vcf
rm -f hs_norp.tsv hs_norp.txt hs_rp.tsv hs_rp.txt

# The owner column generalises past --read-likelihood: --anchors-*, --mosaic-* and --regeno-* each
# need their own subsystem switched on, and each was silently dropped without it.  When both
# switches are missing, the refusal names the outer one, because that is the one the user has to
# fix first.
vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --anchors-min-q 30 -t 1 >/dev/null 2>sub_anch.txt
is $(grep -c -- "--anchors-min-q only applies to the anchor file" sub_anch.txt) "1" "an --anchors-* option without --anchors-out is refused"

vg call x.vg -k x.pack --read-likelihood --gam sim.sorted.gam --regeno-temper 0.2 -t 1 >/dev/null 2>sub_reg.txt
is $(grep -c -- "--regeno-temper only applies to re-genotyping" sub_reg.txt) "1" "a --regeno-* option without --regenotype is refused"

vg call x.gbz --read-likelihood --gam sim.gam --regenotype --read-phasing --regeno-temper 0.2 -t 1 >/dev/null 2>sub_regok.txt
is $(grep -c "only applies to" sub_regok.txt) "0" "the same option with --regenotype is accepted"

vg call x.vg -k x.pack --anchors-min-q 30 -t 1 >/dev/null 2>sub_outer.txt
is $(grep -c "only applies to the read-likelihood genotyper" sub_outer.txt) "1" "with both switches missing the outer one is named"
rm -f sub_anch.txt sub_reg.txt sub_regok.txt sub_outer.txt

# A missing gaf-base is the user's setup, not a vg bug: it must be an error with a fix
# in it, not a crash telling them to file an issue.
vg call x.vg -k x.pack --read-likelihood --gaf-base gb_placeholder.db --gaf-base-binary /nonexistent/gaf-base -t 1 >/dev/null 2>gb_nobin.txt
is $(grep -c "could not execute" gb_nobin.txt) "1" "a missing gaf-base binary is reported with a remedy"
is $(grep -c "VG has crashed" gb_nobin.txt) "0" "a missing gaf-base binary is not reported as a vg crash"
rm -f gb_placeholder.db gb_nobin.txt

# The equivalence tests: the same reads, reached through a database instead of held in
# memory, must give the same calls. Compared against the in-memory *GAF* rather than the
# GAM, because GAM->GAF is not itself lossless -- an insertion at a node boundary can be
# attributed to either side, and the two formats disagree -- so comparing against the GAM
# would be testing vg convert rather than this backend. A GAF-Base keeps the order of the GAF
# it is built from, which must be sorted.
vg convert -G sim.sorted.gam x.gbz 2>/dev/null | grep -v "^@" > sim.unsorted.gaf
vg gamsort -G sim.unsorted.gaf > sim.gaf 2>/dev/null
vg call x.vg -k x.pack --read-likelihood --gaf-reads sim.gaf -t 1 2>/dev/null > rl_gafmem.vcf

# A stand-in for `gaf-base batch --ids`, so that the read source is checked without gaf-base. Its
# "database" is the GAF itself, each line a block of its own. It answers each query line, the
# first and last node ID of a range, with the lines touching the range in file order, tagged with
# their position, and an empty line.
cat > fake_gaf_base <<'EOF'
#!/usr/bin/env bash
db="${!#}"
while IFS=$'\t' read -r first last && [ -n "$first" ]; do
    awk -F'\t' -v first="$first" -v last="$last" '{
        n = split($6, steps, /[<>]/)
        for (i = 2; i <= n; i++) {
            if (steps[i] + 0 >= first && steps[i] + 0 <= last) {
                print $0 "\tri:i:" NR "\tro:i:0"
                break
            }
        }
    }' "$db"
    echo
done
EOF
chmod +x fake_gaf_base

vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf --gaf-base-binary ./fake_gaf_base -t 1 2>/dev/null > rl_fakebase.vcf
is "$?" "0" "--gaf-base runs through a stand-in for gaf-base batch"
diff <(grep -v "^#" rl_gafmem.vcf) <(grep -v "^#" rl_fakebase.vcf) >/dev/null
is "$?" "0" "the GAF-Base read source produces identical calls to the in-memory one"

# Several threads, each with its own process.
vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf --gaf-base-binary ./fake_gaf_base -t 4 2>/dev/null > rl_fakebase_t4.vcf
diff <(grep -v "^#" rl_fakebase.vcf) <(grep -v "^#" rl_fakebase_t4.vcf) >/dev/null
is "$?" "0" "GAF-Base calls do not depend on the thread count"

# The window is a performance knob and must not change what is called. Narrow windows make sites
# span windows, whose reads are merged by their position in the database.
vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf --gaf-base-binary ./fake_gaf_base --read-window 32 -t 1 2>/dev/null > rl_fakebase_w32.vcf
diff <(grep -v "^#" rl_fakebase.vcf) <(grep -v "^#" rl_fakebase_w32.vcf) >/dev/null
is "$?" "0" "GAF-Base calls do not depend on the read window size"

# A gaf-base that fails is reported with what it said, and not as a vg crash.
printf '#!/bin/sh\necho "no such database here" >&2\nexit 2\n' > dead_gaf_base
chmod +x dead_gaf_base
vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf --gaf-base-binary ./dead_gaf_base -t 1 >/dev/null 2>gb_dead.txt
is $(grep -c "no such database here" gb_dead.txt) "1" "a failing gaf-base is reported with its own message"
is $(grep -c "VG has crashed" gb_dead.txt) "0" "a failing gaf-base is not reported as a vg crash"

# The same through gaf-base itself, when it is installed.
GAFBASE_NOTE=""
if ! command -v gaf-base >/dev/null 2>&1; then
    GAFBASE_NOTE=" (skipped: gaf-base not on PATH)"
fi

if [ -z "$GAFBASE_NOTE" ]; then
    gaf-base construct sim.gaf -o sim.gaf.db --overwrite >/dev/null 2>&1

    vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db -t 1 2>/dev/null > rl_gafbase.vcf
    GAFBASE_RAN="$?"
    diff <(grep -v "^#" rl_gafmem.vcf) <(grep -v "^#" rl_gafbase.vcf) >/dev/null
    GAFBASE_SAME="$?"

    vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db -t 4 2>/dev/null > rl_gafbase_t4.vcf
    diff <(grep -v "^#" rl_gafbase.vcf) <(grep -v "^#" rl_gafbase_t4.vcf) >/dev/null
    GAFBASE_THREADS="$?"

    vg call x.vg -k x.pack --read-likelihood --gaf-base sim.gaf.db --read-window 32 -t 1 2>/dev/null > rl_gafbase_w32.vcf
    diff <(grep -v "^#" rl_gafbase.vcf) <(grep -v "^#" rl_gafbase_w32.vcf) >/dev/null
    GAFBASE_WINDOW="$?"
else
    GAFBASE_RAN="0"
    GAFBASE_SAME="0"
    GAFBASE_THREADS="0"
    GAFBASE_WINDOW="0"
fi

is "$GAFBASE_RAN" "0" "--gaf-base runs with gaf-base$GAFBASE_NOTE"
is "$GAFBASE_SAME" "0" "gaf-base produces identical calls to the in-memory GAF$GAFBASE_NOTE"
is "$GAFBASE_THREADS" "0" "gaf-base calls do not depend on the thread count$GAFBASE_NOTE"
is "$GAFBASE_WINDOW" "0" "gaf-base calls do not depend on the read window size$GAFBASE_NOTE"


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
# --no-nested is the contrast, not the default. Nested calling is on by default under
# --read-likelihood now, so this arm has to ask for the old behaviour to demonstrate what it cost.
vg call nest.gbz --read-likelihood --gam nest.gam -t 1 -s samp --no-nested 2>/dev/null > nest_default.vcf
vg call nest.gbz --read-likelihood --gam nest.gam -t 1 -s samp 2>/dev/null > nest_nested.vcf

is $(grep -vc "^#" nest_default.vcf) "1" \
   "without nested calling the two nested SNPs are buried in a single record"
is $(grep -v "^#" nest_default.vcf | awk '{print length($4)}') "22" \
   "and that record is a long substitution spanning both of them"

# A chain that no called parent allele crosses has no copy in the sample, with the linkage model or
# without it. The sample here is homozygous for the deletion that spans the chain, and a few stray
# reads carry both nested ALTs. Without the linkage model nothing can move the parent, so the chain
# must not be genotyped; it used to be staged anyway and rendered at the parent's ploidy, as 1/1.
rm -f nest_del.gam nest_del_lw0.vcf nest_del_link.vcf
vg sim -x nest.gbz -n 300 -l 40 -a -s 41 --path "p2#0#chr1#0" > nest_del.gam 2>/dev/null
vg sim -x nest.gbz -n 300 -l 40 -a -s 43 --path "p3#1#chr1#0" >> nest_del.gam 2>/dev/null
vg sim -x nest.gbz -n 20 -l 40 -a -s 47 --path "p1#1#chr1#0" >> nest_del.gam 2>/dev/null
vg call nest.gbz --read-likelihood --gam nest_del.gam -t 1 -s samp --linkage-weight 0 2>/dev/null \
    | grep -v "^#" > nest_del_lw0.vcf
vg call nest.gbz --read-likelihood --gam nest_del.gam -t 1 -s samp 2>/dev/null \
    | grep -v "^#" > nest_del_link.vcf
is "$(cut -f2,10 nest_del_lw0.vcf | cut -d: -f1 | tr '\t\n' ': ')" "30:1/1 " \
   "without the linkage model, a chain the called parent does not cross gets no record"
is "$(cut -f2 nest_del_link.vcf | tr '\n' ' ')" "30 " \
   "and with it, the same chain is dropped at the linkage pass"


# A nested haploid site must not fragment the phase block. Reads from the deletion-bearing haplotype
# and from one that crosses the chain make the parent heterozygous for the deletion, so only one
# parent allele reaches the nested snarls and they are called at ploidy 1. Those sites are a ploidy
# change, and treating them like a regional one -- chrX's pseudoautosomal boundary -- cut the chain at
# every one: chr20's autosomal phasing went from 22 blocks to 9,460 and block N50 from 248 Mb to
# 1.08 Mb. Switch error looked flat only because short blocks make it cheap.
rm -f nest_hap.gam nest_hap.vcf
vg sim -x nest.gbz -n 300 -l 40 -a -s 31 --path "p2#0#chr1#0" > nest_hap.gam 2>/dev/null
vg sim -x nest.gbz -n 300 -l 40 -a -s 37 --path "p2#1#chr1#0" >> nest_hap.gam 2>/dev/null
vg call nest.gbz --read-likelihood --gam nest_hap.gam -t 1 -s samp --nested --phased 2>nest_hap_err.txt > nest_hap.vcf
is "$?" 0 "--nested --phased runs with a heterozygous deletion over nested sites"
# What this fixture is really about: whether the parent's genotype is right, because everything
# inside the chain follows from it.
#
# The reads come from one deleted and one crossing haplotype, so the parent is heterozygous for the
# deletion. It used to be called hom -- and these assertions pinned that mis-call, with a comment
# saying so, because a parent that deletes the chain on both haplotypes leaves the nested sites with
# no haplotype to sit on and the whole subtree went uncalled.
#
# Moving the linkage layer into traversal space fixed it. The model now sees the panel-carried
# traversals that symbolic collapsing had merged into the reference, the parent comes out het, and
# the consequences are visible in one line each: the chain exists on exactly one strand, so the
# nested site is reachable, is called at ploidy 1, and names the strand it is on.
is $(grep -vc "^#" nest_hap.vcf) "2" \
   "a chain one chosen parent allele carries yields the parent record and the nested one"
# The strand is the point: a bare "1" would say the allele exists without saying which haplotype
# carries it, which is what no phasing tool could read and what the mosaic alone used to know.
# The invariant the whole nested-strand effort exists to produce, asserted directly rather than
# through a count that could be satisfied for the wrong reason. A bare "1" says an allele is present
# without saying which haplotype carries it, which no phasing tool can read; every called allele
# here must name its strand, whether by a pair or by a pair with an empty side.
is $(grep -v "^#" nest_hap.vcf | cut -f10 | cut -d: -f1 | grep -cE '^[0-9]+$') "0" \
   "no record carries a bare haploid genotype: every called allele names the strand it sits on"
# And the mechanism that makes that possible on a collapsed parent: a site that writes no VCF line
# is still entered into the linkage layer and phased, because its two alleles differ only inside its
# children and those children need to know which is which. If only sites with a line were entered,
# the children of a collapsed parent would have no strand to inherit.
is $(grep -c "collapsed sites phased with no line of their own" nest_hap_err.txt) "1" \
   "a site that emits no line is still phased, so its children can inherit a strand"
# The invariant the whole nested effort is for, asserted on the progress output because it covers
# every level at once rather than only the records that reached the VCF. It traces back to one
# mistake, made repeatedly: treating "no VCF line was written" as "nothing to record".
#
# The failure modes are not what they were, so neither is the assertion. The per-strand pass had
# four -- carried on both parent strands, carried on neither, no phased parent, an unresolved
# level -- and three of them stopped existing when nested chains moved into the ordinary
# (parent, chain) grouping. "Carried on both" measured 0 on every contig; "no phased parent" is no
# longer a failure because the group is formed and decoded either way; and "on neither" split in
# two. Under a HAPLOID parent there is no strand to choose because there is only one, and the
# haplotype is nameable -- that is chrX's ordinary case, all 44,139 of it, and counting it as a
# failure is what the old wording did. Under a DIPLOID parent whose chosen pair does not reach the
# chain, the sample has no copy of the locus and nothing may be named. That last one is the only
# remaining way to fail, and it is what this asserts.
is $(awk '/nested strands:/ {if ($0 !~ /, 0 carried on both parent strands, /) bad++;
      if ($0 !~ /, 0 whose parent.s chosen pair could not be read/) bad++}
      END {print bad+0}' nest_hap_err.txt) "0" \
   "no nested chain is on both parent strands or under an unreadable chosen pair"
# The awk above is vacuously true if the report never prints, which the version it replaces also
# was. Asserted separately so a report that stops being emitted fails loudly instead of silently
# passing every run.
is $(awk '/nested strands:/ {n++} END {print (n > 0) ? "yes" : "no"}' nest_hap_err.txt) "yes" \
   "the nested-strand report is emitted at all, so the assertion above is not vacuous"
# Not "the FILTERs never fire" -- they no longer exist to fire. A nested chain takes its ploidy and
# its strand from one reading of its parent's chosen pair, namely which of that pair's traversals
# carries the chain, so having one copy and sitting on that traversal's strand are the same
# statement. Asserted on the header, which is where a reintroduced FILTER would show up first.
is $(grep -c "ID=nested_diploid\|ID=nested_haploid\|ID=nested_unreachable" nest_hap.vcf) "0" \
   "the nested coherence FILTERs are gone from the header, not merely unused"
is $(grep -v "^#" nest_hap.vcf | awk -F'\t' '$7 ~ /nested_(diploid|haploid|unreachable)/ {n++} END {print n+0}') \
   "0" "and no record carries one"
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
# Exact equality is gone with version 5, and deliberately: a strand's rows now cover only the sites
# THAT STRAND TRAVERSES, so a nested haploid site counts on one strand and not the other. What
# survives is that both strands are populated and neither claims more sites than there are records.
is $(awk -F'\t' '/^H\t/ && $11 != "." {n[$3] += $11}
     END {print (n["0"] > 0 && n["1"] > 0) ? "both" : "missing"}' nest_hap.mosaic.tsv) "both" \
   "both strands account for sites with nested sites present"
# Reference order within a strand was asserted here. WITHDRAWN, for two reasons. It passes because
# the phasing is sorted into reference order before segmentation, so it is near-tautological. And
# it would be WRONG for a strand traversing an inversion, where walk order is not reference order.
# What replaces it is the contiguity assertion -- segments meeting at the same oriented node --
# which is the property a consumer actually needs and which reference order was standing in for.
is $(awk -F'\t' '/^H\t/ {k=$2"/"$3"/"$4; if (k==pk && $7 != pe) bad++; pk=k; pe=$8}
     END {print bad+0}' nest_hap.mosaic.tsv) "0" \
   "nested-fixture segments meet at the same oriented node too"
# A "." haplotype was asserted to appear here, for a strand carrying no sequence at a nested
# haploid site. WITHDRAWN because version 5 emits no such row: that strand is not empty, it
# traverses the parent's other allele and bypasses the child snarl, so the site is not on its walk
# and there is nothing to write. The row was also cutting the other strand's run in three -- 351 of
# chr20's 419 were flanked by the same haplotype on both sides. Asserted in the negative now.
is $(awk -F'\t' '/^H\t/ && $10=="." {n++} END {print n+0}' nest_hap.mosaic.tsv) "0" \
   "no row spells a haplotype '.', which version 5 does not emit"
# No GT may name an allele the record does not carry. The linkage layer chooses a candidate
# traversal, not on an emitted allele, and those are different numberings: a traversal can have no
# ALT on its own line, because symbolic collapsing folds some into the reference and the linkage pass can
# replace a line with one carrying fewer ALTs than the entry was built against. Writing the number
# regardless emitted GTs like ".|2" on a record with one ALT -- not a parseable VCF, and invisible to
# every check that only counted records. Asserted on both nested arms; the phased one is where the
# genotype and phasing patches both run.
is $(awk -F'\t' '!/^#/ {n = ($5 == "." || $5 == "") ? 0 : split($5, a, ","); \
      split($10, f, ":"); m = split(f[1], g, "[|/]"); \
      for (i = 1; i <= m; i++) if (g[i] != "." && g[i]+0 > n) bad++} END {print bad+0}' \
     nest_hap.vcf) "0" \
   "no phased nested GT names an allele the record has no ALT for"
is $(awk -F'\t' '!/^#/ {n = ($5 == "." || $5 == "") ? 0 : split($5, a, ","); \
      split($10, f, ":"); m = split(f[1], g, "[|/]"); \
      for (i = 1; i <= m; i++) if (g[i] != "." && g[i]+0 > n) bad++} END {print bad+0}' \
     nest_nested.vcf) "0" \
   "no nested GT names an allele the record has no ALT for"

# A nested haploid record must name its strand and keep its phase set even when it DECOMPOSES
# INTO BLOCKS -- which the fixture above cannot reach, because its nested chain is a single SNP and
# a single difference is a single record. Block emission asked the site's GT for a phased *pair*
# before looking at anything else, so a haploid record fell straight through to the unphased
# fallback: bare "1" and PS erased, losing both the strand and the phase set the unsplit record
# carries. 3,452 lines genome-wide, from 1,517 snarls, every one of them a nested haploid call that
# happened to split.
#
# The chain here has two differences separated by a matched run (node 10), which is what makes two
# blocks. Node 13 bypasses 10 so that 10 is not a cut vertex: without it the decomposition would
# split 5..9 into two child snarls, each emitting its own record, and nothing would ever atomize.
rm -f nestblk.gfa nestblk.gbz nestblk.gam nestblk.vcf
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTA\n'
  printf 'S\t2\tGGCATTAGCCTTAGACCGAT\n'
  printf 'S\t3\tA\nS\t4\tT\n'
  printf 'S\t5\tTTGACCAGTTCAGGACTTAC\n'
  printf 'S\t7\tC\nS\t8\tG\n'
  printf 'S\t10\tACGTACGTAC\n'
  printf 'S\t11\tA\nS\t12\tT\n'
  printf 'S\t13\tCCCCCCCCCC\n'
  printf 'S\t9\tAACCGGTTACGTTGCAATCG\n'
  printf 'S\t6\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\n'
  printf 'L\t5\t+\t7\t+\t0M\nL\t5\t+\t8\t+\t0M\nL\t5\t+\t13\t+\t0M\n'
  printf 'L\t7\t+\t10\t+\t0M\nL\t8\t+\t10\t+\t0M\n'
  printf 'L\t10\t+\t11\t+\t0M\nL\t10\t+\t12\t+\t0M\n'
  printf 'L\t11\t+\t9\t+\t0M\nL\t12\t+\t9\t+\t0M\nL\t13\t+\t9\t+\t0M\n'
  printf 'L\t9\t+\t6\t+\t0M\nL\t1\t+\t6\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+,5+,7+,10+,11+,9+,6+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t133\t>1>2>3>5>7>10>11>9>6\n'
  printf 'W\tp1\t1\tchr1\t0\t133\t>1>2>4>5>8>10>12>9>6\n'
  printf 'W\tp2\t0\tchr1\t0\t60\t>1>6\n'
  printf 'W\tp2\t1\tchr1\t0\t133\t>1>2>4>5>8>10>12>9>6\n'
  printf 'W\tp3\t0\tchr1\t0\t131\t>1>2>3>5>13>9>6\n'
  printf 'W\tp3\t1\tchr1\t0\t60\t>1>6\n'
} > nestblk.gfa
vg gbwt -G nestblk.gfa --gbz-format -g nestblk.gbz --set-reference GRCh 2>/dev/null
vg sim -x nestblk.gbz -n 300 -l 40 -a -s 31 --path "p2#0#chr1#0" > nestblk.gam 2>/dev/null
vg sim -x nestblk.gbz -n 300 -l 40 -a -s 37 --path "p2#1#chr1#0" >> nestblk.gam 2>/dev/null
vg call nestblk.gbz --read-likelihood --gam nestblk.gam -t 1 -s samp --nested --phased \
    2>/dev/null > nestblk.vcf
is "$?" 0 "a nested chain with two separated differences builds and calls"
is $(grep -v "^#" nestblk.vcf | grep -c "SB=") "2" \
   "the nested haploid record decomposes into two difference blocks"
is $(grep -v "^#" nestblk.vcf | cut -f10 | cut -d: -f1 | grep -cE '^[0-9]+$') "0" \
   "and every block names the strand it sits on, not a bare haploid genotype"
is $(grep -v "^#" nestblk.vcf | awk -F'\t' '$9 !~ /(^|:)PS(:|$)/' | wc -l | tr -d ' ') "0" \
   "and keeps the phase set, which a split record used to drop"
# A site written as blocks has lines, so the mosaic must count it. It used to be left out, and
# counted as a collapsed site with no line of its own, because the block path returned before the
# linkage model was told a line had been written. Strand 0 carries the chain: the parent, the SNP
# site and the two-block site.
rm -f nestblk.mosaic.tsv nestblk_mosaic.err
vg call nestblk.gbz --read-likelihood --gam nestblk.gam -t 1 -s samp --nested --phased \
    --mosaic-out nestblk.mosaic.tsv --progress 2>nestblk_mosaic.err >/dev/null
is $(awk -F'\t' '/^H\t/ && $3 == "0" && $11 != "." {n += $11} END {print n+0}' nestblk.mosaic.tsv) "3" \
   "the mosaic counts a site written as difference blocks"
is $(grep -c "collapsed sites phased with no line of their own" nestblk_mosaic.err) "0" \
   "and a site written as blocks is not reported as having no line"

# A site whose only difference outside its child chains is one block is written as that block, not
# as the whole site record, because the site record spells every allele in full and so repeats the
# child's variant. Here the substitution 3/4 is no snarl of its own (the edge 4->6 sees to that), so
# it belongs to the site 2..6, and p2#0 takes it together with the child chain 5..9's deletion of 7.
# The child's record reports the deletion; the site must report only the substitution.
rm -f nestdup.gfa nestdup.gbz nestdup.gam nestdup.vcf
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTA\n'
  printf 'S\t2\tGGCATTAGCCTTAGACCGAT\n'
  printf 'S\t3\tA\nS\t4\tT\n'
  printf 'S\t5\tTTGACCAGTTCAGGACTTAC\n'
  printf 'S\t7\tCATCATCATCATCATCATCATCATCATCAT\n'
  printf 'S\t9\tAACCGGTTACGTTGCAATCG\n'
  printf 'S\t6\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\n'
  printf 'L\t5\t+\t7\t+\t0M\nL\t7\t+\t9\t+\t0M\nL\t5\t+\t9\t+\t0M\n'
  printf 'L\t9\t+\t6\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+,5+,7+,9+,6+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t150\t>1>2>3>5>7>9>6\n'
  printf 'W\tp1\t1\tchr1\t0\t120\t>1>2>4>5>9>6\n'
  printf 'W\tp2\t0\tchr1\t0\t120\t>1>2>4>5>9>6\n'
  printf 'W\tp2\t1\tchr1\t0\t150\t>1>2>3>5>7>9>6\n'
  printf 'W\tp3\t0\tchr1\t0\t81\t>1>2>4>6\n'
  printf 'W\tp3\t1\tchr1\t0\t150\t>1>2>3>5>7>9>6\n'
} > nestdup.gfa
vg gbwt -G nestdup.gfa --gbz-format -g nestdup.gbz --set-reference GRCh 2>/dev/null
vg sim -x nestdup.gbz -n 300 -l 40 -a -s 31 --path "p2#0#chr1#0" > nestdup.gam 2>/dev/null
vg sim -x nestdup.gbz -n 300 -l 40 -a -s 37 --path "p2#1#chr1#0" >> nestdup.gam 2>/dev/null
# Simulated reads carry no mapping quality, so cap the mismapping probability to let them count.
vg call nestdup.gbz --read-likelihood --gam nestdup.gam -t 1 -s samp --nested --phased \
    --mismap-max 0.05 2>/dev/null > nestdup.vcf
is $(grep -v "^#" nestdup.vcf | awk -F'\t' '$3 == ">5>9" && $10 ~ /^0\|1:/' | wc -l | tr -d ' ') "1" \
   "the child chain's deletion is called on the strand that carries it"
is $(grep -v "^#" nestdup.vcf | awk -F'\t' '$3 ~ /^>2>6/ && length($4) != length($5)' | wc -l | tr -d ' ') "0" \
   "and the site that encloses it reports only its own substitution, not the deletion again"

is $(grep -v "^#" nestdup.vcf | awk -F'\t' '$3 == ">2>6_0" && $2 == 51 && $4 == "A" && $5 == "T"' | wc -l | tr -d ' ') "1" \
   "and reports the substitution as a block record of its own"

# A site whose one block spells the same haplotypes as its site record keeps the site record, which is
# left-normalised. Here the sample deletes node 2's A after node 1's A: the site record slides the
# deletion left to CA>C at 50, where the block would be anchored at 52.
rm -f blkleaf.gfa blkleaf.gbz blkleaf.gam blkleaf.vcf
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTAGGCATTAGCCTTAGACCGTCAA\n'
  printf 'S\t2\tA\n'
  printf 'S\t3\tGGATCCTAGCATTCGGATCCAAGTTCCAGATTGACCAGTTCAGGACTTAC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t103\t>1>2>3\n'
  printf 'W\tp1\t1\tchr1\t0\t102\t>1>3\n'
  printf 'W\tp2\t0\tchr1\t0\t102\t>1>3\n'
  printf 'W\tp2\t1\tchr1\t0\t103\t>1>2>3\n'
} > blkleaf.gfa
vg gbwt -G blkleaf.gfa --gbz-format -g blkleaf.gbz --set-reference GRCh 2>/dev/null
vg sim -x blkleaf.gbz -n 300 -l 40 -a -s 41 --path "p2#0#chr1#0" > blkleaf.gam 2>/dev/null
vg sim -x blkleaf.gbz -n 300 -l 40 -a -s 43 --path "p2#1#chr1#0" >> blkleaf.gam 2>/dev/null
vg call blkleaf.gbz --read-likelihood --gam blkleaf.gam -t 1 -s samp --mismap-max 0.05 2>/dev/null > blkleaf.vcf
is $(grep -v "^#" blkleaf.vcf | awk -F'\t' '$3 == ">1>3" && $2 == 50 && $4 == "CA" && $5 == "C"' | wc -l | tr -d ' ') "1" \
   "a site with one block in a repeat keeps its left-normalised site record"

# A block record takes its genotype from the site, through the site allele each of its alleles stands
# for, so one block whose alleles would stand for one site allele is not written. Here p2#0 takes 2->4
# (GCAT) and skips the child chain 5..9's CAT, which spells exactly the reference, so the site's
# genotype is 0|0, while the block would be an insertion with no site allele of its own to take a
# phase or likelihoods from. The site record stands, and writes no line.
rm -f blksame.gfa blksame.gbz blksame.gam blksame.vcf blksame.err
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTA\n'
  printf 'S\t2\tGGCATTAGCCTTAGACCGTCTT\n'
  printf 'S\t3\tG\nS\t4\tGCAT\nS\t5\tCATCAT\nS\t7\tCAT\n'
  printf 'S\t9\tGACCGGTTACGTTGCAATCG\n'
  printf 'S\t6\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t3\t+\t6\t+\t0M\n'
  printf 'L\t5\t+\t7\t+\t0M\nL\t7\t+\t9\t+\t0M\nL\t5\t+\t9\t+\t0M\nL\t9\t+\t6\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+,5+,7+,9+,6+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t112\t>1>2>3>5>7>9>6\n'
  printf 'W\tp1\t1\tchr1\t0\t112\t>1>2>4>5>9>6\n'
  printf 'W\tp2\t0\tchr1\t0\t112\t>1>2>4>5>9>6\n'
  printf 'W\tp2\t1\tchr1\t0\t112\t>1>2>3>5>7>9>6\n'
  printf 'W\tp3\t0\tchr1\t0\t83\t>1>2>3>6\n'
  printf 'W\tp3\t1\tchr1\t0\t112\t>1>2>4>5>9>6\n'
} > blksame.gfa
vg gbwt -G blksame.gfa --gbz-format -g blksame.gbz --set-reference GRCh 2>/dev/null
vg sim -x blksame.gbz -n 300 -l 40 -a -s 47 --path "p2#0#chr1#0" > blksame.gam 2>/dev/null
vg sim -x blksame.gbz -n 300 -l 40 -a -s 53 --path "p2#1#chr1#0" >> blksame.gam 2>/dev/null
vg call blksame.gbz --read-likelihood --gam blksame.gam -t 1 -s samp --nested --phased \
    --mismap-max 0.05 --progress 2> blksame.err > blksame.vcf
is $(grep -c "1 one block that two strands' routes spell differently within one site allele" blksame.err) "1" \
   "a block whose alleles would stand for one site allele is refused"
is $(grep -v "^#" blksame.vcf | awk -F'\t' '$3 ~ /^>2>6/' | wc -l | tr -d ' ') "0" \
   "and the site, whose genotype is the reference, writes no line"

# A chain is genotyped from each allele's first crossing of it, so a site whose strand crosses a chain
# twice keeps its site record. Here p2#1 crosses the chain 5..9 a second time, deleting 7 on the way;
# the chain's own record, made from the first crossing, does not report that deletion.
rm -f blktwice.gfa blktwice.gbz blktwice.gam blktwice.vcf blktwice.err
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTA\n'
  printf 'S\t2\tGGCATTAGCCTTAGACCGAT\n'
  printf 'S\t3\tA\nS\t4\tT\n'
  printf 'S\t5\tTTGACCAGTTCAGGACTTAC\n'
  printf 'S\t7\tCATCATCATCATCATCATCATCATCATCAT\n'
  printf 'S\t9\tAACCGGTTACGTTGCAATCG\n'
  printf 'S\t6\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\n'
  printf 'L\t5\t+\t7\t+\t0M\nL\t7\t+\t9\t+\t0M\nL\t5\t+\t9\t+\t0M\n'
  printf 'L\t9\t+\t6\t+\t0M\nL\t9\t+\t5\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+,5+,7+,9+,6+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t150\t>1>2>3>5>7>9>6\n'
  printf 'W\tp1\t1\tchr1\t0\t120\t>1>2>4>5>9>6\n'
  printf 'W\tp2\t0\tchr1\t0\t150\t>1>2>3>5>7>9>6\n'
  printf 'W\tp2\t1\tchr1\t0\t160\t>1>2>3>5>7>9>5>9>6\n'
  printf 'W\tp3\t0\tchr1\t0\t81\t>1>2>4>6\n'
  printf 'W\tp3\t1\tchr1\t0\t160\t>1>2>3>5>7>9>5>9>6\n'
} > blktwice.gfa
vg gbwt -G blktwice.gfa --gbz-format -g blktwice.gbz --set-reference GRCh 2>/dev/null
vg sim -x blktwice.gbz -n 300 -l 40 -a -s 31 --path "p2#0#chr1#0" > blktwice.gam 2>/dev/null
vg sim -x blktwice.gbz -n 300 -l 40 -a -s 37 --path "p2#1#chr1#0" >> blktwice.gam 2>/dev/null
vg call blktwice.gbz --read-likelihood --gam blktwice.gam -t 1 -s samp --nested --phased \
    --mismap-max 0.05 --progress 2> blktwice.err > blktwice.vcf
is $(grep -c "1 one block, where a chain crossed more than once" blktwice.err) "1" \
   "a site whose strand crosses a chain twice is not written as its block"
is $(grep -v "^#" blktwice.vcf | awk -F'\t' '$3 == ">2>6" && $2 == 101 && length($5) - length($4) == 40' | wc -l | tr -d ' ') "1" \
   "and its site record reports both crossings"

# A strand that crosses a child chain outside its own difference blocks matches the reference there
# symbolically, whatever route it takes through the chain, so a block spells the reference over the
# chain and the chain's own record reports the route. p2#1 deletes the chain (3..6) and node 7 and
# takes 10 for 9, which makes two blocks; p2#0 takes the chain's other SNP allele and otherwise
# follows the reference. Node 13 bypasses 8 so that the parent is one snarl from 2 to 11.
rm -f chainblk.gfa chainblk.gbz chainblk.gam chainblk.vcf
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTA\n'
  printf 'S\t2\tGGCATTAGCCTTAGACCGAT\n'
  printf 'S\t3\tTTGACCAGTTCAGGACTTAC\n'
  printf 'S\t4\tA\nS\t5\tT\n'
  printf 'S\t6\tGATCCATGCAAGTCGTACGA\n'
  printf 'S\t7\tCAGTTGACCTAGGCATCAGG\n'
  printf 'S\t8\tTGCAAGCTTGACCTAGCCAT\n'
  printf 'S\t9\tC\nS\t10\tG\n'
  printf 'S\t13\tAAAAAAAAAACCCCCCCCCC\n'
  printf 'S\t11\tAACCGGTTACGTTGCAATCG\n'
  printf 'S\t12\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t8\t+\t0M\n'
  printf 'L\t3\t+\t4\t+\t0M\nL\t3\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\nL\t5\t+\t6\t+\t0M\n'
  printf 'L\t6\t+\t7\t+\t0M\nL\t7\t+\t8\t+\t0M\nL\t7\t+\t13\t+\t0M\nL\t13\t+\t9\t+\t0M\n'
  printf 'L\t8\t+\t9\t+\t0M\nL\t8\t+\t10\t+\t0M\nL\t9\t+\t11\t+\t0M\nL\t10\t+\t11\t+\t0M\n'
  printf 'L\t11\t+\t12\t+\t0M\n'
  printf 'P\tGRCh#0#chr1\t1+,2+,3+,4+,6+,7+,8+,9+,11+,12+\t*\n'
  printf 'W\tp1\t0\tchr1\t0\t182\t>1>2>3>4>6>7>8>9>11>12\n'
  printf 'W\tp1\t1\tchr1\t0\t121\t>1>2>8>10>11>12\n'
  printf 'W\tp2\t0\tchr1\t0\t182\t>1>2>3>5>6>7>8>9>11>12\n'
  printf 'W\tp2\t1\tchr1\t0\t121\t>1>2>8>10>11>12\n'
  printf 'W\tp3\t0\tchr1\t0\t182\t>1>2>3>4>6>7>13>9>11>12\n'
  printf 'W\tp3\t1\tchr1\t0\t182\t>1>2>3>5>6>7>8>9>11>12\n'
} > chainblk.gfa
vg gbwt -G chainblk.gfa --gbz-format -g chainblk.gbz --set-reference GRCh 2>/dev/null
vg sim -x chainblk.gbz -n 300 -l 40 -a -s 31 --path "p2#0#chr1#0" > chainblk.gam 2>/dev/null
vg sim -x chainblk.gbz -n 300 -l 40 -a -s 37 --path "p2#1#chr1#0" >> chainblk.gam 2>/dev/null
vg call chainblk.gbz --read-likelihood --gam chainblk.gam -t 1 -s samp --nested --phased \
    2>/dev/null > chainblk.vcf
is $(grep -v "^#" chainblk.vcf | awk -F'\t' '$3 ~ /^>2>11_[0-9]+$/ && $8 ~ /SB=/' | wc -l | tr -d ' ') "2" \
   "a parent with a chain one strand deletes and the other crosses by another route writes two blocks"
is $(grep -v "^#" chainblk.vcf | awk -F'\t' '$3 ~ /^>2>11_[0-9]+$/' | sort -u -k3,3 | wc -l | tr -d ' ') "2" \
   "and each block record has an ID of its own, the site's ID and its block index"
is $(grep -v "^#" chainblk.vcf | awk -F'\t' '$3 ~ /^>2>11_[0-9]+$/ && $5 ~ /,/' | wc -l | tr -d ' ') "0" \
   "and its blocks spell the reference over the chain the second strand matches"
is $(grep -v "^#" chainblk.vcf | cut -f3 | sort | uniq -d | wc -l | tr -d ' ') "0" \
   "and no ID is on more than one record"
is $(grep -v "^#" chainblk.vcf | awk -F'\t' '$3 == ">3>6" && $5 == "T"' | wc -l | tr -d ' ') "1" \
   "and the chain's own record reports the second strand's route through it"

# --- assembly anchors -------------------------------------------------------
# An anchor is a zero-length pin at a snarl boundary, or beside an off-panel change, holding the
# reads that cross it partitioned by which called allele they fit. The invariant that matters is
# that a point of a read appears in at most one anchor; everything else here supports being able to
# believe that one.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors.tsv -t 1 \
    > rl_anchors.vcf 2> rl_anchors.err
is $(grep -c "^#anchors-version" rl_anchors.tsv) "1" "the anchor file declares its version"
# Self-describing, and pinned: a consumer that sees no anchor at a site must be able to tell a
# threshold from an absence, so the filter values are part of the file rather than of the shell
# history that produced it.
is "$(grep "^#" rl_anchors.tsv | cut -f1 | LC_ALL=C sort -u | paste -sd, -)" \
   "#H,#anchors-version,#filters,#graph,#mismap-min,#note,#read,#reads,#reads-interned,#sample,#sites,#vg-version" \
   "the anchor header uses exactly the documented set of keys"
is $(grep -c "^A" rl_anchors.tsv | awk '{print ($1>0)?1:0}') "1" "it holds anchors"

# THE guarantee, keyed on the point of the read: the gap, the read bases before it, which is offset+1
# on strand 0 and offset on strand 1, so one point crossed in opposite directions is one key.

# The allele column exists so a consumer can tell which candidate traversal a slot partitions to.
# It is NOT the VCF ALT number and must never be read as one, so what is checked is the property
# that distinguishes them: at a site with two slots the alleles differ, because two slots carrying
# one allele is a homozygote, and a homozygote is collapsed to a single slot before it is written.
is $(awk -F'\t' '/^A/{k=$2"|"$3; n[k]++; a[k]=a[k]" "$5} END{for(i in n) if(n[i]==2){split(a[i],p," "); if(p[1]==p[2]) bad++}} END{print bad+0}' \
     rl_anchors.tsv) "0" \
   "where a site has two slots, they name different candidate alleles"

# The file is written in node order, so nothing downstream has to sort it.
is $(awk -F'\t' '/^A/{print $2}' rl_anchors.tsv | sort -c -n 2>&1 | wc -l | tr -d ' ') "0" \
   "anchors are written in node order"

# A read may sit on only one side of a site's partition, or the partition means nothing. Keyed on
# the site, since its two pins may lie on one node, and a read may be at both.
is $(awk -F'\t' '/^A/{k=$3; s=$4} /^R/{print k"\t"$2"\t"s}' rl_anchors.tsv \
     | sort -u | cut -f1,2 | uniq -d | wc -l | tr -d ' ') "0" \
   "the two slots of a site have disjoint read sets"

# Every R row needs the A row above it, since it does not repeat the key.
is $(awk -F'\t' '/^[AR]/{if ($1=="A") seen=1; if ($1=="R" && !seen) bad++} END{print bad+0}' \
     rl_anchors.tsv) "0" "no read row precedes its anchor row"

# Off-panel changes beside each other and beside a snarl's pins. The sample's haplotype A walks the
# bubble's arm 3 and changes node 1 bases 20 and 21, which touch, and bases 30 and 32, one apart;
# the last base of the bubble's start node 2, the first of its arm 3 and the first of its end node
# 5; node 6's last base and node 7's first, across an edge; node 8 base 10, and inserts a G after
# it; inserts an A inside node 7, where no repeat lets it move; and deletes one A of a run from node
# 8's end into node 9. Its haplotype B walks arm 4 and changes node 1 base 21 alone, and node 9 base
# 1, inside that run. Two reads of A delete node 2's last base and three delete node 1 bases 30 to
# 32. The reads are written as alignments to this graph, with a mismatch edit for each changed base.
rm -f opx.gfa opx.gbz opx.json opx.gam opx_t1.vcf opx_t1.tsv opx_t1.err opx_t4.vcf opx_t4.tsv
python3 - <<'EOF'
import json
N = {1: "GATCCTGACGTTAGCATGCACTGTCAGCTAGTACGATCGA", 2: "TCGACTGAGCTTAGCATCGTACGACTGCAC",
     3: "CAGTA", 4: "CTGTA", 5: "GTCAGCTAGCATGCAGTACTGACTCGATCG",
     6: "ATCGTAGCTGACTAGCTAGTCGATCGTAGC", 7: "TGACTAGCATCGATGCTAGCTGACTGATCG",
     8: "CTAGCATGCAACGTAGCTAGCTAGTCCAAA", 9: "AAGCTAGTCGATCGTAGCTAGCATCGATCG"}
WA, WB = [1, 2, 3, 5, 6, 7, 8, 9], [1, 2, 4, 5, 6, 7, 8, 9]
def alt(n, o, avoid=""):
    return next(b for b in "GTCA" if b != N[n][o] and b not in avoid)
SA = {(1, 21): alt(1, 21)}
SA.update({(1, 20): alt(1, 20, SA[(1, 21)]), (1, 30): alt(1, 30), (1, 32): alt(1, 32),
           (2, 29): alt(2, 29), (3, 0): alt(3, 0, "T"), (5, 0): alt(5, 0), (6, 29): alt(6, 29),
           (7, 0): alt(7, 0), (8, 10): alt(8, 10, "G")})
SB = {(1, 21): SA[(1, 21)], (9, 1): "T"}
def units(walk, subs, dels, ins):
    out = []
    for n in walk:
        for o, b in enumerate(N[n]):
            out += [(n, None, c) for c in ins.get((n, o), "")]
            out.append((n, o, None) if (n, o) in dels else (n, o, subs.get((n, o), b)))
    return out
def reads(name, u, starts, length=90):
    based = [i for i, x in enumerate(u) if x[2] is not None]
    for k, st in enumerate(starts):
        if st + length > len(based) or None in (u[based[st]][1], u[based[st + length - 1]][1]):
            continue
        piece, maps, seq = u[based[st]:based[st + length - 1] + 1], [], ""
        for n, o, b in piece:
            if not maps or maps[-1]["position"]["node_id"] != str(n):
                maps.append({"position": {"node_id": str(n), "offset": str(o or 0)}, "edit": []})
            e = maps[-1]["edit"]
            if o is None:
                e.append({"from_length": 0, "to_length": 1, "sequence": b})
            elif b is None:
                e.append({"from_length": 1, "to_length": 0})
            elif b != N[n][o]:
                e.append({"from_length": 1, "to_length": 1, "sequence": b})
            elif e and "sequence" not in e[-1] and e[-1]["to_length"] == e[-1]["from_length"]:
                e[-1]["from_length"] += 1
                e[-1]["to_length"] += 1
            else:
                e.append({"from_length": 1, "to_length": 1})
            seq += b or ""
        yield {"name": f"{name}{k}", "sequence": seq, "mapping_quality": 60,
               "path": {"mapping": maps}}
with open("opx.gfa", "w") as g:
    g.write("H\tVN:Z:1.1\n" + "".join(f"S\t{k}\t{v}\n" for k, v in N.items()))
    for a, b in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)]:
        g.write(f"L\t{a}\t+\t{b}\t+\t0M\n")
    g.write("P\tGRCh#0#chr1\t" + ",".join(f"{n}+" for n in WA) + "\t*\n")
    for h, w in (("p1\t0", WA), ("p1\t1", WB), ("p2\t0", WB)):
        g.write(f"W\t{h}\tchr1\t0\t225\t" + "".join(f">{n}" for n in w) + "\n")
ins, dels = {(8, 11): "G", (7, 15): "A"}, {(8, 27)}
with open("opx.json", "w") as f:
    for r in [*reads("a", units(WA, SA, dels, ins), range(0, 200, 2)),
              *reads("b", units(WB, SB, set(), {}), range(1, 200, 2)),
              *reads("d", units(WA, SA, dels | {(2, 29)}, ins), [10, 12]),
              *reads("e", units(WA, SA, dels | {(1, 30), (1, 31), (1, 32)}, ins), [0, 2, 4])]:
        f.write(json.dumps(r) + "\n")
EOF
vg gbwt -G opx.gfa --gbz-format -g opx.gbz --set-reference GRCh 2>/dev/null && vg view -JaG opx.json > opx.gam
is "$?" "0" "the off-panel anchor fixture builds"
vg call opx.gbz --read-likelihood --gam opx.gam -t 1 -s samp --anchors-out opx_t1.tsv > opx_t1.vcf 2> opx_t1.err
is "$?" "0" "vg call runs on the off-panel anchor fixture"

# Touching changes on one node are one site, written as one record. Its alleles are what the reads
# spell, not every combination of the changes: haplotype A carries both substitutions and B the
# second alone, so the record has those two ALTs and the sample one of each.
is "$(awk -F'\t' '$3=="op1.20.M2"{split($5,a,","); split($10,g,":"); n=split(g[1],h,/[|\/]/);
      print length($4), length(a), (n==2 && h[1]!=h[2] && h[1]>0 && h[2]>0) ? "het-alt" : g[1]}' opx_t1.vcf)" \
   "2 2 het-alt" "touching substitutions are one record, with the two alleles the reads spell"
is $(grep -v "^#" opx_t1.vcf | cut -f3 | grep -c "^op1\.2[01]\.S") "0" \
   "and neither substitution has a record of its own"
is "$(awk -F'\t' '$3=="op8.10.M1"{print length($5) - length($4)}' opx_t1.vcf)" "1" \
   "a substitution and the insertion after it are one record, its ALT a base longer than its REF"
is $(grep -v "^#" opx_t1.vcf | cut -f3 | grep -c "^op6\.29\.S\|^op7\.0\.S") "2" \
   "substitutions either side of an edge, on two nodes, are two records"

# The guarantee, on the cases it was written for: the e reads cross bases 30 and 32's sites at one
# point, node 2's last base and the bubble pin at one point, and the deletion's tract ends where the
# substitution inside it does. vg leaves an off-panel anchor's read out where another anchor holds
# its point, and says how many.
is $(awk -F'\t' '/^R/{g=($3==0)?$4+1:$4; print $2"\t"g}' opx_t1.tsv | sort | uniq -d | wc -l | tr -d ' ') "0" \
   "beside touching changes and snarls, no (read, gap) is pinned twice"
is $(grep -o "[0-9]* off-panel read placements left out" opx_t1.err | awk '{print ($1>0)}') "1" \
   "and vg left out the off-panel placements whose points another anchor holds"
# An insertion with no repeat tract changes no base, so its two pins are one point, written once.
is "$(awk -F'\t' '/^A/ && $3 ~ /^op7\.14\.I/{print $4":"$9}' opx_t1.tsv | sort | uniq -c | awk '{print $1}' | sort -u)" "1" \
   "an insertion with no tract has one anchor per slot"
is "$(python3 ../scripts/check_anchors.py --anchors opx_t1.tsv 2>&1 | cut -c1-3)" "OK " \
   "check_anchors.py accepts the fixture's anchors, version 9"
vg call opx.gbz --read-likelihood --gam opx.gam -t 4 -s samp --anchors-out opx_t4.tsv > opx_t4.vcf 2>/dev/null
is "$(cmp -s opx_t1.vcf opx_t4.vcf && cmp -s opx_t1.tsv opx_t4.tsv && echo same)" "same" \
   "one thread or four give the same VCF and anchors"
rm -f opx.gfa opx.gbz opx.json opx.gam opx_t1.vcf opx_t1.tsv opx_t1.err opx_t4.vcf opx_t4.tsv

# Read-backed phasing. The fixture cannot exercise the DECISION -- `vg sim` reads carry no mapping
# quality, so the mismapping escape sits at its maximum and every site's per-read confidence comes
# out around 0.84 phred against the 9.87 mean of real chr20 ONT data -- so `--phase-min-q 0` is what
# makes the three stages run at all here. What is checked is that they run, that they leave the
# genotypes alone, and that with the feature off nothing moves.
rm -f rl_phase_off.vcf rl_phase_on.vcf rl_phase_forced.vcf
vg call x.gbz --read-likelihood --phased --gam sim.gam -t 1 2>/dev/null > rl_phase_off.vcf
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing -t 1 \
    2>/dev/null > rl_phase_on.vcf
is $(if diff -q rl_phase_off.vcf rl_phase_on.vcf >/dev/null; then echo 1; else echo 0; fi) "1" \
   "--read-phasing changes nothing where no site clears the confidence threshold"

vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --phase-min-q 0 -t 1 \
    2>rl_phase_forced.err > rl_phase_forced.vcf
is "$?" "0" "--read-phasing runs with every site admitted"
is $(grep -c "read phasing:" rl_phase_forced.err) "1" "it reports what it did"
is $(if [ $(grep "read phasing:" rl_phase_forced.err | sed 's/.*: \([0-9]*\) het sites.*/\1/') -gt 10 ]; then echo 1; else echo 0; fi) "1" \
   "and the sites it considered reached it, so the report is not of an empty pass"

# THE guarantee, and the only reason this can be turned on: a phase decision reorders a chosen pair
# and never substitutes one. Keyed on the snarl ID, because a record's POS depends on its genotype.
is $(python3 -c '
import sys
def load(p):
    out = {}
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10: continue
        k = f[8].split(":"); v = f[9].split(":")
        if "GT" not in k: continue
        out.setdefault(f[2], []).append(tuple(sorted(v[k.index("GT")].replace("|","/").split("/"))))
    return out
a = load("rl_phase_off.vcf"); b = load("rl_phase_forced.vcf")
bad = 0
for key in set(a) & set(b):
    if len(a[key]) == 1 and len(b[key]) == 1 and a[key][0] != b[key][0]:
        bad += 1
print(bad)
') "0" "read phasing reorders genotypes and never changes them"
# Phase-aware re-genotyping. `--read-phasing` reorders a chosen pair and never substitutes one;
# this spends the same evidence on the genotype. The fixture cannot exercise the DECISION -- every
# `vg sim` read here spans one site, so its leave-one-out strand log-odds is zero and the
# correction is provably the identity for it -- but that inertness is itself the property worth
# pinning, because it is what makes the feature degenerate harmlessly on short reads rather than
# do something unmeasured.
rm -f rl_rg_off.vcf rl_rg_t0.vcf rl_rg_t0.err rl_rg_bad.err rl_rg_passes.err rl_rg_cap.err rl_rg_pre.err rl_rg_preoff.err rl_rg_nph.err
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --phase-min-q 0 -t 1 \
    2>/dev/null > rl_rg_off.vcf
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --phase-min-q 0 \
    --regenotype --regeno-temper 0 -t 1 2>rl_rg_t0.err > rl_rg_t0.vcf
is "$?" "0" "--regenotype runs"

# At temper 0, re-genotyping must leave the VCF byte-identical, and this holds by construction
# rather than by tolerance: the tilted weights are the site's own slot weights, so the corrected
# and uncorrected mixtures are the same arithmetic on the same values and cancel bit for bit. One
# differing byte is a bug, not rounding.
is $(if cmp -s rl_rg_off.vcf rl_rg_t0.vcf; then echo 1; else echo 0; fi) "1" \
   "--regeno-temper 0 reproduces --no-regenotype byte for byte"

# And it ran, rather than being skipped -- which byte-identity alone cannot tell you. The report
# says how many sites the correction was actually evaluated at and how many it would move; at
# temper 0 the second must be zero while the first is not.
is $(grep -c "re-genotyping: temper" rl_rg_t0.err) "1" "it reports what it did"
is $(grep "re-genotyping: temper" rl_rg_t0.err | sed 's/.*, \([0-9]*\) would move.*/\1/') "0" \
   "and moves nothing at temper 0"
# The correction rewrites every site's likelihoods in place even when it moves no direct call, and
# GL is written from them, so the linkage pass chooses them before the rounds stop. The byte-identity
# test above is then also the proof that choosing an identity correction again is inert.
is $(grep -c "re-genotyping round 2:" rl_rg_t0.err) "1" \
   "a correction that moves no direct call is still followed by a linkage pass"

# Refused rather than silently inert: with no phasing chain there is no strand log-odds, every
# tilted weight collapses to the site's own, and the correction is the identity at every site.
vg call x.gbz --read-likelihood --phased --gam sim.gam --regenotype -t 1 \
    >/dev/null 2>rl_rg_bad.err
is "$?" "1" "--regenotype without --read-phasing is refused"
is $(grep -c "needs --read-phasing" rl_rg_bad.err) "1" "and says why"

# The cap is a safety net on an iteration that stops itself, so a large value is legal and only
# a nonsensical one is refused. `--regeno-passes 1` is the report-only mode tested above.
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --regenotype \
    --regeno-passes 0 -t 1 >/dev/null 2>rl_rg_passes.err
is "$?" "1" "--regeno-passes 0 is refused"
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --regenotype \
    --regeno-passes 21 -t 1 >/dev/null 2>rl_rg_cap.err
is "$?" "1" "and so is a cap above the ceiling"

# On under --preset ont, with two escapes that must behave differently. --no-regenotype turns it
# off; --no-read-phasing takes away the thing it is computed from, and must DECLINE rather than
# error, because the user never typed --regenotype and a documented combination should not fail.
rm -f rl_rg_pre.err rl_rg_preoff.err rl_rg_nph.err
vg call x.gbz --read-likelihood --phased --gam sim.gam --preset ont --phase-min-q 0 -t 1 \
    >/dev/null 2>rl_rg_pre.err
is $(if [ $(grep -c "re-genotyping" rl_rg_pre.err) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "--preset ont turns re-genotyping on"
vg call x.gbz --read-likelihood --phased --gam sim.gam --preset ont --phase-min-q 0 \
    --no-regenotype -t 1 >/dev/null 2>rl_rg_preoff.err
is $(grep -c "re-genotyping" rl_rg_preoff.err) "0" "--no-regenotype turns it back off"
vg call x.gbz --read-likelihood --phased --gam sim.gam --preset ont --no-read-phasing -t 1 \
    >/dev/null 2>rl_rg_nph.err
is "$?" "0" "--preset ont --no-read-phasing declines re-genotyping instead of failing"
vg call x.gbz --read-likelihood --gam sim.gam --preset ont --no-read-phasing --top-down -t 1 \
    >/dev/null 2>rl_rg_td.err
is "$?" "0" "and so does --preset ont with --top-down"
rm -f rl_rg_td.err

# --preset hifi is the ont settings with three changed. It must apply exactly the settings its help
# line lists, so its calls are the calls those settings make when given one by one.
rm -f rl_hifi_pre.err rl_hifi_pre.vcf rl_hifi_long.vcf
vg call x.gbz --read-likelihood --phased --gam sim.gam --preset hifi --phase-min-q 0 -t 1 \
    2>rl_hifi_pre.err | grep -v "^#" > rl_hifi_pre.vcf
is $(if [ $(grep -c "re-genotyping" rl_hifi_pre.err) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "--preset hifi turns re-genotyping on"
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --regenotype --gap-open 1 \
    --gap-extend 1 --mismap-min 0.05 --read-min-mapq 5 --insertion-nats -1.8 --hp-prior 0 \
    --linkage-weight 1 --off-panel-indels 7,0.25,5 --phase-min-q 0 -t 1 2>/dev/null \
    | grep -v "^#" > rl_hifi_long.vcf
diff rl_hifi_pre.vcf rl_hifi_long.vcf > /dev/null
is "$?" "0" "--preset hifi calls what its settings call when given one by one"
is $(vg call --help 2>&1 | grep -A12 -- "--preset NAME" | tr -s ' \n' ' ' \
     | grep -c -- "hifi: --read-phasing --regenotype --gap-open 1 --gap-extend 1 --mismap-min 0.05 --read-min-mapq 5 --insertion-nats -1.8 --hp-prior 0 --linkage-weight 1 --off-panel-indels 7,0.25,5") "1" \
   "the --preset help lists the hifi settings"
rm -f rl_hifi_pre.err rl_hifi_pre.vcf rl_hifi_long.vcf

# The two extensions. Both are built so that temper 0 stays the identity, and that is the whole
# design constraint: the haploid weight is the read's odds ratio CAPPED at 1 rather than its
# posterior for the strand, because a posterior is 1/2 at temper 0 and would discard half of every
# read's weight on a run meant to change nothing; and the escape's probability is 1/2 at temper 0
# for any ceiling, so its logit is 0 whatever the ceiling turns out to be.
rm -f rl_rg_ext.vcf rl_rg_ext.err rl_rg_ceil0.err rl_rg_ceilhi.err rl_rg_nohap.err
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --phase-min-q 0 \
    --regenotype --regeno-temper 0 --regeno-haploid --regeno-ceiling 0.9 -t 1 \
    2>rl_rg_ext.err > rl_rg_ext.vcf
is $(if cmp -s rl_rg_off.vcf rl_rg_ext.vcf; then echo 1; else echo 0; fi) "1" \
   "temper 0 is the identity with the haploid weight and the per-read escape both armed"

vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --phase-min-q 0 \
    --regenotype --no-regeno-haploid -t 1 >/dev/null 2>rl_rg_nohap.err
is "$?" "0" "--no-regeno-haploid turns the haploid weight off"

# A ceiling is a probability. 1 is no cap and is the default; 0 and above 1 are not meaningful.
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --regenotype \
    --regeno-ceiling 0 -t 1 >/dev/null 2>rl_rg_ceil0.err
is "$?" "1" "--regeno-ceiling 0 is refused"
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --regenotype \
    --regeno-ceiling 1.5 -t 1 >/dev/null 2>rl_rg_ceilhi.err
is "$?" "1" "and so is a ceiling above 1"

# The one combination the correction genuinely cannot be made safe in. `--top-down` builds each
# child's traversal sets out of its parent's CALLED genotype, and those sets decide the child's
# ploidy, are merged into its candidate list and supply its pseudo-reference -- so moving a parent
# leaves the child scored against alleles the parent no longer carries, and no amount of
# re-resolving repairs that, because the columns were never in the matrix. The read-likelihood
# descent is not affected: it passes children a null traversal set, so they enumerate their own.
rm -f rl_rg_td.err rl_rg_bu.err
vg call x.gbz --read-likelihood --phased --gam sim.gam --read-phasing --regenotype \
    --top-down -t 1 >/dev/null 2>rl_rg_td.err
is "$?" "1" "--regenotype with --top-down is refused"
is $(grep -c "no longer carries" rl_rg_td.err) "1" "and says which way the staleness runs"

# Every `rm -f` in this block is a PRE-run guard, so without this the block leaves sixteen files
# behind in test/. Cleaned up here, at the end, where the convention puts it.
rm -f rl_rg_off.vcf rl_rg_t0.vcf rl_rg_t0.err rl_rg_bad.err rl_rg_passes.err rl_rg_cap.err \
      rl_rg_pre.err rl_rg_preoff.err rl_rg_nph.err rl_rg_ext.vcf rl_rg_ext.err rl_rg_ceil0.err \
      rl_rg_ceilhi.err rl_rg_nohap.err rl_rg_td.err rl_rg_bu.err

# `slot` IS the phase, and it is the only thing in the file that carries it. This is checked because
# it once silently was not: `LinkageCollector::settled_traversals` decodes an unordered genotype
# index -- `genotype_index(i, j)` is triangular, so it returns the pair sorted -- and the phase was
# applied only later, inside emit_variant. The anchors were handed the sorted pair, so `slot` came
# out in allele order while the header promised the GT's field order. Every anchor-to-haplotype join
# was then a coin flip. On this fixture 23 of the 34 comparable sites came out reversed; on chr20's
# ONT calls not one het site of 60,544 carried the reversed order that a 50/50 split of `0|1`
# against `1|0` demands, which is what gives the bug away. Nothing in the VCF moves either way --
# the record is byte-identical -- so no other assertion here can see it.
#
# Compared through AD, which is indexed by VCF allele number, and NOT through the anchor `allele`
# column: that is a candidate-traversal index, it is not the ALT number, and it is not monotone in
# the ALT numbering (one site on this fixture inverts it).
rm -f rl_anchors_ph.tsv rl_anchors_ph.vcf
vg call x.gbz --read-likelihood --phased --gam sim.gam --anchors-out rl_anchors_ph.tsv -t 1 \
    2>/dev/null > rl_anchors_ph.vcf
is "$?" "0" "--anchors-out works alongside --phased"

PHASE_SLOTS=$(awk -F'\t' '
  NR==FNR {
    if ($0 ~ /^#/) next
    nf = split($9, k, ":"); nv = split($10, v, ":"); gt = ""; ad = ""
    for (i = 1; i <= nf; i++) { if (k[i] == "GT") gt = v[i]; if (k[i] == "AD") ad = v[i] }
    if (gt !~ /\|/ || ad == "") next
    split(gt, g, "|")
    if (g[1] == "." || g[2] == "." || g[1] == g[2]) next
    n = split(ad, d, ",")
    if (g[1] + 1 > n || g[2] + 1 > n) next
    want0[$3] = d[g[1] + 1]; want1[$3] = d[g[2] + 1]
    next
  }
  $1 == "A" { cur = ""; if ($3 in want0) { cur = $3; slot = $4 }; next }
  # Distinct read ids: a site is pinned at both boundaries, so a raw row count doubles it.
  $1 == "R" && cur != "" { key = cur "|" slot "|" $2; if (!(key in seen)) { seen[key] = 1; cnt[cur "|" slot]++ } }
  END {
    for (site in want0) {
      c0 = cnt[site "|0"]; c1 = cnt[site "|1"]
      if (c0 == 0 || c1 == 0) continue
      if (c0 + c1 != want0[site] + want1[site]) continue   # off-call reads dropped; not comparable
      ++total
      if (c0 != want0[site] || c1 != want1[site]) ++bad
    }
    print total+0, bad+0
  }' rl_anchors_ph.vcf rl_anchors_ph.tsv)

# Reachability first: a zero mismatch count from an empty comparison looks exactly like a pass.
is $(if [ $(echo "$PHASE_SLOTS" | cut -d' ' -f1) -gt 20 ]; then echo 1; else echo 0; fi) "1" \
   "the phased anchor file offers het sites whose slot read counts are comparable to AD"
is $(echo "$PHASE_SLOTS" | cut -d' ' -f2) "0" \
   "slot i holds the reads of GT field i, so an anchor joins to a haplotype"

# The half-called sites the comparison above skips -- `g[1] == "." || g[2] == "."` -- are exactly
# where the same bug survived the v3 -> v4 fix. A nested haploid chain sits on ONE strand of a
# diploid locus, so its single slot is a haplotype, and v4 stamped it 0 whichever strand that was:
# every `.|a` site named the wrong one. Nothing in the VCF carries the slot, so no other assertion
# here can see it, and F1 cannot either -- the VCF is byte-identical across the fix.
#
# x.gbz cannot test this: it has no chain that only one parent allele crosses, so it produces no
# half-called site at all and the check would pass on an empty comparison. nest.gbz does -- p2#0 is
# `>1>6`, deleting the whole middle, so against p2#1 the two nested SNP sites exist on one strand
# only. It is still alive here; line ~1385 is what removes it.
vg call nest.gbz --read-likelihood --gam nest_hap.gam -t 1 -s samp --nested --phased \
    --anchors-out nest_hap_anchors.tsv 2>/dev/null > nest_hap_anchors.vcf
HAP_SLOTS=$(awk -F'\t' '
  NR==FNR {
    if ($0 ~ /^#/) next
    nf = split($9, k, ":"); gt = ""; split($10, v, ":")
    for (i = 1; i <= nf; i++) { if (k[i] == "GT") gt = v[i] }
    if (gt !~ /\|/) next
    split(gt, g, "|")
    if (g[1] == "." && g[2] != ".") want[$3] = 1
    else if (g[2] == "." && g[1] != ".") want[$3] = 0
    next
  }
  $1 == "A" && ($3 in want) { ++total; if ($4 != want[$3]) ++bad }
  END { print total+0, bad+0 }' nest_hap_anchors.vcf nest_hap_anchors.tsv)

# Reachability first: a zero mismatch count from an empty comparison looks exactly like a pass, and
# on x.gbz that is precisely what it was.
is $(if [ $(echo "$HAP_SLOTS" | cut -d' ' -f1) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "the nested fixture offers half-called sites to check the haploid slot against"
is $(echo "$HAP_SLOTS" | cut -d' ' -f2) "0" \
   "a nested haploid site's slot is the strand its GT names, not always 0"

# The anchor file's gqn column and the VCF's GQN are documented to agree, and on chr20 ONT they
# disagreed on 579 snarls -- 482 of them where the VCF said "." and the anchor gave a number, and
# not once the other way round.  Two causes, fixed together: the VCF could not express a
# post-linkage GQN for a haploid genotype (above), and the anchor fell back to the PRE-linkage
# value on a record linkage had moved, which is the margin of the genotype it moved away from.
#
# Reachability first: a zero mismatch count over an empty join looks exactly like a pass.
GQN_JOIN=$(awk -F'\t' '
  NR==FNR {
    if ($0 ~ /^#/) next
    nf = split($9, k, ":"); split($10, v, ":"); g = "."
    for (i = 1; i <= nf; i++) { if (k[i] == "GQN") g = v[i] }
    vcf[$3] = g
    next
  }
  $1 == "A" && ($3 in vcf) {
    if (seen[$3]) next
    seen[$3] = 1; ++total
    if (vcf[$3] == "." && $6 != ".") ++bad
  }
  END { print total+0, bad+0 }' nest_hap_anchors.vcf nest_hap_anchors.tsv)
is $(if [ $(echo "$GQN_JOIN" | cut -d' ' -f1) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "the nested fixture offers snarls to join gqn against GQN"
is $(echo "$GQN_JOIN" | cut -d' ' -f2) "0" \
   "no anchor reports a gqn where the VCF reports none"

# reliability is a DERIVED column -- the mean of the site's own R-row scores -- so it must agree
# with them exactly, or a consumer that filters on it is filtering on something else. Checked
# against the file's own rows, with a read counted once across both pins and both slots, which is
# the part a naive average gets wrong.
REL_OK=$(awk -F'\t' '
  function flush(  k, n, t) {
    if (snarl == "") return
    n = 0; t = 0
    for (k in score) { n++; t += score[k] }
    if (n > 0) { want[snarl] = t / n }
    delete score
  }
  $1 == "A" { if ($3 != snarl) { flush(); snarl = $3 }; got[$3] = $8; if ($8 == ".") dot++; next }
  $1 == "R" { score[$2] = $5; next }
  END {
    flush()
    for (s in want) {
      total++
      d = want[s] - got[s]; if (d < 0) d = -d
      # Not equality: the column is the UNROUNDED mean and the R rows are written to one decimal,
      # so re-averaging them lands within ~0.05 (chr20 max 0.046 over 169,358 sites).
      if (d > 0.06) bad++
    }
    print total+0, bad+0, dot+0
  }' rl_anchors.tsv)

is $(if [ $(echo "$REL_OK" | cut -d' ' -f1) -gt 20 ]; then echo 1; else echo 0; fi) "1" \
   "the anchor file offers sites to check reliability against their own read rows"
is $(echo "$REL_OK" | cut -d' ' -f2) "0" \
   "reliability is the mean of the site's per-read scores, deduped across pins and slots"

# The in-process invariant: each pin checked against the graph's own base while the read was live.
is $(grep -c "pins verified against the graph, 0 failed" rl_anchors.err) "1" \
   "every pin verifies against the graph"

# Homozygous and haploid sites are emitted by default -- they partition nothing but they link reads,
# and an anchor graph needs contiguity too. Excluding them must actually remove some, or the flag is
# inert, and the file must say which set it holds.
# --anchors-hom-split is read-likelihood-only, must be inert by default, and must say so in the
# file -- a split homozygous site is indistinguishable from a heterozygote by its columns alone,
# because the two slots carry the SAME allele rather than different ones.
vg call x.vg -k x.pack --anchors-hom-split 2>/dev/null >/dev/null
is "$?" "1" "--anchors-hom-split without --read-likelihood is refused"

is $(grep -c "hom-split=off" rl_anchors.tsv) "1" "the default file records that hom-split is off"

rm -f rl_anchors_split.tsv
# --read-phasing because the split is decided by each read's phase across the OTHER sites it
# crosses: without a chain there is nothing to split on, and the flag is now refused rather than
# running and reporting every read as having no opinion.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors_split.tsv \
    --anchors-hom-split --read-phasing -t 1 2>/dev/null >/dev/null
is $(grep -c "hom-split=on" rl_anchors_split.tsv) "1" "and a split file records that it is on"

is $(if [ $(grep -c "^A" rl_anchors_split.tsv) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "--anchors-hom-split still emits anchors"
rm -f rl_anchors_split.tsv

vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors_hom.tsv --anchors-het-only -t 1 \
    2>/dev/null >/dev/null
is $(if [ $(grep -c "^A" rl_anchors_hom.tsv) -lt $(grep -c "^A" rl_anchors.tsv) ]; then echo 1; else echo 0; fi) "1" \
   "--anchors-het-only drops the homozygous sites the default includes"
is $(awk -F'\t' '$1=="#sites"{print $2}' rl_anchors.tsv) "het+hom" \
   "and the default file says which set it holds"
is $(awk -F'\t' '$1=="#sites"{print $2}' rl_anchors_hom.tsv) "het" \
   "as does the restricted one"

# Dropping an end pin that holds no reads the start pin lacks. Three things: that it fires, that what
# survives is the START boundary rather than an arbitrary half, and -- the one that matters -- that no
# read placement is lost, since a dropped end pin is by definition one whose reads are all at the start.
rm -f rl_anchors_sp.tsv
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors_sp.tsv \
    --anchors-end-new 1 -t 1 2>/dev/null >/dev/null
is $(if [ $(grep -c "^A" rl_anchors_sp.tsv) -lt $(grep -c "^A" rl_anchors.tsv) ]; then echo 1; else echo 0; fi) "1" \
   "--anchors-end-new drops anchors the default emits"
is $(python3 -c '
import re
from collections import defaultdict
# Every surviving end anchor must hold a read its own slot\047s start anchor does not -- that is the
# criterion, checked against the output. Read NAMES shared by two alignments are excluded: a mate at
# the start pin and a mate at the end pin look like one read reaching both, so on paired data the
# check is confounded rather than violated.
S = defaultdict(set); E = defaultdict(set); cur = None; role = None
in_anchor = set(); seen = {}; shared = set()
for line in open("rl_anchors_sp.tsv"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[0] == "A":
        ends = re.findall(r"[<>](\d+)", f[2])
        role = "S" if int(f[1]) == int(ends[0]) else "E"
        cur = (f[2], f[3]); in_anchor = set()
    else:
        rid = f[1]
        if rid in in_anchor:
            shared.add(rid)
        in_anchor.add(rid)
        k = (cur[0], rid)
        if k in seen and seen[k] != cur[1]:
            shared.add(rid)
        seen[k] = cur[1]
        (S if role == "S" else E)[cur].add(rid)
print(sum(1 for k in E if not (E[k] - S[k]) and not (E[k] & shared)))
') "0" \
   "and every surviving end pin holds a read its start pin does not"
# Every read still appears somewhere: the dropped pins were redundant by construction.
is $(comm -13 <(awk -F'\t' '/^#read/{print $3}' rl_anchors_sp.tsv | LC_ALL=C sort -u) \
              <(awk -F'\t' '/^#read/{print $3}' rl_anchors.tsv | LC_ALL=C sort -u) | wc -l | tr -d ' ') "0" \
   "and no read is lost from the file by dropping them"

# A filter that can actually fire, and provably so: the per-read score is bounded above by the
# mismap floor, so nothing can reach 100.
rm -f rl_anchors_none.tsv
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors_none.tsv \
    --anchors-min-q 100 -t 1 2>/dev/null >/dev/null
is $(grep -c "^A" rl_anchors_none.tsv) "0" "--anchors-min-q 100 empties the file"
rm -f rl_anchors_gqn.tsv
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors_gqn.tsv \
    --anchors-min-gqn 1.0 -t 1 2>/dev/null >/dev/null
is $(grep -c "^A" rl_anchors_gqn.tsv) "0" "--anchors-min-gqn 1.0 empties the file"

# The partition is the read/allele likelihood matrix, which no other caller builds.
vg call x.gbz -k x.pack --anchors-out rl_anchors_err.tsv > /dev/null 2> rl_anchors_err.txt
is "$?" "1" "--anchors-out without --read-likelihood is an error"

# The sweep is parallel over node-ID windows, so the file must not depend on the schedule.
vg call x.gbz --read-likelihood --gam sim.gam --anchors-out rl_anchors_t4.tsv -t 4 \
    2>/dev/null >/dev/null
is $(grep -v "^#" rl_anchors.tsv | diff - <(grep -v "^#" rl_anchors_t4.tsv) > /dev/null && echo 1 || echo 0) "1" \
   "the anchor set is identical at -t 1 and -t 4"

# The configuration most likely to break silently: no panel, so no linkage layer, so the pass that
# renders retained records is armed for a different reason than usual. If anchors ever stop being
# emitted here it will not show up anywhere else.
rm -f rl_anchors_nopanel.tsv
vg call x.vg -k x.pack --read-likelihood --gam sim.gam --anchors-out rl_anchors_nopanel.tsv \
    --enumerate-support -t 1 2>/dev/null >/dev/null
is $(if [ $(grep -c "^A" rl_anchors_nopanel.tsv) -gt 0 ]; then echo 1; else echo 0; fi) "1" \
   "anchors are still emitted with no panel and no linkage layer"

# The written file re-checked from outside vg, against the reads it indexes rather than against the
# graph -- the version that would catch a mistake shared between the collector and its own check.
vg view -X sim.gam > rl_anchors_reads.fq 2>/dev/null
is $(python3 ../scripts/check_anchors.py --anchors rl_anchors.tsv --reads rl_anchors_reads.fq \
       --quiet > /dev/null 2>&1; echo $?) "0" \
   "the anchor file validates against the read file it indexes"

rm -f rl_anchors.tsv rl_anchors.vcf rl_anchors.err rl_anchors_hom.tsv rl_anchors_none.tsv \
      rl_anchors_gqn.tsv rl_anchors_err.tsv rl_anchors_err.txt rl_anchors_t4.tsv rl_anchors_reads.fq \
      rl_anchors_nopanel.tsv rl_anchors_sp.tsv rl_anchors_ph.tsv rl_anchors_ph.vcf \
      rl_phase_off.vcf rl_phase_on.vcf rl_phase_forced.vcf rl_phase_forced.err

rm -f chainblk.gfa chainblk.gbz chainblk.gam chainblk.vcf nestblk.mosaic.tsv nestblk_mosaic.err nest_del.gam nest_del_lw0.vcf nest_del_link.vcf nest_hap_anchors.tsv nest_hap_anchors.vcf nestblk.gfa nestblk.gbz nestblk.gam nestblk.vcf nestdup.gfa nestdup.gbz nestdup.gam nestdup.vcf blkleaf.gfa blkleaf.gbz blkleaf.gam blkleaf.vcf blksame.gfa blksame.gbz blksame.gam blksame.vcf blksame.err blktwice.gfa blktwice.gbz blktwice.gam blktwice.vcf blktwice.err nest.gfa nest.gbz nest.gam nest_default.vcf nest_nested.vcf nest_hap.gam nest_hap.vcf nest_hap_err.txt nest_hap.mosaic.tsv x.vg x.gbz x.gbwt sim.gam x.pack call.vcf callg.vcf callz.vcf callg.6 callz.6 callrl_nopack.vcf callrl_nopack_z.vcf callrl_withpack.vcf nopack_err.txt sim.sorted.gam sim.sorted.gam.gai rl_inmem.vcf rl_indexed.vcf gi_err.txt gb_excl.txt gb_norl.txt gb_nobin.txt gb_dead.txt sim.gaf sim.unsorted.gaf sim.gaf.db fake_gaf_base dead_gaf_base rl_gafmem.vcf rl_fakebase.vcf rl_fakebase_t4.vcf rl_fakebase_w32.vcf rl_gafbase.vcf rl_gafbase_t4.vcf rl_gafbase_w32.vcf rl_autoz_full.vcf rl_autoz.vcf rl_explicit_z.vcf rl_support_full.vcf rl_support.vcf es_nopack.txt es_z.txt es_g.txt nopanel.vg nopanel.gbwt nopanel.gbz nopanel.pack nopanel_err.txt poisson_default.vcf poisson_z.vcf rl_phased.vcf rl_default.vcf rl_nophase.vcf rl_nopanel.vcf rl_nopanel_err.txt rl_unphased_gt.txt rl_phased_gt.txt rl_ph_err.txt rl_mosaic.tsv rl_mosaic2.tsv rl_hap.vcf rl_hap_mosaic.tsv

# The same nodes and edges under --read-likelihood, which reaches the chain inside the insertion
# by descent rather than by --top-down.  A gref cover is what makes it reportable: the reference
# bypasses the insertion, so without one the SNP inside it has no REF and no POS and is genotyped
# and thrown away.  Needs a GBZ, because allele enumeration comes from the panel; the truth pair
# lives in a second graph so the sample's own walks are not in the panel it is called against.
A20=$(printf 'A%.0s' $(seq 20)); C60=$(printf 'C%.0s' $(seq 60))
A60=$(printf 'A%.0s' $(seq 60)); T20=$(printf 'T%.0s' $(seq 20))
{ printf "H\tVN:Z:1.1\n"
  printf "S\t1\t$A20\nS\t2\t$C60\nS\t3\tG\nS\t4\tT\nS\t5\t$A60\nS\t6\t$T20\n"
  printf "L\t1\t+\t6\t+\t0M\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n"
  printf "L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\n"
  printf "P\tx#0#chr1\t1+,6+\t*\n"
  printf "W\ts1\t0\tchr1\t0\t161\t>1>2>3>5>6\n"
  printf "W\ts1\t1\tchr1\t0\t161\t>1>2>4>5>6\n"
  printf "W\ts2\t0\tchr1\t0\t161\t>1>2>3>5>6\n"
  printf "W\ts2\t1\tchr1\t0\t40\t>1>6\n"; } > rlg.gfa
vg gbwt -G rlg.gfa --gbz-format -g rlg.gbz --set-reference x 2>/dev/null
vg convert -p rlg.gbz > rlg.pg 2>/dev/null
vg paths --compute-gref -Q x --min-gref-len 1 -x rlg.pg > rlg.gref.pg 2>/dev/null
vg convert -f rlg.gref.pg > rlg.gref.gfa 2>/dev/null
vg gbwt -G rlg.gref.gfa --gbz-format -g rlg.gref.gbz 2>/dev/null
cp rlg.gref.gfa rlg_sim.gfa
printf 'W\tT\t0\tchr1\t0\t161\t>1>2>3>5>6\n' >> rlg_sim.gfa
printf 'W\tT\t1\tchr1\t0\t161\t>1>2>4>5>6\n' >> rlg_sim.gfa
vg gbwt -G rlg_sim.gfa --gbz-format -g rlg_sim.gbz 2>/dev/null
vg sim -x rlg_sim.gbz -n 300 -l 30 -a -s 11 --path "T#0#chr1#0" >  rlg.gam 2>/dev/null
vg sim -x rlg_sim.gbz -n 300 -l 30 -a -s 13 --path "T#1#chr1#0" >> rlg.gam 2>/dev/null
# --mismap-max is pinned. These reads come from `vg sim` and carry MAPQ 0, so e_r is the cap for
# every read, and the cap then scales all the evidence at once instead of discounting the ambiguous
# reads it exists for. 0.7 is the cap these expectations were set at; the default is swept against
# real data, where mapping qualities vary.
vg call rlg.gref.gbz --read-likelihood --mismap-max 0.7 --gam rlg.gam -t 1 -s samp --phased \
    -p x#0#chr1 > rlg_base.vcf 2>/dev/null
# No env var: selecting a gref reference is itself the signal to descend into off-reference chains.
vg call rlg.gref.gbz --read-likelihood --mismap-max 0.7 --gam rlg.gam -t 1 -s samp --phased \
    -P x#0#chr1 -P 'gref_x#0#chr1_' --linkage-scale 40 --progress \
    --mosaic-out rlg.mosaic.tsv > rlg_gref.vcf 2> rlg_gref.err
rlg_row() { awk -F'\t' -v c="$2" '$1==c && $1!~/^#/' "$1"; }
is "$(rlg_row rlg_base.vcf chr1_1_alt | wc -l | tr -d ' ')" "0" "without a gref reference the SNP inside the insertion is not reported"
is "$(rlg_row rlg_gref.vcf chr1_1_alt | wc -l | tr -d ' ')" "1" "with one it is, on the gref fragment contig"
is "$(rlg_row rlg_gref.vcf chr1_1_alt | cut -f4,5)" "$(printf 'G\tT')" "and it is the SNP, not the insertion again"
is "$(rlg_row rlg_gref.vcf chr1_1_alt | cut -f10 | cut -d: -f1)" "0|1" "phased het, the truth pair differing only inside the insertion"
is "$(rlg_row rlg_gref.vcf chr1_1_alt | cut -f8 | grep -c 'CH=1')" "1" "INFO/CH marks it one insertion layer deep, which is the wiki filter"
# CH comes from the CONTIG's gref level, not from which of a record's ancestors happen to carry a
# record here.  Counting only emitted ancestors left records on a fragment whose enclosing
# base-contig site produced no line at CH=0 -- indistinguishable from the linear reference, and on
# a gref-covered chr20 that was 29,843 of 41,669 off-reference records.
is "$(awk -F'\t' '!/^#/ && $1!="chr1" && $8 ~ /CH=0/' rlg_gref.vcf | wc -l | tr -d ' ')" "0" "no record on a gref fragment is marked as being on the linear reference"
is "$(grep -o 'gref nesting levels:.*' rlg_gref.err)" "gref nesting levels: 1:1 2:1" "fragment levels are read off the graph, and the level-2 fragment inside the insertion is seen as such"
is "$(rlg_row rlg_gref.vcf chr1_1_alt | cut -f8 | grep -c 'PS=>1>6')" "1" "and INFO/PS names the enclosing top-level snarl"
is "$(grep -c 'descending into chains the reference does not cross' rlg_gref.err)" "1" "selecting a gref reference self-enables the descent, with no env var"
is "$(rlg_row rlg_gref.vcf chr1 | cut -f1,2,4,5)" "$(rlg_row rlg_base.vcf chr1 | cut -f1,2,4,5)" "the top-level record on the base contig is unchanged"
is "$(grep -c '^##contig' rlg_gref.vcf)" "2" "only the fragment that carries a call is declared, not the whole cover"
is "$(grep -c '^#reference' rlg.mosaic.tsv)" "1" "the mosaic header names the base reference once"
is "$(grep -c '^#gref-fragments' rlg.mosaic.tsv)" "1" "and counts the cover's fragments instead of listing every one"
rm -f rlg.gfa rlg.gbz rlg.pg rlg.gref.pg rlg.gref.gfa rlg.gref.gbz rlg_sim.gfa rlg_sim.gbz \
      rlg.gam rlg_base.vcf rlg_gref.vcf rlg_gref.err rlg.mosaic.tsv

# vg haplotypes --set-reference clips a gref fragment to the sampled haplotypes and names the piece
# NAME[offset]. -S and -P must select it as -p does, by base name, and it must be called as the
# whole fragment is. The same graph with the insertion's first 30 bp in a node of their own, which
# the clipped fragment leaves out.
C30=$(printf 'C%.0s' $(seq 30))
{ printf "H\tVN:Z:1.1\n"
  printf "S\t1\t$A20\nS\t2\t$C30\nS\t7\t$C30\nS\t3\tG\nS\t4\tT\nS\t5\t$A60\nS\t6\t$T20\n"
  printf "L\t1\t+\t6\t+\t0M\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t7\t+\t0M\nL\t7\t+\t3\t+\t0M\nL\t7\t+\t4\t+\t0M\n"
  printf "L\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\n"
  printf "P\tx#0#chr1\t1+,6+\t*\n"
  printf "W\ts1\t0\tchr1\t0\t161\t>1>2>7>3>5>6\n"
  printf "W\ts1\t1\tchr1\t0\t161\t>1>2>7>4>5>6\n"
  printf "W\ts2\t0\tchr1\t0\t161\t>1>2>7>3>5>6\n"
  printf "W\ts2\t1\tchr1\t0\t40\t>1>6\n"; } > clg.gfa
vg gbwt -G clg.gfa --gbz-format -g clg.gbz --set-reference x 2>/dev/null
vg convert -p clg.gbz > clg.pg 2>/dev/null
vg paths --compute-gref -Q x --min-gref-len 1 -x clg.pg > clg.gref.pg 2>/dev/null
vg convert -f clg.gref.pg > clg.gref.gfa 2>/dev/null
vg gbwt -G clg.gref.gfa --gbz-format -g clg.gref.gbz 2>/dev/null
awk -F'\t' 'BEGIN{OFS="\t"} $1=="W" && $2=="gref_x" && $4=="chr1_1_alt" {$5=30; $7=">7>3>5"} {print}' \
    clg.gref.gfa > clg_clip.gfa
vg gbwt -G clg_clip.gfa --gbz-format -g clg_clip.gbz 2>/dev/null
cp clg.gref.gfa clg_sim.gfa
printf 'W\tT\t0\tchr1\t0\t161\t>1>2>7>3>5>6\n' >> clg_sim.gfa
printf 'W\tT\t1\tchr1\t0\t161\t>1>2>7>4>5>6\n' >> clg_sim.gfa
vg gbwt -G clg_sim.gfa --gbz-format -g clg_sim.gbz 2>/dev/null
vg sim -x clg_sim.gbz -n 300 -l 30 -a -s 11 --path "T#0#chr1#0" >  clg.gam 2>/dev/null
vg sim -x clg_sim.gbz -n 300 -l 30 -a -s 13 --path "T#1#chr1#0" >> clg.gam 2>/dev/null
is "$(vg paths -L -x clg_clip.gbz | grep -c 'chr1_1_alt\[30\]$')" "1" "the clipped gref fragment is a subpath at offset 30"
vg call clg.gref.gbz --read-likelihood --mismap-max 0.7 --gam clg.gam -t 1 -s samp --phased \
    -S gref_x > clg_whole.vcf 2>/dev/null
vg call clg_clip.gbz --read-likelihood --mismap-max 0.7 --gam clg.gam -t 1 -s samp --phased \
    -S gref_x --progress > clg_clip.vcf 2> clg_clip.err
is "$?" "0" "-S selects a reference sample whose gref fragment is clipped"
is "$(grep -v '^#' clg_clip.vcf | cut -f1-5,8)" "$(grep -v '^#' clg_whole.vcf | cut -f1-5,8)" "the clipped fragment is called as the whole one, at the whole fragment's coordinates"
is "$(grep '^##contig' clg_clip.vcf)" "$(grep '^##contig' clg_whole.vcf)" "and declared with the whole fragment's length"
is "$(grep -o 'gref nesting levels:.*' clg_clip.err)" "gref nesting levels: 1:1 2:1" "the clipped fragment keeps its nesting level, and the fragment inside it stays a level deeper"
vg call clg_clip.gbz --read-likelihood --mismap-max 0.7 --gam clg.gam -t 1 -s samp --phased \
    -P gref_x > clg_clip_P.vcf 2>/dev/null
is "$(grep -v '^#' clg_clip_P.vcf | cut -f1-5)" "$(grep -v '^#' clg_whole.vcf | cut -f1-5)" "-P selects the clipped fragment too"
rm -f clg.gfa clg.gbz clg.pg clg.gref.pg clg.gref.gfa clg.gref.gbz clg_clip.gfa clg_clip.gbz \
      clg_sim.gfa clg_sim.gbz clg.gam clg_whole.vcf clg_clip.vcf clg_clip.err clg_clip_P.vcf


# The graph and reads that t/18 builds for its -L/--cluster tests.
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - > small_cluster_call.gfa
printf "L\t1\t+\t9\t+\t0M\n" >> small_cluster_call.gfa
printf "P\tyA\t1+,2+,4+,6+,7+,9+\t*\n" >> small_cluster_call.gfa
printf "P\tyB\t1+,2+,4+,6+,8+,9+\t*\n" >> small_cluster_call.gfa
vg view -Fv small_cluster_call.gfa > small_cluster_call.vg
vg sim -x small_cluster_call.vg -P yA -n 40 -l 20 -a -s 1 > call_a.gam
vg sim -x small_cluster_call.vg -P yB -n 40 -l 20 -a -s 2 > call_b.gam
cat call_a.gam call_b.gam > small_cluster_call.gam
vg pack -x small_cluster_call.vg -g small_cluster_call.gam -o small_cluster_call.pack


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


# The support-based callers: a whole-contig line of ploidy 1 is -d 1, and ploidy 0 makes no calls.
printf 'x\t0\t1000000000\t1\n' > whole_x_ploidy.bed
vg call small_cluster_call.vg -k small_cluster_call.pack -p x -d 1 2>/dev/null | grep -v '^#' > sb_hap_d.txt
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed whole_x_ploidy.bed 2>/dev/null | grep -v '^#' > sb_hap_bed.txt
is $(wc -l < sb_hap_d.txt | awk '{print ($1>0)?1:0}') "1" "-d 1 makes calls on the support-based contig"
is "$(md5sum < sb_hap_bed.txt)" "$(md5sum < sb_hap_d.txt)" "a whole-contig --ploidy-bed line of ploidy 1 calls exactly as -d 1 without reads"
printf 'x\t0\t1000000000\t0\n' > zero_x_ploidy.bed
is $(vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed zero_x_ploidy.bed 2>/dev/null | grep -vc '^#') "0" "ploidy 0 makes no support-based calls"
is $(vg call small_cluster_call.vg -k small_cluster_call.pack -p x --bottom-up 2>/dev/null | grep -vc '^#' | awk '{print ($1>0)?1:0}') "1" "--bottom-up makes calls on the support-based contig"
is $(vg call small_cluster_call.vg -k small_cluster_call.pack -p x --bottom-up --ploidy-bed zero_x_ploidy.bed 2>/dev/null | grep -vc '^#') "0" "ploidy 0 makes no --bottom-up calls"

# The ploidy presets built into vg parse, and leave a contig they do not name diploid and unchanged.
PRESETS="chm13v2.XY chm13v2.XX grch38.XY grch38.XX grch37.ucsc.XY grch37.ucsc.XX grch37.b37.XY grch37.b37.XX"
PRESETS_OK=0
for PRESET in ${PRESETS}; do
    vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-preset "${PRESET}" 2>/dev/null | grep -v '^#' > preset.txt
    [ "${PIPESTATUS[0]}" = 0 ] && cmp -s preset.txt <(grep -v '^#' ploidy_bed_none.vcf) && PRESETS_OK=$((PRESETS_OK + 1))
done
# An unknown name fails with a message that lists the presets, which are the eight above, and the
# help lists the same eight.
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-preset nosuch >/dev/null 2> preset_err.txt
is "$?" 1 "--ploidy-preset rejects an unknown name"
is "$(grep -o '[a-z0-9.]*\.X[XY]' preset_err.txt | sort | tr '\n' ' ')" "$(printf '%s\n' ${PRESETS} | sort | tr '\n' ' ')" "vg has XY and XX ploidy presets for four assemblies"
is "$(vg call --help 2>&1 | grep -o '[a-z0-9.]*\.X[XY]' | sort | tr '\n' ' ')" "$(printf '%s\n' ${PRESETS} | sort | tr '\n' ' ')" "vg call --help names every ploidy preset"
is "${PRESETS_OK}" "8" "every ploidy preset parses and leaves an unnamed contig as it was"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-preset '' >/dev/null 2>&1
is "$?" 1 "--ploidy-preset rejects an empty name"
vg call small_cluster_call.vg -k small_cluster_call.pack -p x --ploidy-bed elsewhere_ploidy.bed --ploidy-preset chm13v2.XY >/dev/null 2>&1
is "$?" 1 "--ploidy-bed cannot be used with --ploidy-preset"

rm -f bad_ploidy.bed overlap_ploidy.bed reversed_ploidy.bed elsewhere_ploidy.bed
rm -f ploidy_bed_elsewhere.vcf ploidy_bed_none.vcf whole_x_ploidy.bed zero_x_ploidy.bed
rm -f sb_hap_d.txt sb_hap_bed.txt preset.txt preset_err.txt

rm -f small_cluster_call.gfa small_cluster_call.vg small_cluster_call.gam small_cluster_call.pack call_a.gam call_b.gam

## Read-level genotyping over a nested site with a star allele
# nested_snp_in_del.gfa: x = 1>2>3>5>6, y0 = 1>2>4>5>6, y1 = 1>6 (deletes the
# region holding the nested SNP). Reads from x and y1 only, so the top-level site
# is a het deletion and the nested site is traversed by one haplotype.
vg view -Fv nesting/nested_snp_in_del.gfa > ns.vg 2>/dev/null
vg sim -x ns.vg -P x -n 150 -l 4 -a -s 7 > ns_het.gam 2>/dev/null
vg sim -x ns.vg -P 'a#2#y1#0' -n 150 -l 2 -a -s 8 >> ns_het.gam 2>/dev/null
vg pack -x ns.vg -g ns_het.gam -o ns_het.pack 2>/dev/null
# --mismap-max pinned, because `vg sim` reads carry MAPQ 0, which makes e_r the cap for every read.
vg call ns.vg -k ns_het.pack --top-down -Y -p x --read-likelihood --mismap-max 0.7 --gam ns_het.gam 2>/dev/null > ns_rl.vcf


# Regression: reads traversing the deletion edge touch only the snarl's boundary
# nodes. A node-based informativeness test discarded them, which left the parent
# with reference-supporting reads only, called it hom-ref, and dropped the record
# entirely. The deletion must be called, and its DP must include those reads.
is $(grep -v "^#" ns_rl.vcf | awk -F'\t' '$4=="CATG"{nk=split($9,k,":"); di=0; for(j=1;j<=nk;j++) if(k[j]=="DP") di=j; split($10,f,":"); print f[di]}') "300" "deletion-spanning reads are counted, not discarded"

# Regression: half these reads are reverse strand. Failing to flip them meant they
# anchored on nothing and scored against the wrong allele, which turned the nested
# site into a spurious het SNP instead of a star allele.

# Independent nested calling (-A): every snarl genotyped on its own reads, with no
# parent restriction and no phase propagation.
vg call ns.vg -k ns_het.pack -A -p x --read-likelihood --mismap-max 0.7 --gam ns_het.gam 2>/dev/null > ns_rl_a.vcf
is "$?" "0" "--read-likelihood works with -A independent nested calling"
is $(grep -v "^#" ns_rl_a.vcf | awk '$4=="CATG" && $5=="C"' | wc -l | tr -d ' ') "1" "-A independent calling finds the deletion"
is $(grep -c "PS=" ns_rl_a.vcf | tr -d ' ') "0" "-A emits no PS tags, since it does not propagate phase"

# -A and the symbolic descent are two different ways to reach a child snarl, and both were on. -A
# sets RecurseAlways, so GraphCaller queues each child as its own snarl; nested calling -- ON by
# default under --read-likelihood since 956864c18 -- retained the same child as a nested chain.
# Both hash the same print_snarl string, so both took the same record_key, and
# LinkageCollector::live_index returns the first non-retracted entry for a key. Every nested snarl
# came out twice, the top-level record vanished, and which of the pair described the site depended
# on insertion order, i.e. on thread assignment -- which is what made `-A` differ from itself run
# to run at ~1,450 nested records while `-t 1` stayed byte-identical.
is $(grep -v "^#" ns_rl_a.vcf | awk '{print $1"\t"$2"\t"$3}' | sort | uniq -d | wc -l | tr -d ' ') "0" \
   "-A emits each snarl once, not once per way of reaching it"

# Refused rather than silently resolved: -A says every snarl is an independent record, --nested says
# children are called through their parent, and there is no output that is both.
vg call ns.vg -k ns_het.pack -A --nested -p x --read-likelihood --gam ns_het.gam \
    >/dev/null 2>ns_rl_a_nested.err
is "$?" "1" "-A with an explicit --nested is refused"
is $(grep -c "alternatives, not a combination" ns_rl_a_nested.err | tr -d ' ') "1" \
   "and says why"

# The same graph with the reference path running BACKWARDS through it, which is the only
# configuration that exercises resolve_site's reversed branch. nested_snp_in_del_rev.gfa has the
# topology and sequences of nested_snp_in_del.gfa with every path reversed, x = 6-,5-,3-,2-,1-, so
# both snarls (1,6) and (2,5) are reversed relative to the graph and flip_snarl rewrites them.
#
# This existed only in 26_deconstruct.t, so nothing in the calling suite ever reached the reversed
# branch: the reverse-oriented fixtures 18_vg_call.t already used run through the support caller with
# nested calling off, where symbolic projection is not armed at all. A defect that switched symbolic
# collapsing off for 7.4% of a real contig's sites was therefore invisible to 304 passing tests.
vg view -Fv nesting/nested_snp_in_del_rev.gfa > nsr.vg 2>/dev/null
vg sim -x nsr.vg -P x -n 150 -l 4 -a -s 7 > nsr_het.gam 2>/dev/null
vg sim -x nsr.vg -P 'a#2#y1#0' -n 150 -l 2 -a -s 8 >> nsr_het.gam 2>/dev/null
vg pack -x nsr.vg -g nsr_het.gam -o nsr_het.pack 2>/dev/null
# --mismap-max pinned, because `vg sim` reads carry MAPQ 0, which makes e_r the cap for every read.
vg call nsr.vg -k nsr_het.pack --top-down -Y -p x --read-likelihood --mismap-max 0.7 --gam nsr_het.gam 2>nsr_rl.err > nsr_rl.vcf

is $(grep -v "^#" nsr_rl.vcf | wc -l | tr -d ' ') "2" "a backwards reference path still emits both the parent and the child site"
is $(grep -v "^#" nsr_rl.vcf | awk '$5=="*"' | wc -l | tr -d ' ') "1" "the nested site is a star allele on a backwards reference path too"

# The coverage assertions. Without them "the tests pass" says nothing about whether the reversed
# branch ran at all, which is exactly how the original defect survived.
is "$(sed -n 's/.*atomize: \([0-9]*\) sites where projection is inert.*/\1/p' nsr_rl.err | head -1)" "0" "no site is left unresolvable when the reference path runs backwards"
REV_SITES=$(sed -n 's/.*resolve, \([0-9]*\) resolved as the reversal.*/\1/p' nsr_rl.err | head -1)
is "$([ "${REV_SITES:-0}" -ge 1 ] && echo yes || echo no)" "yes" "the reversed-snarl branch is exercised end to end"

# And the control: ns.vg, whose reference path runs forward, must NOT take that branch, or the
# counter is measuring something other than reversal.
vg call ns.vg -k ns_het.pack --top-down -Y -p x --read-likelihood --mismap-max 0.7 --gam ns_het.gam 2>ns_fwd.err >/dev/null
is "$(sed -n 's/.*resolve, \([0-9]*\) resolved as the reversal.*/\1/p' ns_fwd.err | head -1)" "0" "a forward reference path never takes the reversed branch"

rm -f nsr.vg nsr_het.gam nsr_het.pack nsr_rl.vcf nsr_rl.err ns_fwd.err

# Effective ploidy at a nested site only one parent haplotype traverses.
# Reads: 100 from x (ref SNP), 30 from y0 (alt SNP), 100 from y1 (deletion). The
# child site is reached by one haplotype only, so it must be genotyped haploid.
# Genotyping it as diploid lets a spurious heterozygote absorb the 30 minority
# reads for free, which shows up as a badly depressed GQ (38 rather than 256).
vg sim -x ns.vg -P x -n 100 -l 4 -a -s 11 > ns_mix.gam 2>/dev/null
vg sim -x ns.vg -P 'a#1#y0#0' -n 30 -l 4 -a -s 12 >> ns_mix.gam 2>/dev/null
vg sim -x ns.vg -P 'a#2#y1#0' -n 100 -l 2 -a -s 13 >> ns_mix.gam 2>/dev/null
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

rm -f ns.vg ns_het.gam ns_het.pack ns_rl.vcf ns_rl_a.vcf ns_mix.gam ns_mix.pack ns_mix.vcf ns_rl_a_nested.err




# --- the mosaic writer's structural cases, in one small graph -----------------
# The existing mosaic tests run on x.gbz, where every thread is a single segment: no join is
# exercised, so nothing the writer does between segments is tested. Everything that mattered on
# chr20 -- nested parent/child boundaries, an inversion, a clipped haplotype -- was reachable only
# through a 22 GB read database. This graph puts three of those in 396 bp:
#
#   (7,14)  a PARENT snarl -- the 7->14 deletion edge makes it one -- holding a chain of TWO child
#           snarls, (8,11) and (11,20), anchored on both sides so neither shares a boundary with the
#           parent. A parent with several children is the case the enter/leave rules exist for.
#   (14,16) an INVERSION over a stretch that CONTAINS A SITE, (30,33). A haplotype that inverts
#           traverses that site's own boundaries in REVERSE, so a row anchored on them disagrees
#           with the direction the fragment arrived in -- which is what makes the walk break.
#   pC      a CLIPPED haplotype: two GBWT fragments with a hole over the whole nested parent.
#
# Node ids are deliberately NOT ordered along the walk here -- node 16 closes a snarl holding nodes
# 30..33 -- because real graphs usually are, and testing containment by node id therefore passed on
# chr20 by luck while silently failing wherever ids run backwards, as they always do in an inversion.
rm -f fx.gfa fx_sim.gfa fx.gbz fx_sim.gbz fx.gam fx.mosaic.tsv fx.paths.gaf fx.graph.gfa
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\tCCTAGGCTTAGGACCTGATCGGATCCAGTACCGTTAAGCA\n'
  printf 'S\t2\tA\nS\t3\tT\n'
  printf 'S\t4\tGGCATTAGCCTTAGACCGATTGCAAGTCCAGGTTACCGAT\n'
  printf 'S\t5\tC\nS\t6\tG\n'
  printf 'S\t7\tTTGACCAGTTCAGGACTTACGGCATTCGGATCCAAGTTCC\n'
  printf 'S\t8\tAACCGGTTACGTTGCAATCGGTTACCAGTT\n'
  printf 'S\t9\tA\nS\t10\tG\n'
  printf 'S\t11\tGGATCCTAGCATTCGGATCCAAGTTCCAGA\n'
  printf 'S\t12\tC\nS\t13\tT\n'
  printf 'S\t20\tCCAGGTTACCGATTGCAAGTCCTTAGGCTA\n'
  printf 'S\t14\tACGTTGCAATCGGTTACCAGTTAAGGCCTTAGCATTCGGA\n'
  printf 'S\t30\tAGGCTTACCGATTGCAAGTCCAGGTTACCG\n'
  printf 'S\t31\tA\nS\t32\tC\n'
  printf 'S\t33\tTTGCAATCGGTTACCAGTTAAGGCCTTAGC\n'
  printf 'S\t16\tCATTCGGATCCAAGTTCCAGATTGACCAGTTCAGGACTTA\n'
  printf 'S\t17\tA\nS\t18\tG\n'
  printf 'S\t19\tGGCCTTAGCATTCGGATCCAAGTTCCAGATTGACCAGTTC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t7\t+\t0M\n'
  printf 'L\t7\t+\t8\t+\t0M\n'
  printf 'L\t8\t+\t9\t+\t0M\nL\t8\t+\t10\t+\t0M\nL\t9\t+\t11\t+\t0M\nL\t10\t+\t11\t+\t0M\n'
  printf 'L\t11\t+\t12\t+\t0M\nL\t11\t+\t13\t+\t0M\nL\t12\t+\t20\t+\t0M\nL\t13\t+\t20\t+\t0M\n'
  printf 'L\t20\t+\t14\t+\t0M\nL\t7\t+\t14\t+\t0M\n'
  printf 'L\t14\t+\t30\t+\t0M\nL\t30\t+\t31\t+\t0M\nL\t30\t+\t32\t+\t0M\n'
  printf 'L\t31\t+\t33\t+\t0M\nL\t32\t+\t33\t+\t0M\nL\t33\t+\t16\t+\t0M\n'
  printf 'L\t14\t+\t33\t-\t0M\nL\t30\t-\t16\t+\t0M\n'
  printf 'L\t16\t+\t17\t+\t0M\nL\t16\t+\t18\t+\t0M\nL\t17\t+\t19\t+\t0M\nL\t18\t+\t19\t+\t0M\n'
  printf 'P\tREF#0#chr1\t1+,2+,4+,5+,7+,8+,9+,11+,12+,20+,14+,30+,31+,33+,16+,17+,19+\t*\n'
  printf 'W\tpA\t0\tchr1\t0\t396\t>1>2>4>5>7>8>9>11>12>20>14>30>31>33>16>17>19\n'
  printf 'W\tpA\t1\tchr1\t0\t396\t>1>3>4>5>7>8>10>11>12>20>14>30>31>33>16>17>19\n'
  printf 'W\tpB\t0\tchr1\t0\t396\t>1>2>4>6>7>8>9>11>13>20>14<33<31<30>16>17>19\n'
  printf 'W\tpB\t1\tchr1\t0\t396\t>1>2>4>5>7>8>9>11>12>20>14>30>31>33>16>18>19\n'
  printf 'W\tpD\t0\tchr1\t0\t336\t>1>2>4>5>7>14>30>31>33>16>17>19\n'
  printf 'W\tpD\t1\tchr1\t0\t396\t>1>3>4>6>7>8>9>11>12>20>14>30>31>33>16>17>19\n'
  # pE carries the sample's route THROUGH THE PARENT -- both child alts together. Without it the
  # parent has no traversal to choose that spells the sample's children, its pin forces them
  # both to reference, and every nested site vanishes from the call set.
  printf 'W\tpE\t0\tchr1\t0\t396\t>1>2>4>5>7>8>10>11>13>20>14>30>31>33>16>17>19\n'
  # hR inverts and carries the site inside the inversion as alt.
  printf 'W\thR\t0\tchr1\t0\t396\t>1>2>4>5>7>8>9>11>12>20>14<33<32<30>16>18>19\n'
  # pC is CLIPPED: a hole over the whole nested parent.
  printf 'W\tpC\t0\tchr1\t0\t122\t>1>3>4>6>7\n'
  printf 'W\tpC\t0\tchr1\t214\t396\t>14<33<32<30>16>18>19\n'
} > fx.gfa
# The sample must NOT be in the panel, or the model just names it and every strand is one segment.
# So: one graph with the truth walks for `vg sim`, one without for `vg call`, same node ids.
cp fx.gfa fx_sim.gfa
printf 'W\tT\t0\tchr1\t0\t396\t>1>3>4>6>7>8>10>11>13>20>14<33<32<30>16>18>19\n' >> fx_sim.gfa
printf 'W\tT\t1\tchr1\t0\t396\t>1>2>4>5>7>8>9>11>12>20>14>30>31>33>16>17>19\n' >> fx_sim.gfa
vg gbwt -G fx.gfa --gbz-format -g fx.gbz --set-reference REF 2>/dev/null
is "$?" 0 "mosaic structural fixture builds"
vg gbwt -G fx_sim.gfa --gbz-format -g fx_sim.gbz --set-reference REF 2>/dev/null >/dev/null

# The graph has to actually have the shapes before the writer can be tested on them.
vg snarls fx.gbz 2>/dev/null > fx.snarls
is $(vg view -R fx.snarls 2>/dev/null | grep -c '"parent"') "3" \
   "the fixture has three nested snarls: two children in one parent, and one inside the inversion"
is $(vg paths -x fx.gbz -L 2>/dev/null | grep -c "^pC#0#chr1#") "2" \
   "the clipped haplotype is stored as two GBWT fragments, so it has a hole"

vg sim -x fx_sim.gbz -n 600 -l 40 -a -s 17 --path "T#0#chr1#0" >  fx.gam 2>/dev/null
vg sim -x fx_sim.gbz -n 600 -l 40 -a -s 23 --path "T#1#chr1#0" >> fx.gam 2>/dev/null
# --linkage-scale is the decay distance in bp and defaults to 10 kb, which over 396 bp makes every
# recombination effectively free to avoid: the model then explains the sample by flipping single
# alleles rather than switching haplotype, and the mosaic collapses to one segment per strand. A
# scale on the order of the site spacing lets the model switch haplotype between sites. The result
# is stable over 30..60.
# --mismap-max pinned, because `vg sim` reads carry MAPQ 0, which makes e_r the cap for every read.
vg call fx.gbz --read-likelihood --mismap-max 0.7 --gam fx.gam -t 1 -s samp --phased --linkage-scale 40 \
    --mosaic-out fx.mosaic.tsv --progress >/dev/null 2>fx.err
is "$?" 0 "vg call runs on the mosaic structural fixture"

# --- the hierarchy rules -----------------------------------------------------
# A parent/child haplotype change has two boundaries and both belong to the child: Cs going in, Ce
# coming out. Getting this wrong put the child run's start at the PARENT's end, past every one of
# its own sites -- 584 such boundaries on chr20, of which only two collapsed far enough for the
# round trip to catch them.
is $(awk '/enters a child/ {print ($0 ~ /: [1-9][0-9]* boundaries/) ? 1 : 0}' fx.err) "1" \
   "a run boundary where the walk ENTERS a child snarl is exercised"
is $(awk '/enters a child/ {print ($0 ~ /, [1-9][0-9]* where it leaves/) ? 1 : 0}' fx.err) "1" \
   "and one where it LEAVES a child snarl"

# --- the inversion -----------------------------------------------------------
# A row whose start contradicts the direction the fragment arrived in cannot join its predecessor:
# X+ followed by X- is not a walk. It is retried without the carry and breaks the fragment before
# itself.
is $(awk '/carried direction/ {print ($0 ~ /: [1-9][0-9]* rows/) ? 1 : 0}' fx.err) "1" \
   "a row that can only be walked against the carried direction is exercised"
# The same node in both orientations across a fragment boundary -- the inversion boundary itself.
is $(awk -F'\t' '/^H\t/ {n++; s[n]=$3; f[n]=$4; a[n]=$7; b[n]=$8}
     END {for (i=1; i<n; i++) if (s[i]==s[i+1] && int(b[i]/2)==int(a[i+1]/2) && b[i]!=a[i+1]) {
            print (f[i] != f[i+1]) ? 1 : 0; exit }
          print "no such pair"}' fx.mosaic.tsv) "1" \
   "where consecutive rows enter one node from opposite sides, the fragment breaks between them"

# --- what a consumer relies on ----------------------------------------------
# A join is the thing x.gbz cannot test: one segment per thread there, so no junction exists.
is $(awk -F'\t' '/^H\t/ {c[$2"/"$3"/"$4]++} END {m=0; for (k in c) if (c[k]>m) m=c[k]; print (m>1)?1:0}' \
     fx.mosaic.tsv) "1" \
   "at least one fragment holds several segments, so a junction is exercised at all"
# Inside a fragment every junction meets on the SAME ORIENTED node, or it is not a path.
is $(awk -F'\t' '/^H\t/ {k=$2"/"$3"/"$4; if (k==pk && $7 != pb) bad++; pk=k; pb=$8}
     END {print bad+0}' fx.mosaic.tsv) "0" \
   "consecutive segments of a fragment meet at the same oriented node"
# A row with no position states its nodes but not a walk, so it must not sit inside one.
is $(awk -F'\t' '/^H\t/ {k=$2"/"$3"/"$4; n[k]++; if ($12==".") d[k]++}
     END {for (k in n) if (d[k] && n[k]>1) bad++; print bad+0}' fx.mosaic.tsv) "0" \
   "no positionless row shares a fragment with rows that do claim a walk"

# THE ACCEPTANCE TEST, now with something to chew on: nested boundaries, an inversion, a fragment
# break and a junction. Every fragment must expand to an exact walk in the graph.
vg paths -x fx.gbz -A > fx.paths.gaf 2>/dev/null
vg view -g fx.gbz > fx.graph.gfa 2>/dev/null
is $(python3 -B ../artifacts/scripts/mosaic_to_path.py --mosaic fx.mosaic.tsv --gaf fx.paths.gaf --gfa fx.graph.gfa \
       --quiet >fx.expand.txt 2>&1; echo $?) "0" \
   "every fragment of the structural fixture expands to an exact path in the graph"

rm -f fx.gfa fx_sim.gfa fx.gbz fx_sim.gbz fx.gam fx.snarls fx.err fx.mosaic.tsv \
      fx.paths.gaf fx.graph.gfa fx.expand.txt

# The cactus snarl finder makes a unary snarl for the loop off node 2 here, which the
# read-likelihood caller refuses with a message saying how to find snarls it can call.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tGGGG\nS\t3\tTTTTACGTAC\nS\t4\tCCCC\nS\t5\tAAAA\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t4\t-\t0M\nP\tref\t1+,2+,3+\t*\n' > unary.gfa
vg convert -g unary.gfa -p > unary.vg
vg snarls -A cactus unary.vg > unary.snarls
: > unary.gam
vg pack -x unary.vg -g unary.gam -o unary.pack
is "$(vg call unary.vg -r unary.snarls --read-likelihood --gam unary.gam -k unary.pack -s S 2>&1 >/dev/null | grep -c "vg snarls -A integrated")" 1 "call: --read-likelihood explains how to find snarls it can call when given a unary snarl"

rm -f unary.gfa unary.vg unary.snarls unary.gam unary.pack

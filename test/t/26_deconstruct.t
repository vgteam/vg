#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 119

vg mod -U 10 msgas/hla_v.vg | vg mod -c - > hla_v.vg
vg index hla_v.vg -x hla.xg

vg deconstruct hla.xg -p "gi|157734152:29563108-29564082" > hla_decon.vcf
is $(grep -v "#" hla_decon.vcf | wc -l) 17 "deconstructed hla vcf has correct number of sites"
is $(grep -v "#" hla_decon.vcf | grep 822 | awk '{print $4 "-" $5}') "C-CGCGGGCGCCGTGGATGGAGCA" "deconstructed hla vcf has correct insertion"
vg deconstruct hla.xg -p "gi|568815592:29791752-29792749" > hla_decon.vcf
is $(grep -v "#" hla_decon.vcf | wc -l) 17 "deconstructed hla vcf with other path has correct number of sites"
is $(grep -v "#" hla_decon.vcf | grep 824 | awk '{print $4 "-" $5}') "CGCGGGCGCCGTGGATGGAGCA-C" "deconstructed hla vcf has correct deletion"

vg deconstruct hla.xg -p "gi|568815592:29791752-29792749" -e > hla_decon_path.vcf
grep -v "#" hla_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' | sort > hla_decon.tsv
grep -v "#" hla_decon_path.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' | sort > hla_decon_path.tsv
diff hla_decon.tsv hla_decon_path.tsv
is "$?" 0 "path-based and exhaustive decontruction give equivalent sites when expected"

# want to extract a sample, but bcftools -s doesn't seem to work on travis.  so we torture it out with awk
SAMPLE_COL=$(grep CHROM hla_decon_path.vcf | tr '\t' '\n' | nl | grep "528476637" | awk '{print $1}')
is $(grep -v "#" hla_decon_path.vcf | awk -v x="$SAMPLE_COL" '{print $x}' | uniq) 1 "path that differs from reference in every alt has correct genotype"

SAMPLE_COL=$(grep CHROM hla_decon_path.vcf | tr '\t' '\n' | nl | grep "568815564" | awk '{print $1}')
is $(grep -v "#" hla_decon_path.vcf | awk -v x="$SAMPLE_COL" '{print $x}' | uniq) 0 "path that is same as reference in every alt has correct genotype"

is $(grep "#" hla_decon_path.vcf | grep "568815592") "##contig=<ID=gi|568815592:29791752-29792749,length=998>" "reference contig correctly written"


rm -f hla_decon.vcf hla_decon_path.vcf  hla_decon.tsv hla_decon_path.tsv hla_v.vg hla.xg

cp sv/x.inv.gfa inv.gfa
printf "P\ty\t1+,2-,3+\t9M,20M,21M\n" >> inv.gfa
vg view -Fv inv.gfa > inv.vg
vg index inv.vg -x inv.xg
vg deconstruct inv.xg -p x -e > inv_decon.vcf
grep -v "#" inv_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > inv_decon.tsv
printf "x\t10\tCTTGGAAATTTTCTGGAGTT\tAACTCCAGAAAATTTCCAAG\t1\n" > inv_truth.tsv
diff inv_decon.tsv inv_truth.tsv
is "$?" 0 "deconstruct correctly handles a simple inversion"

rm -f inv_decon.vcf inv_decon.tsv inv_truth.tsv

vg deconstruct inv.xg -p y -e > inv_decon.vcf
grep -v "#" inv_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > inv_decon.tsv
printf "y\t10\tAACTCCAGAAAATTTCCAAG\tCTTGGAAATTTTCTGGAGTT\t1\n" > inv_truth.tsv
diff inv_decon.tsv inv_truth.tsv
is "$?" 0 "deconstruct correctly handles a simple inversion when the reference contains the reversing edge"

vg gbwt -G inv.gfa -g inv.chop.gbz --gbz-format --max-node 5
vg deconstruct inv.chop.gbz -p x -O > inv.chop.decon
vg gbwt -G inv.gfa -g inv.gbz --gbz-format --max-node 1025
vg deconstruct inv.gbz -p x -O > inv.decon
diff inv.decon inv.chop.decon
is "$?" 0 "deconstruct applies translation from gbz with -O"

rm -f inv.gfa inv.vg inv.xg inv_decon.vcf inv_decon.tsv inv_truth.tsv inv.gbz inv.chop.gbz inv.gbz inv.chop.decon inv.decon


vg construct -v small/x.vcf.gz -r small/x.fa | vg view -g - > cyclic_small.gfa
 printf "L\t33\t+\t30\t+\t0M\n" >> cyclic_small.gfa
 printf "P\ty\t1+,3+,5+,6+,8+,9+,11+,12+,14+,15+,17+,18+,20+,21+,23+,24+,26+,27+,29+,30+,32+,33+,30+,31+,33+,35+,36+,38+,40+,41+,43+,44+,46+,47+,49+,50+,52+,53+,55+,56+,57+,60+,61+,62+,64+,65+,67+,68+,70+,71+,73+,74+,76+,78+,79+,81+,83+,84+,86+,87+,89+,90+,92+,94+,95+,97+,98+,100+,101+,102+,103+,104+,106+,107+,109+,110+,112+,114+,115+,117+,118+,120+,122+,123+,124+,126+,127+,129+,130+,132+,133+,135+,136+,137+,139+,141+,142+,144+,145+,147+,148+,149+,151+,152+,154+,155+,157+,158+,159+,160+,162+,163+,165+,166+,168+,169+,171+,172+,174+,176+,177+,179+,181+,182+,184+,185+,187+,188+,190+,191+,193+,195+,196+,198+,199+,201+,202+,204+,205+,206+,207+,209+,210+,211+,212+,214+,215+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,12M,1M,6M,32M,9M,1M,2M,1M,18M,1M,19M,1M,17M,19M,1M,17M,1M,12M,5M,1M,8M,1M,12M,1M,3M,1M,13M,1M,2M,1M,32M,18M,1M,1M,3M,1M,9M,1M,6M,1M,2M,1M,14M,1M,1M,32M,1M,1M,5M,1M,19M,1M,23M,1M,1M,6M,1M,2M,1M,32M,26M,1M,6M,1M,14M,1M,11M,22M,1M,9M,1M,20M,12M,1M,32M,30M,1M,17M,1M,6M,1M,14M,1M,1M,11M,6M,1M,9M,1M,15M,1M,32M,10M,1M,17M,1M,1M,1M,12M,1M,4M,1M,17M,1M,9M,1M,3M,1M,30M,1M,1M,1M,16M,1M,10M,1M,10M,1M,1M,1M,14M,1M,1M,5M,1M,1M,1M,3M,1M,10M,1M,4M,1M,27M,2M,25M,1M,1M\n" >> cyclic_small.gfa
vg view -Fv cyclic_small.gfa > cyclic_small.vg
vg index cyclic_small.vg -x cyclic_small.xg
vg deconstruct cyclic_small.xg -p y -e > cyclic_small_decon.vcf
grep -v "#" cyclic_small_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > cyclic_small_decon.tsv
printf "y\t121\tTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAACGTTTGACAATCTATCACCAGGGGTAATGTGGGGAAA\tTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA\t1\n" > cyclic_small_truth.tsv
diff cyclic_small_decon.tsv cyclic_small_truth.tsv
is "$?" 0 "deconstruct correctly handles a cycle in the reference path when contained inside snarl"

rm -f cyclic_small_decon.vcf cyclic_small_decon.tsv cyclic_small_truth.tsv

vg deconstruct cyclic_small.xg -p x -e > cyclic_small_decon.vcf
grep -v "#" cyclic_small_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > cyclic_small_decon.tsv
printf "x\t121\tTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA\tTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAACGTTTGACAATCTATCACCAGGGGTAATGTGGGGAAA\t1\n" > cyclic_small_truth.tsv
diff cyclic_small_decon.tsv cyclic_small_truth.tsv
is "$?" 0 "deconstruct correctly handles a cycle in the alt path"

rm -f cyclic_small_decon.vcf cyclic_small_decon.tsv cyclic_small_truth.tsv cyclic_small.gfa cyclic_small.vg cyclic_small.xg

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view -g - > cyclic_tiny.gfa
printf "L\t12\t+\t9\t+\t0M\n" >> cyclic_tiny.gfa
printf "P\ty\t1+,3+,5+,6+,8+,9+,11+,12+,9+,10+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,19M,1M,4M,1M,11M\n" >> cyclic_tiny.gfa
vg view -Fv cyclic_tiny.gfa > cyclic_tiny.vg
vg index cyclic_tiny.vg -x cyclic_tiny.xg
vg find -x cyclic_tiny.xg  -n 10 -n 11 -n 12 -n 13 -n 14 -n 15 -c 1 > cycle.vg
# TODO: Make deconstruct see through subpaths to the base path
vg view cycle.vg | sed 's/\([xy]\)\[[-0-9]*\]/\1/g' >cycle-asfullpaths.gfa
vg index cycle-asfullpaths.gfa -x cycle.xg
vg deconstruct cycle.xg -p y -e -t 1 > cycle_decon.vcf
is $(grep -v "#" cycle_decon.vcf | wc -l) 1 "cyclic reference deconstruction has correct number of variants"
grep -v "#" cycle_decon.vcf | grep 20 | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > cycle_decon.tsv
grep -v "#" cycle_decon.vcf | grep 44 | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' >> cycle_decon.tsv
printf "y\t44\tA\tT\t1\n" > cycle_decon_truth.tsv
diff cycle_decon.tsv cycle_decon_truth.tsv
is "$?" 0 "deconstruct correctly handles cycle in the reference path that spans snarl"

rm -f cyclic_tiny_decon.vcf cyclic_tiny_decon.tsv cyclic_tiny_truth.tsv cyclic_tiny.gfa cyclic_tiny.vg cyclic_tiny.xg
rm -f cycle.vg cycle.xg cycle_decon.vcf cycle_decon.tsv cycle_decon_truth.tsv

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view -g - > tiny_names.gfa
printf "P\tref.1\t1+,3+,5+,6+,8+,9+,11+,12+,14+,15+\t*\n" >> tiny_names.gfa
printf "P\talt1.1\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t*\n" >> tiny_names.gfa
printf "P\talt1.2\t1+,2+,4+,6+,7+,9+,11+,12+,14+,15+\t*\n" >> tiny_names.gfa
printf "P\talt2.3\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t*\n" >> tiny_names.gfa
printf "P\talt2.4\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t*\n" >> tiny_names.gfa
vg view -Fv tiny_names.gfa > tiny_names.vg
vg index tiny_names.vg -x tiny_names.xg
vg deconstruct tiny_names.xg -P ref -H . -e -d 1 | sort > tiny_names_decon.vcf
is $(grep -v "#" tiny_names_decon.vcf | wc -l) 2 "-P -H options return correct number of variants"
is $(grep -v "#" tiny_names_decon.vcf | grep ref.1 | wc -l) 2 "-P -H options use correct reference name"
is $(grep -v "#" tiny_names_decon.vcf | grep ref.1 | grep 14 | grep "CONFLICT=alt1" | wc -l) 0 "-P -H does not find conflict in alt1 in second variant"
vg deconstruct tiny_names.vg -P ref -H . -e -d 1 | sort > tiny_names_decon_vg.vcf
diff tiny_names_decon.vcf tiny_names_decon_vg.vcf
is "$?" 0 "deconstructing vg graph gives same output as xg graph"

rm -f tiny_names.gfa tiny_names.vg tiny_names.xg tiny_names_decon.vcf tiny_names_decon_vg.vcf

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg
vg index -x x.xg x.vg
vg gbwt -v small/x.vcf.gz -o x.gbwt -x x.vg
vg deconstruct x.xg -g x.gbwt | bgzip > x.decon.vcf.gz
tabix -f -p vcf  x.decon.vcf.gz
cat small/x.fa |  bcftools consensus small/x.vcf.gz -s 1 -H 1 > small.s1.h1.fa
cat small/x.fa |  bcftools consensus small/x.vcf.gz -s 1 -H 2 > small.s1.h2.fa
cat small/x.fa |  bcftools consensus x.decon.vcf.gz -s 1 -H 1 > decon.s1.h1.fa
cat small/x.fa |  bcftools consensus x.decon.vcf.gz -s 1 -H 2 > decon.s1.h2.fa
diff small.s1.h1.fa decon.s1.h1.fa
is "$?" 0 "haplotype 1 preserved when deconstructing small test with gbwt"
diff small.s1.h2.fa decon.s1.h2.fa
is "$?" 0 "haplotype 2 preserved when deconstructing small test with gbwt"

vg autoindex -r small/x.fa -v small/x.vcf.gz -w giraffe -p x
vg deconstruct x.giraffe.gbz > x.gbz.decon.vcf
gzip -dc x.decon.vcf.gz > x.decon.vcf
diff x.decon.vcf x.gbz.decon.vcf
is "$?" 0 "gbz deconstruction gives same output as gbwt deconstruction"

rm -f x.vg x.xg x.gbwt x.decon.vcf.gz x.decon.vcf.gz.tbi x.decon.vcf x.gbz.decon.vcf x.giraffe.gbz x.min x.dist small.s1.h1.fa small.s1.h2.fa decon.s1.h1.fa decon.s1.h2.fa

# todo you could argue merging shouldn't happen here because there's no child snarl
# this check should come into play with nesting support
vg construct -r small/x.fa -v small/x.vcf.gz | vg view - > small_cluster.gfa
printf "L\t1\t+\t9\t+\t0M\n" >> small_cluster.gfa
printf "P\ty\t1+,2+,4+,6+,7+,9+\t*\n" >> small_cluster.gfa
printf "P\tz\t1+,9+\t*\n" >> small_cluster.gfa
vg deconstruct small_cluster.gfa -p x > small_cluster_0.vcf
vg deconstruct small_cluster.gfa -p x -L 0.3 --cluster-min-len 0 > small_cluster_3.vcf
is "$(tail -1 small_cluster_0.vcf | awk '{print $5}')" "GATTTGA,G" "cluster-free deconstruction finds all alt alleles"
is "$(tail -1 small_cluster_3.vcf | awk '{print $5}')" "G" "clustered deconstruction finds fewer alt alleles"
is "$(tail -1 small_cluster_3.vcf | awk '{print $10}')" "0:0.333:0" "clustered deconstruction finds correct allele info"

# --cluster-min-len gates when the site's core length is below the threshold.
# small_cluster's only snarl has core length 6 bp.
vg deconstruct small_cluster.gfa -p x -L 0.3 --cluster-min-len 0 > small_cluster_3_off.vcf
is "$(tail -1 small_cluster_3_off.vcf | awk '{print $5}')" "G" "--cluster-min-len 0 is equivalent to clustering everywhere"
vg deconstruct small_cluster.gfa -p x -L 0.3 --cluster-min-len 50 > small_cluster_3_50.vcf
is "$(tail -1 small_cluster_3_50.vcf | awk '{print $5}')" "GATTTGA,G" "--cluster-min-len 50 gates clustering on a small site"
vg deconstruct small_cluster.gfa -p x -L 0.3 --cluster-min-len 6 > small_cluster_3_6.vcf
is "$(tail -1 small_cluster_3_6.vcf | awk '{print $5}')" "G" "--cluster-min-len at the site length still clusters"
vg deconstruct small_cluster.gfa -p x --cluster-min-len 50 2>/dev/null > small_cluster_min_only.vcf
diff small_cluster_0.vcf small_cluster_min_only.vcf
is "$?" 0 "--cluster-min-len without -L is a no-op"

rm -f small_cluster.gfa small_cluster_0.vcf small_cluster_3.vcf small_cluster_3_off.vcf small_cluster_3_50.vcf small_cluster_3_6.vcf small_cluster_min_only.vcf

# --cluster-min-len gates on CORE LENGTH -- the longest allele once the prefix and suffix shared by
# every allele are stripped -- not on the raw snarl interior.  Measuring the interior would gate a
# 1bp SNP on the size of whatever snarl happens to contain it, and would disagree with vg call,
# whose records are flattened down to an anchor base.
core_nalt() { vg deconstruct nesting/$1.gfa -p x -L 0.6 --cluster-min-len $2 2>/dev/null | awk -F'\t' '$1!~/^#/{printf "%d", split($5,a,",")}'; }

is "$(core_nalt core_snp_in_flanks 50)" "2" "a 1bp SNP in a large snarl is not clustered at --cluster-min-len 50"
is "$(core_nalt core_sv60 50)"         "1" "a 60bp SV is still clustered at --cluster-min-len 50"
is "$(core_nalt core_del61 50)"        "1" "a 61bp deletion is still clustered (the reference allele is measured)"
is "$(core_nalt core_ins49 50)"        "2" "a 49bp insertion is not clustered at 50 (the anchor base is not counted)"
is "$(core_nalt core_ins50 50)"        "1" "a 50bp insertion is clustered at 50"
# 50 is the DEFAULT, so -L on its own must gate identically.  Every metric assertion in this file
# passes --cluster-min-len 0 for exactly that reason.
default_nalt() { vg deconstruct nesting/$1.gfa -p x -L 0.6 2>/dev/null | awk -F'\t' '$1!~/^#/{printf "%d", split($5,a,",")}'; }
is "$(default_nalt core_snp_in_flanks)" "$(core_nalt core_snp_in_flanks 50)" "-L alone gates a small variant like --cluster-min-len 50"
is "$(default_nalt core_sv60)"          "$(core_nalt core_sv60 50)"          "-L alone still clusters an SV"

# A PURE DELETION -- a traversal straight from snarl start to snarl end -- used to keep its boundary
# handles, which made its handle set disjoint from every other allele's by construction.  Its
# similarity was therefore structurally 0 and no -L could ever merge it.  It is now scored against
# the site: two deletions that both remove 59 of a 60bp site score 59/60.
# The ladder is -L 0.99 / 0.983334 / 0.983333 / 0.6 / 0.1, printing the ALT count at each.
# --cluster-min-len 0 so this ladder measures the METRIC alone: the gate now defaults to 50, which
# would decide del1_vs_snp1 (a 1bp site) before the similarity was ever consulted.  The gate has its
# own ladder in core_nalt below.
dladder() { for L in 0.99 0.983334 0.983333 0.6 0.1 ; do
    vg deconstruct nesting/$1.gfa -p x -L $L --cluster-min-len 0 2>/dev/null |
      awk -F'\t' '$1!~/^#/{n=split($5,a,","); if($5=="."||$5=="")n=0; printf "%d",n}' ; done ; }

is "$(dladder del59_vs_del60)"     "22111" "a 59bp and a 60bp deletion cluster, flipping at 59/60"
is "$(dladder del60_vs_del59ins1)" "22111" "a deletion clusters with a deletion carrying a novel base"
is "$(dladder del60_vs_snp)"       "22222" "a 60bp deletion never clusters with a 60bp substitution"
is "$(dladder inv60_vs_del60)"     "22222" "a 60bp inversion never clusters with a 60bp deletion"
# only a PURE DELETION is scored against the site; scaling every pair by it would make unrelated
# alleles in a big snarl look alike, which these two pin
is "$(dladder inv60_in_2kb)"       "22222" "a 2kb snarl does not make an inversion clusterable"
is "$(dladder unrelated10_in_2kb)" "22222" "a 2kb snarl does not make two unrelated 10bp alleles clusterable"
# the small-variant case: deleting a base and substituting it are different events, so scoring a
# pure deletion against the site must not make them look alike
is "$(dladder del1_vs_snp1)"       "22222" "a 1bp deletion never clusters with a 1bp SNP"

# The flip point itself is pinned by dladder above.  What is left to pin here is that the clustering
# machinery is entirely invisible at the default: -L 1.0 must match no -L byte for byte, and neither
# may carry the TS/TL FORMAT fields that any -L < 1.0 adds (see the -L 0.99 run, which does not merge
# but does emit them).
vg deconstruct nesting/del59_vs_del60.gfa -p x 2>/dev/null > deldefault.vcf
vg deconstruct nesting/del59_vs_del60.gfa -p x -L 1.0 2>/dev/null > delnoop.vcf
diff deldefault.vcf delnoop.vcf
is "$?" 0 "-L 1.0 is byte-identical to no -L on a pure-deletion site"
is "$(grep -cE 'GT:TS:TL' deldefault.vcf)" "0" "the default emits no TS/TL, so -L leaves no trace when off"
is "$(vg deconstruct nesting/del59_vs_del60.gfa -p x -L 0.99 2>/dev/null | grep -cE 'GT:TS:TL')" "1" "-L < 1.0 emits TS/TL even when nothing merges"
rm -f deldefault.vcf delnoop.vcf

# -L in a nested run is lossy at the parent and precise at the children, exactly as in vg call:
# the parent gives the collapsed view of a large variant, the child records the detail.
vg deconstruct nesting/nested_snp_in_ins.gfa -p x -a -L 0.6 --cluster-min-len 0 >/dev/null 2>&1
is "$?" 0 "-L is allowed in a nested run (-a)"
# ...but not with -R.  A "*" in a child record means "an upstream deletion covers this site", and
# clustering can absorb the allele that deletion came from -- which one survives is decided by
# traversal order, so by path names -- leaving the "*" referring to nothing in the file.  That is a
# malformed record rather than a lossy one.  vg call refuses -L with -Y for the same reason.
vg deconstruct nesting/nested_snp_in_ins.gfa -p x -a -R -L 0.6 >/dev/null 2>&1
is "$?" 1 "-L is rejected with -R/--star-allele"
vg deconstruct nesting/nested_snp_in_ins.gfa -p x -a -R >/dev/null 2>&1
is "$?" 0 "-R without -L still works"

# Star alleles (-R) are a marker, not sequence.  They must survive verbatim as "*" in every
# orientation and at every site type: complement['*'] is 'N', so reverse-complementing one turns it
# into a real ambiguous base, and prepending the anchor base turns it into "A*"/"AN".
star_alt() { awk -F'\t' -v id="$2" '$3 == id {print $2" "$4" "$5}' "$1"; }

vg deconstruct nesting/nested_snp_in_del.gfa         -p x -a -R 2>/dev/null > star_fwd.vcf
vg deconstruct nesting/nested_snp_in_del_rev.gfa     -p x -a -R 2>/dev/null > star_rev.vcf
vg deconstruct nesting/nested_del_star_indel.gfa     -p x -a -R 2>/dev/null > star_indel.vcf
vg deconstruct nesting/nested_del_star_indel_rev.gfa -p x -a -R 2>/dev/null > star_rev_indel.vcf
vg deconstruct nesting/nested_mnp_star.gfa           -p x -a -R 2>/dev/null > star_mnp.vcf
vg deconstruct nesting/nested_ins_star_only.gfa      -p x -a -R 2>/dev/null > star_emptyref.vcf
vg deconstruct nesting/insertion_with_three_snps.gfa -p x -a -R 2>/dev/null > star_three.vcf
vg deconstruct nesting/star_cluster.gfa              -p x -a -R 2>/dev/null > star_cluster.vcf
vg deconstruct nesting/star_allele_cluster.gfa       -p x -a -R 2>/dev/null > star_allele_cluster.vcf

# The *_rev fixtures are keyed on '<5<2', not '>2>5': a snarl ID is spelled in the orientation the
# reference traverses it, so a reverse-oriented site reads backwards.  Only the ID differs -- the
# POS/REF/ALT asserted here are the same shape as the forward case.
is "$(star_alt star_fwd.vcf       '>2>5')" "3 T A,*"      "star is * at a forward substitution site"
is "$(star_alt star_rev.vcf       '<5<2')" "3 A T,*"      "star is not reverse complemented into N"
is "$(star_alt star_indel.vcf     '>2>5')" "2 ATTTT AA,*" "star is not padded with the previous base"
is "$(star_alt star_rev_indel.vcf '<5<2')" "2 CAAAA CT,*" "star survives reversal and padding together"
is "$(star_alt star_mnp.vcf       '>2>5')" "3 TT AA,*"    "a star does not force indel padding of the real alleles"
is "$(star_alt star_emptyref.vcf  '>2>5')" "2 A *"        "empty REF keeps its anchor base with a star-only ALT"

# -R adds star alleles; it must not move or respell the real ones
diff <(vg deconstruct nesting/nested_mnp_star.gfa -p x -a -R 2>/dev/null | grep -v '^#' | cut -f1,2,4) \
     <(vg deconstruct nesting/nested_mnp_star.gfa -p x -a    2>/dev/null | grep -v '^#' | cut -f1,2,4)
is "$?" 0 "-R does not change CHROM/POS/REF of the real alleles"

# Negative sweep over every fixture above: an AT entry of "." identifies a star allele, and AT is
# Number=R so the mapping is by index, independent of ALT order.
star_fields() {
    grep -hv "^#" "$@" | awk -F'\t' '{
        n = split($5, alts, ",");
        if (match($8, /AT=[^;]*/)) {
            split(substr($8, RSTART + 3, RLENGTH - 3), ats, ",");
            for (i = 1; i <= n; i++) if (ats[i + 1] == ".") print alts[i];
        }}'
}
STAR_VCFS="star_fwd.vcf star_rev.vcf star_indel.vcf star_rev_indel.vcf star_mnp.vcf star_emptyref.vcf star_three.vcf star_cluster.vcf star_allele_cluster.vcf"
is "$(star_fields $STAR_VCFS | grep -c '')"      "8" "the star allele sweep is not vacuous"
is "$(star_fields $STAR_VCFS | grep -vc '^\*$')" "0" "every allele whose AT entry is . is spelled exactly *"

rm -f $STAR_VCFS

# Nesting tests now use a two-step process:
# 1. Compute gref cover with vg paths
# 2. Run vg deconstruct against the gref reference with -a to use the gref paths

# Test: SNP inside deletion
vg paths --compute-gref --min-gref-len 0 -x nesting/nested_snp_in_del.gfa -Q x > nested_snp_in_del.gref.pg
vg deconstruct nested_snp_in_del.gref.pg -P gref_x -a > nested_snp_in_del.vcf
grep -v ^# nested_snp_in_del.vcf | awk '{print $4 "\t" $5 "\t" $10}' > nested_snp_in_del.tsv
printf "CATG\tCAAG,C\t1|2\n" > nested_snp_in_del_truth.tsv
printf "T\tA\t1|.\n" >> nested_snp_in_del_truth.tsv
diff nested_snp_in_del.tsv nested_snp_in_del_truth.tsv
is "$?" 0 "nested deconstruction gets correct allele for snp inside deletion"

# The other half of the CH story.  The reference allele of the deletion spells out the
# deleted sequence, so the nested SNP has coordinates on the same contig: it stays at LV=1
# and takes no contig hop.  nested_snp_in_ins below is the same SNP inside an insertion,
# where it lands on a gref contig at LV=0 with CH=1.  Nothing else distinguishes the two.
is "$(grep -v ^# nested_snp_in_del.vcf | grep -o 'LV=[0-9]*' | tr '\n' ' ')" "LV=0 LV=1 " "a SNP inside a deletion stays nested on its own contig"
is "$(grep -v ^# nested_snp_in_del.vcf | grep -o 'CH=[0-9]*' | tr '\n' ' ')" "CH=0 CH=0 " "a SNP inside a deletion takes no contig hop"
is $(grep -v ^# nested_snp_in_del.vcf | cut -f1 | sort -u | wc -l) 1 "both records are on the base contig"

rm -f nested_snp_in_del.gref.pg nested_snp_in_del.vcf nested_snp_in_del.tsv nested_snp_in_del_truth.tsv

# Test: SNP inside insertion with LV field checks
vg paths --compute-gref --min-gref-len 0 -x nesting/nested_snp_in_ins.gfa -Q x > nested_snp_in_ins.gref.pg
vg deconstruct nested_snp_in_ins.gref.pg -P gref_x -a > nested_snp_in_ins.vcf
is $(grep -c "^##INFO=<ID=LV," nested_snp_in_ins.vcf) 1 "LV is declared in the header"
is $(grep -c "^##INFO=<ID=CH," nested_snp_in_ins.vcf) 1 "CH is declared in the header"
grep -v ^# nested_snp_in_ins.vcf | awk '{print $4 "\t" $5 "\t" $10}' > nested_snp_in_ins.tsv
# With -P gref_x, nested variants are on gref contigs (gref_x_1_alt), parent on gref_x
# So order is: insertion (on gref_x) then SNP (on gref_x_1_alt)
printf "C\tCAAG,CATG\t1|2\n" > nested_snp_in_ins_truth.tsv
printf "A\tT\t0|1\n" >> nested_snp_in_ins_truth.tsv
diff nested_snp_in_ins.tsv nested_snp_in_ins_truth.tsv
is "$?" 0 "nested deconstruction gets correct allele for snp inside insert"

# A SNP inside an insertion has no coordinates on the base contig at all, so its record
# lands on a gref contig, and LV counts only ancestors on the record's own contig.  So the
# SNP is top-level where it lives (LV=0) and CH=1 records that the step crossed into
# non-reference sequence.  Compare with nested_snp_in_del, where the same SNP stays on the
# base contig: LV=1, CH=0.
is $(grep -v ^# nested_snp_in_ins.vcf | grep -c "LV=0") 2 "both sites are top-level on their own contig"
is $(grep -v ^# nested_snp_in_ins.vcf | grep -c "CH=1") 1 "the nested allele is one contig hop from the base reference"
# The LV=0 record on the gref contig must keep PS: vcfbub's rescue of the children of
# popped bubbles is keyed on it, and it is the only in-VCF link back to the base site.
is $(grep -v ^# nested_snp_in_ins.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=1/ && $8 ~ /PS=/' | wc -l) 1 "a per-contig top-level site on a gref contig still carries PS"

# The cover selects three contigs here (gref_x, gref_x_1_alt, gref_x_2_alt) but
# gref_x_2_alt carries no variant, and a reference contig that produced no record is not
# worth declaring.  So the header has two.
is $(grep -c "^##contig" nested_snp_in_ins.vcf) 2 "contigs with no records are not declared"
is $(grep -c "^##contig=<ID=gref_x_2_alt," nested_snp_in_ins.vcf) 0 "the undeclared contig is the one with no records"
is $(grep -v "^#" nested_snp_in_ins.vcf | cut -f1 | sort -u | wc -l) 2 "every declared contig has records, and every contig with records is declared"
# The contig set comes from the records, not the snarl tree, so it does not need -a.
vg deconstruct nested_snp_in_ins.gref.pg -P gref_x > nested_snp_in_ins.flat.vcf
is $(grep -c "^##contig" nested_snp_in_ins.flat.vcf) 1 "pruning works without -a"
rm -f nested_snp_in_ins.flat.vcf

rm -f nested_snp_in_ins.gref.pg nested_snp_in_ins.vcf nested_snp_in_ins.tsv nested_snp_in_ins_truth.tsv nested_snp_in_ins_contigs.tsv

# Test: Double-nested SNP (run with -t 1, -t 2, and default threads to verify determinism)
for thread_opt in "-t 1" "-t 2" ""; do
    thread_label=${thread_opt:- default}
    vg paths --compute-gref --min-gref-len 0 $thread_opt -x nesting/nested_snp_in_nested_ins.gfa -Q x > nested_snp_in_nested_ins.gref.pg
    vg deconstruct nested_snp_in_nested_ins.gref.pg -P gref_x -a > nested_snp_in_nested_ins.vcf
    is $(grep -v ^# nested_snp_in_nested_ins.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/' | wc -l) 1 "one absolutely top-level site in double-nested SNP ($thread_label)"
    is $(grep -v ^# nested_snp_in_nested_ins.vcf | grep "LV=" | wc -l) 3 "all nested sites found in double-nested SNP ($thread_label)"
    is "$(grep -v ^# nested_snp_in_nested_ins.vcf | grep -o 'LV=[0-9]*' | tr '\n' ' ')" "LV=0 LV=0 LV=1 " "per-contig levels in double-nested SNP ($thread_label)"
    rm -f nested_snp_in_nested_ins.gref.pg nested_snp_in_nested_ins.vcf
done

# Test: Nested site with cycle
vg paths --compute-gref --min-gref-len 0 -x nesting/nested_snp_in_ins_cycle.gfa -Q x > nested_snp_in_ins_cycle.gref.pg
vg deconstruct nested_snp_in_ins_cycle.gref.pg -P gref_x -a > nested_snp_in_ins_cycle.vcf
is $(grep -v ^# nested_snp_in_ins_cycle.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/' | wc -l) 1 "one absolutely top-level site found in nested cycle"
is $(grep -v ^# nested_snp_in_ins_cycle.vcf | grep -c "CH=1") 1 "the nested site is one contig hop out, in nested cycle"
rm -f nested_snp_in_ins_cycle.gref.pg nested_snp_in_ins_cycle.vcf

# Test: MNP handling
vg paths --compute-gref --min-gref-len 0 -x nesting/mnp.gfa -Q x > mnp.gref.pg
vg deconstruct mnp.gref.pg -P gref_x -a > mnp.vcf
printf "gref_x\t3\t>2>7\tTCAT\tATTT\n" > mnp_truth.tsv
grep -v ^# mnp.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > mnp.tsv
diff  mnp_truth.tsv mnp.tsv
is "$?" 0 "nested deconstruction handles mnp"

rm -f mnp.gref.pg mnp.vcf mnp_truth.tsv mnp.tsv

# Test 1: Deep nesting (3+ levels) - triple nested SNP
vg paths --compute-gref --min-gref-len 0 -x nesting/triple_nested.gfa -Q x > triple_nested.gref.pg
vg deconstruct triple_nested.gref.pg -P gref_x -a > triple_nested.vcf
is $(grep -v ^# triple_nested.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/' | wc -l) 1 "one absolutely top-level site in triple nested"
is $(grep -v ^# triple_nested.vcf | grep "LV=" | wc -l) 5 "all nested sites found in triple nested"
# One hop into the insertion, then four levels of nesting all on that one gref contig.
is "$(grep -v ^# triple_nested.vcf | grep -o 'LV=[0-9]*' | tr '\n' ' ')" "LV=0 LV=0 LV=1 LV=2 LV=3 " "LV is the per-contig level in triple nested"
is "$(grep -v ^# triple_nested.vcf | grep -o 'CH=[0-9]*' | tr '\n' ' ')" "CH=0 CH=1 CH=1 CH=1 CH=1 " "CH counts steps into non-reference sequence"
rm -f triple_nested.gref.pg triple_nested.vcf

# Test 2: Multiple children at same level - insertion with 2 nested SNPs
vg paths --compute-gref --min-gref-len 0 -x nesting/insertion_with_three_snps.gfa -Q x > insertion_with_three_snps.gref.pg
vg deconstruct insertion_with_three_snps.gref.pg -P gref_x -a > multi_child.vcf
is $(grep -v ^# multi_child.vcf | grep "LV=" | wc -l) 3 "expected number of sites with LV field"
is $(grep -v ^# multi_child.vcf | grep LV=1 | wc -l) 2 "two child SNPs found at level 1"
rm -f insertion_with_three_snps.gref.pg multi_child.vcf

# Test 3: NestingInfo field propagation - verify LV field
vg paths --compute-gref --min-gref-len 0 -x nesting/nested_snp_in_ins.gfa -Q x > field_check.gref.pg
vg deconstruct field_check.gref.pg -P gref_x -a > field_check.vcf
is $(grep -v ^# field_check.vcf | awk -F'\t' '$8 ~ /LV=0/ && $8 ~ /CH=0/' | wc -l) 1 "LV/CH present and zero for the top-level site"
rm -f field_check.gref.pg field_check.vcf

# Test 4: Multiple reference traversals with nesting (should reduce to one)
vg paths --compute-gref --min-gref-len 0 -x nesting/cyclic_ref_nested.gfa -Q x > cyclic_ref_nested.gref.pg
vg deconstruct cyclic_ref_nested.gref.pg -P gref_x -a > cyclic_ref_nested.vcf
is $(grep -v ^# cyclic_ref_nested.vcf | wc -l) 1 "cyclic reference with nesting produces single variant"
rm -f cyclic_ref_nested.gref.pg cyclic_ref_nested.vcf

# Test 5: Cyclic reference outputs multiple variants (one per reference traversal)
# Reference x visits the snarl twice (via node 2 then 3), alt a visits twice (via node 3 then 4)
# With -c 0 to disable context-jaccard (which doesn't work well on tiny graphs), we should get 2 variants
vg deconstruct nesting/cyclic_ref_multiple_variants.gfa -p x -a -c 0 > cyclic_ref_multi.vcf
is $(grep -v ^# cyclic_ref_multi.vcf | wc -l) 2 "cyclic reference with -a outputs variant for each reference traversal"
# Both records share the snarl ID, because the ID names the snarl and not the traversal of
# it.  Each must still report its own reference interval; a name-keyed lookup used to hand
# both of them whichever one was written last.
is "$(grep -v ^# cyclic_ref_multi.vcf | cut -f3 | sort -u | wc -l)" "1" "the two records share one snarl ID"
is "$(grep -v ^# cyclic_ref_multi.vcf | grep -o 'RS=[0-9]*' | tr '\n' ' ')" "RS=20 RS=44 " "each record keeps its own RS, not the other's"
rm -f cyclic_ref_multi.vcf

# A gref cover writes the reference twice: under its own name and in the gref
# namespace.  Whichever of the two you deconstruct against, the other one is not a
# sample -- it is the same sequence, and genotyping it inflates AC/AF/AN/NS.  A second
# reference assembly (CHM13 here) has no gref copy and must survive in both.
vg paths --compute-gref --min-gref-len 1 -x nesting/base_and_gref.gfa -Q GRCh38 > base_and_gref.pg
vg deconstruct base_and_gref.pg -P GRCh38 -a > base_ref.vcf
vg deconstruct base_and_gref.pg -P gref_GRCh38 -a > gref_ref.vcf

is $(grep "^#CHROM" base_ref.vcf | cut -f10- | tr '\t' '\n' | grep -c "^gref_") 0 "base-reference vcf has no gref sample columns"
is $(grep "^#CHROM" gref_ref.vcf | cut -f10- | tr '\t' '\n' | grep -c "^GRCh38") 0 "gref vcf has no base-reference sample column"
is $(grep "^#CHROM" base_ref.vcf | cut -f10- | wc -w) 3 "base-reference vcf keeps every genuine sample"
is $(grep "^#CHROM" gref_ref.vcf | cut -f10- | wc -w) 3 "gref vcf keeps every genuine sample"
is $(grep -v "^#" base_ref.vcf | grep -c "AN=3") 1 "base-reference vcf allele counts are not inflated by the gref copy"
is $(grep -v "^#" gref_ref.vcf | grep -c "AN=3") 2 "gref vcf allele counts are not inflated by the base reference"

rm -f base_and_gref.pg base_ref.vcf gref_ref.vcf

# A subranged reference must survive into the gref VCF with its coordinates: every region
# deconstructed, at the same POS as the base VCF.
vg paths --compute-gref --min-gref-len 1 -x nesting/subranged_ref.gfa -Q GRCh38 > subranged.pg
vg deconstruct subranged.pg -P GRCh38 -a | grep -v "^#" | cut -f2,4,5 > subranged_base.tsv
vg deconstruct subranged.pg -P gref_GRCh38 -a | grep -v "^#" | cut -f2,4,5 > subranged_gref.tsv
is $(cat subranged_gref.tsv | wc -l) 2 "every subpath of a subranged reference is deconstructed in the gref vcf"
diff subranged_base.tsv subranged_gref.tsv
is $? 0 "gref vcf keeps the subrange offsets, so positions match the base vcf"

rm -f subranged.pg subranged_base.tsv subranged_gref.tsv

# Selecting one contig's gref reference must not let the base sample back in through its
# other contigs: it would contribute an all-reference column and inflate AN.
vg paths --compute-gref --min-gref-len 1 -x nesting/two_contig_gref.gfa -Q GRCh38 > two_contig.pg
vg deconstruct two_contig.pg -P gref_GRCh38#0#chr1 -a > one_contig.vcf
is $(grep "^#CHROM" one_contig.vcf | cut -f10- | tr '\t' '\n' | grep -c "^GRCh38$") 0 "base reference is not a sample when only its gref contig is selected"
is $(grep -v "^#" one_contig.vcf | grep -c "AN=2") 1 "allele counts are not inflated when only one gref contig is selected"

rm -f two_contig.pg one_contig.vcf

# With no -p/-P, every reference-sense path is a reference, including both views of the
# gref pair.  The record belongs on the base contig: a gref name sorts before the path it
# was copied from (gref_x < x), so name order alone would put the derived name on it.
vg paths --compute-gref --min-gref-len 1 -x nesting/nested_snp_in_ins.gfa -Q x > default_ref.pg
vg deconstruct default_ref.pg -a > default_ref.vcf
# LV=0 alone is not enough: under per-contig LV several records are top-level.  The one
# that is top-level in the whole file is the one that also took no contig hop.
is $(grep -v "^#" default_ref.vcf | awk '$8 ~ /LV=0/ && $8 ~ /CH=0/ {print $1}') "x" "top-level record goes on the base contig, not its gref copy"
is $(grep -v "^#" default_ref.vcf | grep -c "RC=x;") 2 "nested records point back at the base contig"

rm -f default_ref.pg default_ref.vcf

# The base/gref link has to be recognised even when the reference is haplotype sense --
# any GFA without an RS header.  Going from a gref name back to the base path by dropping
# the prefix cannot work there (the phase block is not recoverable), which left the base
# sample in the VCF and inflated the counts.
vg paths --compute-gref --min-gref-len 1 -x nesting/two_contig_gref_nors.gfa -Q GRCh38#0#chr1 > nors.pg
vg deconstruct nors.pg -P gref_GRCh38#0#chr1 -a > nors.vcf
is $(grep "^#CHROM" nors.vcf | cut -f10- | tr '\t' '\n' | grep -c "^GRCh38$") 0 "base reference is excluded even when it is haplotype sense"

# ... and the sample set must not depend on whether the haplotypes come from a GBWT
vg gbwt -E -x nors.pg -o nors.gbwt
diff <(vg deconstruct nors.pg -P gref_GRCh38#0#chr1 -a | grep "^#CHROM") \
     <(vg deconstruct nors.pg -g nors.gbwt -P gref_GRCh38#0#chr1 -a | grep "^#CHROM")
is $? 0 "the gbwt sample scan excludes the other view of the reference too"

rm -f nors.pg nors.vcf nors.gbwt

# =============================================================================
# RC, RS, RD tag tests (reference coordinate tags for nested snarls)
# =============================================================================

# Test: RC, RS, RD headers and tags in deconstruct output
vg paths --compute-gref --min-gref-len 0 -x nesting/triple_nested.gfa -Q x > rc_decon_test.gref.pg
vg deconstruct rc_decon_test.gref.pg -P gref_x -a > rc_decon_test.vcf

# A gref fragment whose enclosing snarl only the reference and its own gref copy span.  That
# parent has no variant, so it never reaches the VCF, and the walk that looks for an ancestor to
# take RC/RS/RD from finds nothing.  It used to fall back to this record's own contig and
# position, which on a gref contig is not a reference coordinate at all -- and is indistinguishable
# from the legitimate case where a record really is top-level on a reference contig.
vg paths --compute-gref --min-gref-len 1 -x nesting/gref_island_no_parent_record.gfa -Q REF > island.pg
vg deconstruct island.pg -P gref_REF -a -C > island.vcf
is $(grep -v "^#" island.vcf | awk '$8 ~ /LV=0/ && $8 ~ /CH=0/ {print $1}') "chr1_1_alt" "the island record has no ancestor in the VCF"
is $(grep -v "^#" island.vcf | grep -c "RC=chr1;") 1 "a suppressed parent still gives its child a reference coordinate"
is $(grep -v "^#" island.vcf | grep -c "RC=chr1_1_alt") 0 "no record names its own gref contig as its reference"
rm -f island.pg island.vcf

# Check for RC, RS, RD headers
is $(grep -c "##INFO=<ID=RC" rc_decon_test.vcf) 1 "deconstruct: RC header is present in VCF"
is $(grep -c "##INFO=<ID=RS" rc_decon_test.vcf) 1 "deconstruct: RS header is present in VCF"
is $(grep -c "##INFO=<ID=RD" rc_decon_test.vcf) 1 "deconstruct: RD header is present in VCF"

# Check that all variants have RC, RS, RD tags
# Use grep -c "" instead of wc -l: macOS wc -l pads output with leading
# whitespace which breaks string comparison in is assertions.
RC_DECON_COUNT=$(grep -v "^#" rc_decon_test.vcf | grep -c "")
RC_DECON_TAG=$(grep -v "^#" rc_decon_test.vcf | grep -c "RC=")
is "$RC_DECON_TAG" "$RC_DECON_COUNT" "deconstruct: all variants have RC tag"

RS_DECON_TAG=$(grep -v "^#" rc_decon_test.vcf | grep -c "RS=")
is "$RS_DECON_TAG" "$RC_DECON_COUNT" "deconstruct: all variants have RS tag"

RD_DECON_TAG=$(grep -v "^#" rc_decon_test.vcf | grep -c "RD=")
is "$RD_DECON_TAG" "$RC_DECON_COUNT" "deconstruct: all variants have RD tag"

# Check that nested variants point to top-level's contig
NESTED_DECON_RC=$(grep -v "^#" rc_decon_test.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | grep -c "RC=gref_x")
NESTED_DECON_COUNT=$(grep -v "^#" rc_decon_test.vcf | awk -F'\t' '$8 ~ /LV=[1-9]/' | grep -c "")
is "$NESTED_DECON_RC" "$NESTED_DECON_COUNT" "deconstruct: all nested variants have RC=gref_x (top-level contig)"

rm -f rc_decon_test.gref.pg rc_decon_test.vcf

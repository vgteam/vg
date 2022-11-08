#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 24

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > tiny.vg
vg index tiny.vg -x tiny.xg
vg deconstruct tiny.xg -p x -t 1 > tiny_decon.vcf
# we pop out that GC allele because it gets invented by the adjacent snps in the graph
gzip -dc tiny/tiny.vcf.gz | tail -3 | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > tiny_orig.tsv
cat tiny_decon.vcf | grep -v "#" | grep -v GC | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > tiny_dec.tsv
diff tiny_orig.tsv tiny_dec.tsv
is "$?" 0 "deconstruct retrieved original VCF (modulo adjacent snp allele)"
grep '>1>6' tiny_decon.vcf | awk '{print $8}' > allele_travs
printf "AT=>1>3>5>6,>1>3>4>6,>1>2>5>6,>1>2>4>6\n" > expected_allele_travs
diff allele_travs expected_allele_travs
is "$?" 0 "deconstruct produced expected AT field"

rm -f tiny.vg tiny.xg tiny_decon.vcf tiny_orig.tsv tiny_dec.tsv allele_travs expected_allele_travs

vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 | vg mod -U 10 - | vg mod -c - > hla.vg
vg index hla.vg -x hla.xg

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


rm -f hla_decon.vcf hla_decon_path.vcf  hla_decon.tsv hla_decon_path.tsv hla.vg hla.xg

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

printf "P\ty\t1+,2-,3+\t*\n" >> inv.gfa
vg gbwt -G inv.gfa -g inv.chop.gbz --gbz-format --max-node 5
vg deconstruct inv.gbz -p x > inv.chop.decon
vg gbwt -G inv.gfa -g inv.gbz --gbz-format --max-node 1025
vg deconstruct inv.gbz -p x > inv.decon
diff inv.decon inv.chop.decon
is "$?" 0 "deconstruct automatically applies translation from gbz"

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
vg index cycle.vg -x cycle.xg
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
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg
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



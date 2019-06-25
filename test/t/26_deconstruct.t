#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 16

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > tiny.vg
vg deconstruct tiny.vg -p x -t 1 > tiny_decon.vcf
# we pop out that GC allele because it gets invented by the adjacent snps in the graph
gzip -dc tiny/tiny.vcf.gz | tail -3 | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > tiny_orig.tsv
cat tiny_decon.vcf | grep -v "#" | grep -v GC | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > tiny_dec.tsv
diff tiny_orig.tsv tiny_dec.tsv
is "$?" 0 "deconstruct retrieved original VCF (modulo adjacent snp allele)"

rm -f tiny.vg tiny_decon.vcf tiny_orig.tsv tiny_dec.tsv

vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 | vg mod -U 10 - | vg mod -c - > hla.vg

vg deconstruct hla.vg -p "gi|157734152:29563108-29564082" > hla_decon.vcf
is $(grep -v "#" hla_decon.vcf | wc -l) 17 "deconstructed hla vcf has correct number of sites"
is $(grep -v "#" hla_decon.vcf | grep 822 | awk '{print $4 "-" $5}') "C-CGCGGGCGCCGTGGATGGAGCA" "deconstructed hla vcf has correct insertion"
vg deconstruct hla.vg -p "gi|568815592:29791752-29792749" > hla_decon.vcf
is $(grep -v "#" hla_decon.vcf | wc -l) 17 "deconstructed hla vcf with other path has correct number of sites"
is $(grep -v "#" hla_decon.vcf | grep 824 | awk '{print $4 "-" $5}') "CGCGGGCGCCGTGGATGGAGCA-C" "deconstructed hla vcf has correct deletion"

vg deconstruct hla.vg -p "gi|568815592:29791752-29792749" -e > hla_decon_path.vcf
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


rm -f hla_decon.vcf hla_decon_path.vcf  hla_decon.tsv hla_decon_path.tsv hla.vg

cp sv/x.inv.gfa inv.gfa
printf "P\ty\t1+,2-,3+\t9M,20M,21M\n" >> inv.gfa
vg view -Fv inv.gfa > inv.vg
vg deconstruct inv.vg -p x -e > inv_decon.vcf
grep -v "#" inv_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > inv_decon.tsv
printf "x\t10\tCTTGGAAATTTTCTGGAGTT\tAACTCCAGAAAATTTCCAAG\t1\n" > inv_truth.tsv
diff inv_decon.tsv inv_truth.tsv
is "$?" 0 "deconstruct correctly handles a simple inversion"

rm -f inv_decon.vcf inv_decon.tsv inv_truth.tsv

vg deconstruct inv.vg -p y -e > inv_decon.vcf
grep -v "#" inv_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > inv_decon.tsv
printf "y\t10\tAACTCCAGAAAATTTCCAAG\tCTTGGAAATTTTCTGGAGTT\t1\n" > inv_truth.tsv
diff inv_decon.tsv inv_truth.tsv
is "$?" 0 "deconstruct correctly handles a simple inversion when the reference contains the reversing edge"

rm -f inv.gfa inv.vg inv_decon.vcf inv_decon.tsv inv_truth.tsv


vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view -g - > cyclic_tiny.gfa
printf "L\t12\t+\t9\t+\t0M\n" >> cyclic_tiny.gfa
printf "P\ty\t1+,3+,5+,6+,8+,9+,11+,12+,9+,10+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,19M,1M,4M,1M,11M\n" >> cyclic_tiny.gfa
vg view -Fv cyclic_tiny.gfa > cyclic_tiny.vg
vg deconstruct cyclic_tiny.vg -p y -e > cyclic_tiny_decon.vcf
grep -v "#" cyclic_tiny_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > cyclic_tiny_decon.tsv
printf "y\t13\tGGAAATTTTCTGGAGTTCTATTATATAAATTTTCTGGAGTTCTATAATATT\tGGAAATTTTCTGGAGTTCTATTATATT\t1\n" > cyclic_tiny_truth.tsv
diff cyclic_tiny_decon.tsv cyclic_tiny_truth.tsv
is "$?" 0 "deconstruct correctly handles a cycle in the reference path"

rm -f cyclic_tiny_decon.vcf cyclic_tiny_decon.tsv cyclic_tiny_truth.tsv

vg deconstruct cyclic_tiny.vg -p x -e > cyclic_tiny_decon.vcf
grep -v "#" cyclic_tiny_decon.vcf | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $10}' > cyclic_tiny_decon.tsv
printf "x\t13\tGGAAATTTTCTGGAGTTCTATTATATT\tGGAAATTTTCTGGAGTTCTATTATATAAATTTTCTGGAGTTCTATAATATT\t1\n" > cyclic_tiny_truth.tsv
diff cyclic_tiny_decon.tsv cyclic_tiny_truth.tsv
is "$?" 0 "deconstruct correctly handles a cycle in the alt path"

rm -f cyclic_tiny_decon.vcf cyclic_tiny_decon.tsv cyclic_tiny_truth.tsv cyclic_tiny.gfa cyclic_tiny.vg


vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa | vg view -g - > tiny_names.gfa
printf "P\tref.1\t1+,3+,5+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
printf "P\talt1.1\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
printf "P\talt1.2\t1+,2+,4+,6+,7+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
printf "P\talt2.3\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
printf "P\talt2.3\t1+,2+,4+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M\n" >> tiny_names.gfa
vg view -Fv tiny_names.gfa > tiny_names.vg
vg deconstruct tiny_names.vg -P ref -A alt1,alt2 -e > tiny_names_decon.vcf
is $(grep -v "#" tiny_names_decon.vcf | wc -l) 2 "-P -A options return correct number of variants"
is $(grep -v "#" tiny_names_decon.vcf | grep ref.1 | wc -l) 2 "-P -A options use correct reference name"
is $(grep -v "#" tiny_names_decon.vcf | grep ref.1 | grep 14 | grep "CONFLICT=alt1" | wc -l) 1 "-P -A identifies conflict in alt1 in second variant"




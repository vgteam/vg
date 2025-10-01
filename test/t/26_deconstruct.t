#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 40

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
vg deconstruct small_cluster.gfa -p x -L 0.3 > small_cluster_3.vcf
is "$(tail -1 small_cluster_0.vcf | awk '{print $5}')" "GATTTGA,G" "cluster-free deconstruction finds all alt alleles"
is "$(tail -1 small_cluster_3.vcf | awk '{print $5}')" "G" "clustered deconstruction finds fewer alt alleles"
is "$(tail -1 small_cluster_3.vcf | awk '{print $10}')" "0:0.333:0" "clustered deconstruction finds correct allele info"

rm -f small_cluster.gfa small_cluster_0.vcf small_cluster_3.vcf

vg deconstruct nesting/nested_snp_in_del.gfa -p x -nR > nested_snp_in_del.vcf
grep -v ^# nested_snp_in_del.vcf | awk '{print $4 "\t" $5 "\t" $10}' > nested_snp_in_del.tsv
printf "CATG\tCAAG,C\t1|2\n" > nested_snp_in_del_truth.tsv
printf "T\tA,*\t1|2\n" >> nested_snp_in_del_truth.tsv
diff nested_snp_in_del.tsv nested_snp_in_del_truth.tsv
is "$?" 0 "nested deconstruction gets correct star-allele for snp inside deletion"

is $(grep PA=0 nested_snp_in_del.vcf | wc -l) 2 "PA=0 correctly set at both sites of neseted deletion"

rm -f nested_snp_in_del.vcf nested_snp_in_del.tsv nested_snp_in_del_truth.tsv

vg deconstruct nesting/nested_snp_in_del.gfa -p x -n > nested_snp_in_del.vcf
grep -v ^# nested_snp_in_del.vcf | awk '{print $4 "\t" $5 "\t" $10}' > nested_snp_in_del.tsv
printf "CATG\tCAAG,C\t1|2\n" > nested_snp_in_del_truth.tsv
printf "T\tA\t1|.\n" >> nested_snp_in_del_truth.tsv
diff nested_snp_in_del.tsv nested_snp_in_del_truth.tsv
is "$?" 0 "nested deconstruction gets correct allele for snp inside deletion (without -R)"

rm -f nested_snp_in_del.vcf nested_snp_in_del.tsv nested_snp_in_del_truth.tsv

vg deconstruct nesting/nested_snp_in_ins.gfa -p x -n > nested_snp_in_ins.vcf
grep -v ^# nested_snp_in_ins.vcf | awk '{print $4 "\t" $5 "\t" $10}' > nested_snp_in_ins.tsv
printf "A\tT\t0|1\n" > nested_snp_in_ins_truth.tsv
printf "C\tCAAG,CATG\t1|2\n" >> nested_snp_in_ins_truth.tsv
diff nested_snp_in_ins.tsv nested_snp_in_ins_truth.tsv
is "$?" 0 "nested deconstruction gets correct allele for snp inside insert"

is $(grep LV=0 nested_snp_in_ins.vcf | head -1 | grep PA=0 | wc -l) 1 "PA=0 set for base allele of nested insertion"
is $(grep LV=1 nested_snp_in_ins.vcf | tail -1 | grep PA=1 | wc -l) 1 "PA=1 set for nested allele of nested insertion"

grep ^##contig nested_snp_in_ins.vcf > nested_snp_in_ins_contigs.tsv
printf "##contig=<ID=a#1#y0,length=5>\n" > nested_snp_in_ins_contigs_truth.tsv
printf "##contig=<ID=x,length=2>\n" >> nested_snp_in_ins_contigs_truth.tsv
diff nested_snp_in_ins_contigs.tsv nested_snp_in_ins_contigs_truth.tsv
is "$?" 0 "nested deconstruction gets correct contigs in vcf header for snp inside insert"

rm -f nested_snp_in_ins.tsv nested_snp_in_ins.tsv nested_snp_in_ins_truth.tsv  nested_snp_in_ins_contigs.tsv nested_snp_in_ins_contigs_truth.tsv

vg deconstruct nesting/nested_snp_in_ins2.gfa -p x -nR | grep -v ^# | awk '{print $4 "\t" $5 "\t" $10 "\t" $11}' > nested_snp_in_ins2.tsv
printf "A\tT,*\t0|1\t0|2\n" > nested_snp_in_ins2._truth.tsv
printf "C\tCAAG,CATG\t1|2\t1|0\n" >> nested_snp_in_ins2._truth.tsv
diff nested_snp_in_ins2.tsv nested_snp_in_ins2._truth.tsv
is "$?" 0 "nested deconstruction gets correct star allele for snp ins2.ide ins2.ert"

rm -f nested_snp_in_ins2.tsv nested_snp_in_ins2._truth.tsv

# todo: the integrated snarl finder doesnt anchor to the reference
#       probably not an issue on most real graphs from vg construct / minigraph cacuts
#       but seems like something that needs reviewing
vg snarls nesting/nested_snp_in_nested_ins.gfa -A cactus > nested_snp_in_nested_ins.snarls
vg deconstruct nesting/nested_snp_in_nested_ins.gfa -r nested_snp_in_nested_ins.snarls -P x -n > nested_snp_in_nested_ins.vcf
is $(grep -v ^# nested_snp_in_nested_ins.vcf | grep LV=0 | awk '{print $8}') "AC=1,1;AF=0.5,0.5;AN=2;AT=>1>6,>1>2>3>31>33>34>35>5>6,>1>2>3>31>32>34>35>5>6;NS=1;PA=0;PL=1;PR=1;RC=x;RD=1;RL=1;RS=1;LV=0" "INFO tags correct for level-0 site of double-nested SNP"
is $(grep -v ^# nested_snp_in_nested_ins.vcf | grep LV=1 | awk '{print $8}') "AC=1;AF=0.5;AN=2;AT=>2>3>31>33>34>35>5,>2>3>31>32>34>35>5;NS=1;PA=1;PL=8;PR=1;RC=x;RD=1;RL=8;RS=1;LV=1;PS=>1>6" "INFO tags correct for level-1 site of double-nested SNP"
is $(grep -v ^# nested_snp_in_nested_ins.vcf | grep LV=2 | awk '{print $8}') "AC=1;AF=0.5;AN=2;AT=>31>33>34,>31>32>34;NS=1;PA=0;PL=5;PR=5;RC=x;RD=1;RL=8;RS=1;LV=2;PS=>2>5" "INFO tags correct for level-2 site of double-nested SNP"

rm -f nested_snp_in_nested_ins.snarls nested_snp_in_nested_ins.vcf

vg deconstruct nesting/nested_snp_in_ins_cycle.gfa -P x -n > nested_snp_in_ins_cycle.vcf
printf "A#1#y0\t3\t>2>5\tA\tT\tAC=2;AF=0.5;AN=4;AT=>2>4>5,>2>3>5;NS=2;PA=1;PL=7;PR=1;RC=x;RD=1;RL=7;RS=1;LV=1;PS=>1>6\t0|1\t0|1\n" >  nested_snp_in_ins_cycle_truth.tsv
printf "x\t1\t>1>6\tC\tCAAGAAG,CATG,CAAG\tAC=1,2,1;AF=0.25,0.5,0.25;AN=4;AT=>1>6,>1>2>4>5>2>4>5>6,>1>2>3>5>6,>1>2>4>5>6;NS=2;PA=0;PL=1;PR=1;RC=x;RD=1;RL=1;RS=1;LV=0\t1|2\t3|2\n" >> nested_snp_in_ins_cycle_truth.tsv
grep -v ^#  nested_snp_in_ins_cycle.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $10 "\t" $11}' > nested_snp_in_ins_cycle.tsv
diff nested_snp_in_ins_cycle.tsv nested_snp_in_ins_cycle_truth.tsv
is "$?" 0 "nested deconstruction handles cycle"

rm -f nested_snp_in_ins_cycle.vcf nested_snp_in_ins_cycle_truth.tsv nested_snp_in_ins_cycle_truth.tsv

vg snarls nesting/mnp.gfa --algorithm cactus > mnp.snarls
vg deconstruct nesting/mnp.gfa -r mnp.snarls -p x -n -f mnp.fa > mnp.vcf
printf "x\t3\t>2>7\tTCAT\tATTT\n" > mnp_truth.tsv
grep -v ^# mnp.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > mnp.tsv
diff  mnp_truth.tsv mnp.tsv
is "$?" 0 "nested deconstruction handles mnp"
printf "a#1#y0#0\t2\t6\t>8>11\tx\t2\t6\n" > mnp.nesting.truth.tsv
diff mnp.fa.nesting.tsv mnp.nesting.truth.tsv
is "$?" 0 "nested deconstruction makes correct mnp tsv"
printf ">a#1#y0#0[2-6]\nATTT\n"  > mnp.fa.truth
diff mnp.fa mnp.fa.truth
is "$?" 0 "nested deconstruction makes correct fasta"

rm -f mnp.snarls  mnp.vcf mnp_truth.tsv mnp.tsv mnp.nesting.truth.tsv  mnp.fa.truth





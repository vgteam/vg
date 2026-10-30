#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order

plan tests 53

# A triallelic SNP with the haplotypes stored as W lines at offset 1000 of contig c,
# as Minigraph-Cactus stores haplotype fragments
awk -v OFS='\t' '$1 == "P" && $2 ~ /^hap/ {
    n = split($3, steps, ",")
    walk = ""
    for (i = 1; i <= n; i++) {
        walk = walk (substr(steps[i], length(steps[i])) == "+" ? ">" : "<") substr(steps[i], 1, length(steps[i]) - 1)
    }
    print "W", "s", ($2 == "hapA" ? 1 : 2), "c", 1000, 1201, walk
    next
} { print }' nesting/core_snp_in_flanks.gfa > lift.gfa
vg paths -x lift.gfa -Q x --compute-gref > lift.gref.vg
vg convert -f lift.gref.vg > lift.gref.gfa
vg gbwt -G lift.gref.gfa --gbz-format -g lift.gbz
vg unitigs -r -k 11 -g lift.tigs.gaf lift.gbz > lift.tigs.fa

# The same 5 bases on the gref and on both junction tigs, a stranded feature across the
# SNP, a feature with no BED name, and a duplicate line
printf 'gref_x\t92\t97\tg5\t0\t+\nu1\t2\t7\tg5\t0\t+\nu2\t2\t7\tg5\t0\t+\n' > lift.bed
printf 'u1\t8\t13\tsnpA\t5\t-\nu2\t8\t13\tsnpB\nu2\t8\t13\tsnpB\ngref_x\t0\t5\n' >> lift.bed
vg lift -u lift.tigs.gaf -a lift.asm.bed lift.gbz lift.bed > lift.gaf
is $? 0 "vg lift runs on a cover of gref paths and junction tigs"

is "$(awk 'END { print NR }' lift.gaf)" "4" "copies of a feature on several cover sequences are kept once"

is "$(grep -c $'\tsp:i:3\tgr:Z:gref_x:93-97:+:5=\tpc:Z:gref_x:93-97:+\tan:Z:g5$' lift.gaf)" "1" "a feature on the gref is placed at its gref coordinates, which are also its path coordinate"

is "$(grep -c $'<6<4<2\t.*\tgr:Z:gref_x:99-103:-:2=1X2=\tpc:Z:s#1#c:1099-1103:-\tan:Z:snpA$' lift.gaf)" "1" "a feature on the - strand of a tig is lifted to the reverse walk, with the SNP as a substitution"

is "$(awk '$4 == "g5"' lift.asm.bed | cut -f 1-3,6 | sort | tr '\t\n' ' ;')" "s#1#c 1092 1097 +;s#2#c 1092 1097 +;x 92 97 +;" "a haplotype place is on its contig, past the fragment's offset"

is "$(awk '$4 == "snpA"' lift.asm.bed | cut -f 1-6 | tr '\t' ' ')" "s#1#c 1098 1103 snpA 5 -" "a stranded feature keeps its fields and gets the strand of the path"

vg paths -x lift.gbz -F -Q 's#1#c' | awk 'NR > 1' | tr -d '\n' | cut -c 99-103 > lift.place.txt
awk '$1 == ">u1" { getline; print substr($0, 9, 5) }' lift.tigs.fa > lift.feature.txt
is "$(cat lift.place.txt) $(awk '{ print length($0) }' lift.place.txt)" "$(cat lift.feature.txt) 5" "the bases at a place on a haplotype are those of the feature on the tig"

is "$(awk '$4 == "." && $5 == 0 && $6 == "."' lift.asm.bed | awk 'END { print NR }')" "3" "a feature without a name, score or strand gets ., 0 and ."

join -t $'\t' <(cut -f 1,13 lift.gaf | sed 's/sp:i://' | sort) <(awk -v OFS='\t' '{ print $NF }' lift.asm.bed | sort | uniq -c | awk -v OFS='\t' '{ print $2, $1 }') > lift.support.txt
is "$(awk -F '\t' '$2 == $3 { n++ } END { print n + 0 }' lift.support.txt) $(cut -f 1 lift.gaf | sort | uniq -d | awk 'END { print NR }')" "4 0" "each feature has a unique ID and its support counts its places on paths"

vg lift -u lift.tigs.gaf -a lift.s1.bed -Q 's#1' lift.gbz lift.bed > /dev/null
is "$(cut -f 1 lift.s1.bed | sort -u)" "s#1#c" "places can be limited to paths by name prefix"

vg lift -u lift.tigs.gaf -t 1 -a lift.t1.bed lift.gbz lift.bed > lift.t1.gaf
vg lift -u lift.tigs.gaf -t 4 -a lift.t4.bed lift.gbz lift.bed > lift.t4.gaf
cmp -s lift.t1.gaf lift.t4.gaf && cmp -s lift.t1.bed lift.t4.bed
is $? 0 "the output does not depend on the thread count"

# Path coordinates: a feature on the gref and on tigs, features off the gref on tigs, and
# features on a path, along the gref and off it
cp lift.bed lift.names.input.bed
printf 'u2\t9\t11\tcg\nx\t95\t105\talonggref\ns#1#c\t1098\t1103\toffgref\n' >> lift.names.input.bed
vg lift -u lift.tigs.gaf -c lift.names.bed lift.gbz lift.names.input.bed > lift.names.gaf 2> lift.names.err
is "$(grep -c $'\tpc:Z:s#2#c:1100-1101:+\tan:Z:cg$' lift.names.gaf)" "1" "a feature off the gref is named on the path that carries it, in contig coordinates"

is "$(grep -c $'\tpc:Z:gref_x:96-105:+\tan:Z:alonggref$' lift.names.gaf)" "1" "a feature read from a path along the gref is named on the gref"

is "$(grep -c $'\tpc:Z:s#1#c:1099-1103:+\tan:Z:offgref$' lift.names.gaf) $(grep -c . lift.names.err)" "1 0" "a feature read from a path off the gref is named on a path that carries it"

is "$(cut -f 1-4 lift.names.bed | tr '\t\n' ' ;')" "gref_x 0 5 .;gref_x 92 97 g5;gref_x 95 105 alonggref;s#1#c 1098 1103 offgref;s#1#c 1098 1103 snpA;s#2#c 1098 1103 snpB;s#2#c 1099 1101 cg;" "the canonical BED has each feature once, at its path coordinate, in priority order"

grep -v '^@' lift.tigs.gaf > lift.tigs.noheader.gaf
vg lift -u lift.tigs.noheader.gaf lift.gbz lift.names.input.bed 2> /dev/null | cmp -s - lift.names.gaf
is $? 0 "a tig GAF without a header gives the same output"

printf 'u99\t0\t5\n' > lift.unknown.bed
vg lift -u lift.tigs.gaf lift.gbz lift.unknown.bed > /dev/null 2>&1
is $? 1 "a feature on an unknown sequence is an error"

printf 'u1\t10\t30\n' > lift.long.bed
vg lift -u lift.tigs.gaf lift.gbz lift.long.bed 2>&1 > /dev/null | grep -q "past the end"
is $? 0 "a feature past the end of its sequence is an error"

vg lift lift.gbz lift.bed 2>&1 > /dev/null | grep -q "use -u"
is $? 0 "a feature on a tig needs the tig GAF"

# A plain cover, without a gref
vg gbwt -G lift.gfa --gbz-format -g lift.plain.gbz
vg unitigs -k 11 -g lift.plain.gaf lift.plain.gbz > lift.plain.fa
printf 'u1\t0\t5\tp\nu3\t9\t14\tq\n' > lift.plain.bed
vg lift -u lift.plain.gaf -a lift.plain.asm.bed lift.plain.gbz lift.plain.bed > lift.plain.out.gaf
is "$(grep -c 'gr:Z:' lift.plain.out.gaf)" "0" "features on a plain cover have no gref placement"

is "$(awk 'END { print NR }' lift.plain.asm.bed)" "4" "features on a plain cover are lifted to every path that visits them"

# Path to path: features on paths, in full-path coordinates
printf 'x\t80\t120\tspan\t0\t+\ns#1#c\t1080\t1120\tfromA\t0\t-\n' > lift.paths.bed
vg lift -a lift.paths.exact.bed lift.gbz lift.paths.bed > /dev/null
is "$(awk '$4 == "span"' lift.paths.exact.bed | cut -f 1-3)" "$(printf 'x\t80\t120')" "a feature across a SNP has only its own exact place"

vg lift -m 0.9 -a lift.paths.lift.bed lift.gbz lift.paths.bed > /dev/null
is "$(awk '$4 == "span" { print $1, $2, $3, $6, $NF }' lift.paths.lift.bed | tr '\n' ';')" "x 80 120 + 1;s#1#c 1080 1120 + 0.975;s#2#c 1080 1120 + 0.975;" "a liftover reaches the haplotypes with a SNP, at their contig coordinates"

is "$(awk '$4 == "fromA" && $1 == "x" { print $2, $3, $6 }' lift.paths.lift.bed)" "80 120 -" "a feature on a haplotype contig lifts to the reference on its strand"

vg lift -m 0.9 -b -a lift.paths.blocks.bed lift.gbz lift.paths.bed > /dev/null
is "$(awk '$4 == "span" && $1 == "s#1#c" { print $7, $8, $10, $11, $12 }' lift.paths.blocks.bed)" "1080 1120 2 20,19, 0,21," "blocks leave out the SNP base"

printf 's#1#c\t10\t20\tbefore\n' > lift.outside.bed
vg lift lift.gbz lift.outside.bed 2>&1 > /dev/null | grep -q "1 features on paths are not within one fragment"
is $? 0 "a feature outside the fragments of its path is skipped with a warning"

# A 50 bp insertion on the haplotypes, after base 70 of x
vg gbwt -G nesting/core_ins50.gfa --gbz-format -g lift.ins.gbz
printf 'x\t50\t110\twide\n' > lift.ins.bed
vg lift -m 0.9 -b -a lift.ins.blocks.bed lift.ins.gbz lift.ins.bed > /dev/null
is "$(awk '$1 == "hapA" { print $2, $3, $10, $11, $12, $NF }' lift.ins.blocks.bed)" "50 160 2 20,40, 0,70, 1" "a liftover spans an insertion shorter than the feature"

printf 'x\t60\t100\tnarrow\n' > lift.ins.bed
vg lift -m 0.7 -a lift.ins.lift.bed lift.ins.gbz lift.ins.bed > /dev/null
is "$(awk '$1 == "hapA" { print $2, $3, $NF }' lift.ins.lift.bed)" "120 150 0.75" "an insertion longer than the feature splits it"

printf 'x\t0\t10\tshared\n' > lift.shared.bed
is "$(vg lift lift.ins.gbz lift.shared.bed | grep -o 'pc:Z:[^[:space:]]*') $(vg lift -P x lift.ins.gbz lift.shared.bed | grep -o 'pc:Z:[^[:space:]]*')" "pc:Z:hapA:1-10:+ pc:Z:x:1-10:+" "a feature carried by several paths is named on the first in priority order"

vg unitigs -k 11 -g lift.ins.tigs.gaf lift.ins.gbz > /dev/null
vg lift -u lift.ins.tigs.gaf lift.gbz lift.bed 2>&1 > /dev/null | grep -q "not compatible"
is $? 0 "tigs of another graph are rejected"

# VCF input: the 50 bp insertion as a gref alt path, with the haplotypes as W lines at
# offset 1000 of contig c. Haplotype 2 has a SNP at the last base of the insertion.
awk -v OFS='\t' '$1 == "P" && $2 ~ /^hap/ {
    n = split($3, steps, ",")
    walk = ""
    for (i = 1; i <= n; i++) {
        walk = walk (substr(steps[i], length(steps[i])) == "+" ? ">" : "<") substr(steps[i], 1, length(steps[i]) - 1)
    }
    print "W", "s", ($2 == "hapA" ? 1 : 2), "c", 1000, 1190, walk
    next
} { print }' nesting/core_ins50.gfa > lift.vcf.gfa
vg paths -x lift.vcf.gfa -Q x --compute-gref > lift.vcf.gref.vg
vg convert -f lift.vcf.gref.vg > lift.vcf.gref.gfa
vg gbwt -G lift.vcf.gref.gfa --gbz-format -g lift.vcf.gbz
vg deconstruct -P x -a lift.vcf.gbz > lift.decon.vcf 2> /dev/null
vg lift -c lift.decon.bed lift.vcf.gbz lift.decon.vcf > lift.decon.lifted.vcf
is $? 0 "vg lift reads the allele walks of a VCF from vg deconstruct"

is "$(grep '^x' lift.decon.lifted.vcf | cut -f 8 | tr ';' '\n' | grep -E '^(GR|PC)=' | tr '\n' ' ')" "GR=gref_x_1_alt:1-50:+:50=,gref_x_1_alt:1-49:+:49=1I PC=gref_x_1_alt:1-50:+,s#2#c:1071-1120:+ " "each ALT allele, its walk without the first and last nodes, is placed and named on the first path that carries it"

is "$(grep -c -E '^##INFO=<ID=(GR|PC),Number=A,' lift.decon.lifted.vcf)" "2" "the lifted VCF declares GR and PC"

is "$(grep -v '^##' lift.decon.lifted.vcf | grep -v '^gref' | cut -f 1-7,9-)" "$(grep -v '^##' lift.decon.vcf | cut -f 1-7,9-)" "lifting a VCF keeps its other columns, samples included"

is "$(grep '^gref' lift.decon.lifted.vcf | cut -f 1-5,10 | tr '\t' ' ') $(grep '^gref' lift.decon.lifted.vcf | cut -f 8 | tr ';' '\n' | grep -E '^(AT|AN|CH|RC|RS)=' | tr '\n' ' ')" "gref_x_1_alt 50 gref_x_1_alt-50-SNV->3>5>6-1 A G 0|1 AT=>3>4>6,>3>5>6 AN=2 CH=1 RC=x RS=70 " "a SNP at the end of an insertion, where the walk leaves the alt contig, is broken out as a phased record on the contig"

is "$(grep -c '^##contig=<ID=gref_x_1_alt,length=50>' lift.decon.lifted.vcf)" "1" "the lifted VCF declares the alt contigs of its nested records"

is "$(cut -f 1-6 lift.decon.bed | tr '\t\n' ' ;')" "gref_x_1_alt 0 50 >2>6:1 0 +;s#2#c 1070 1120 >2>6:2 0 +;s#2#c 1119 1120 gref_x_1_alt-50-SNV->3>5>6-1 0 +;" "the canonical BED names each allele by its record ID and number, and each nested variant by its ID"

# PanGenie calls: the ALT walk is in the ID, and there is no AT
INS1=$(grep -v '^#' lift.decon.vcf | cut -f 5 | cut -d , -f 1)
INS2=$(grep -v '^#' lift.decon.vcf | cut -f 5 | cut -d , -f 2)
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > lift.pg.vcf
printf 'x\t70\tx-70-INS->2>3>4>6-50\tT\t%s\t.\tPASS\tAC=1\n' "$INS1" >> lift.pg.vcf
printf 'x\t70\tx-70-INS->2>3>5>6-50\tT\t%s\t.\tPASS\t.\n' "$INS2" >> lift.pg.vcf
printf 'x\t70\tx-70-DEL->2>6-1\tTA\tT\t.\tPASS\t.\nx\t70\tnowalk\tT\tTA\t.\tPASS\t.\n' >> lift.pg.vcf
vg lift lift.vcf.gbz lift.pg.vcf > lift.pg.lifted.vcf 2> lift.pg.err
is "$(grep '^x' lift.pg.lifted.vcf | cut -f 3,8 | tr '\t\n' ' ;')" "x-70-INS->2>3>4>6-50 AC=1;GR=gref_x_1_alt:1-50:+:50=;PC=gref_x_1_alt:1-50:+;x-70-INS->2>3>5>6-50 GR=gref_x_1_alt:1-49:+:49=1I;PC=s#2#c:1071-1120:+;x-70-DEL->2>6-1 .;nowalk .;" "PanGenie IDs give the walks of calls without AT, and a deletion has no bases to place"

is "$(grep '^gref' lift.pg.lifted.vcf | awk -F '\t' '{ print $1, $2, $4, $5, NF }')" "gref_x_1_alt 50 A G 8" "calls without samples give nested records without genotypes"

is "$(grep -c '1 ALT alleles have no walk' lift.pg.err)" "1" "alleles without walks are counted"

bgzip -c lift.pg.vcf > lift.pg.vcf.gz
vg lift lift.vcf.gbz lift.pg.vcf.gz 2> /dev/null | cmp -s - lift.pg.lifted.vcf
is $? 0 "a compressed VCF gives the same output"

vg lift lift.vcf.gbz lift.pg.lifted.vcf 2> /dev/null | cmp -s - lift.pg.lifted.vcf
is $? 0 "lifting a lifted VCF replaces its tags"

printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nx\t70\tx-70-INS->2>3>4>6-50\tT\tTAAAA\t.\tPASS\t.\n' > lift.foreign.vcf
vg lift lift.vcf.gbz lift.foreign.vcf 2>&1 > /dev/null | grep -q "is the VCF from this graph"
is $? 0 "an allele that is not in the sequence of its walk gets a warning"

printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nx\t70\tx-70-INS->2>99>6-1\tT\tTA\t.\tPASS\t.\n' > lift.foreign.vcf
vg lift lift.vcf.gbz lift.foreign.vcf 2>&1 > /dev/null | grep -q "node 99 is not in the graph"
is $? 0 "a walk through a node that is not in the graph is an error"

# Nested variants: a 63 bp insertion with three SNP sites, as PanGenie calls of two
# samples. Haplotype s#1 has the insertion that becomes the gref alt contig, s#2 has all
# three SNPs, t#2 has the middle one, and t#1 lacks the insertion.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTGACCAGTAGGCAT\nS\t3\tA\nS\t13\tG\nS\t4\tCCGATTAGCATGACG\n' > lift.snp3.gfa
printf 'S\t5\tC\nS\t15\tT\nS\t6\tGGTACCTAGGATCCA\nS\t7\tG\nS\t17\tA\nS\t8\tTAGCTAGGCTTACGA\nS\t9\tCATGCATGCATGCATGCATG\n' >> lift.snp3.gfa
for edge in 1,2 1,9 2,3 2,13 3,4 13,4 4,5 4,15 5,6 15,6 6,7 6,17 7,8 17,8 8,9; do
    printf 'L\t%s\t+\t%s\t+\t0M\n' ${edge%,*} ${edge#*,} >> lift.snp3.gfa
done
printf 'P\tx\t1+,9+\t*\nW\ts\t1\tc\t1000\t1103\t>1>2>3>4>5>6>7>8>9\nW\ts\t2\tc\t1000\t1103\t>1>2>13>4>15>6>17>8>9\n' >> lift.snp3.gfa
printf 'W\tt\t1\tc\t1000\t1040\t>1>9\nW\tt\t2\tc\t1000\t1103\t>1>2>3>4>15>6>7>8>9\n' >> lift.snp3.gfa
vg paths -x lift.snp3.gfa -Q x --compute-gref > lift.snp3.gref.vg
vg convert -f lift.snp3.gref.vg > lift.snp3.gref.gfa
vg gbwt -G lift.snp3.gref.gfa --gbz-format -g lift.snp3.gbz
# The ALT sequence of a walk, after the padding base T
ins() {
    for node in "$@"; do
        awk -v node=$node '$1 == "S" && $2 == node { printf "%s", $3 }' lift.snp3.gfa
    done
}
call() {
    printf 'x\t20\tx-20-INS-%s-63\tT\tT%s\t.\tPASS\t.\tGT\t%s\t%s\n' "$1" "$(ins $(echo $1 | tr '<>' '  ' | cut -d ' ' -f 3-10))" "$2" "$3"
}
printf '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' > lift.snp3.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\tt\n' >> lift.snp3.vcf
call '>1>2>3>4>5>6>7>8>9' 0/1 0/0 >> lift.snp3.vcf
call '>1>2>13>4>15>6>17>8>9' 0/1 0/0 >> lift.snp3.vcf
call '>1>2>3>4>15>6>7>8>9' 0/0 0/1 >> lift.snp3.vcf
vg lift lift.snp3.gbz lift.snp3.vcf > lift.snp3.lifted.vcf
is "$(grep '^gref' lift.snp3.lifted.vcf | cut -f 1,2,4,5,10,11 | tr '\t\n' ' ;')" "gref_x_1_alt 16 A G 0/1 0/.;gref_x_1_alt 32 C T 0/1 1/.;gref_x_1_alt 48 G A 0/1 0/.;" "SNPs in a called insertion are broken out as records on its alt contig, with . for haplotypes without the contig"

is "$(grep '^gref' lift.snp3.lifted.vcf | cut -f 8 | tr ';' '\n' | grep -E '^(AC|AN)=' | tr '\n' ' ')" "AC=1 AN=3 AC=2 AN=3 AC=1 AN=3 " "nested allele counts are among the haplotypes that carry the contig"

sed 's/x-20-INS->1>2>13>4>15>6>17>8>9-63/x-20-INS-<9<8<17<6<15<4<13<2<1-63/' lift.snp3.vcf > lift.snp3.rev.vcf
vg lift lift.snp3.gbz lift.snp3.rev.vcf 2> /dev/null | grep '^gref' | cmp -s - <(grep '^gref' lift.snp3.lifted.vcf)
is $? 0 "a walk that reads the alt contig backwards gives the same nested records"

awk -F '\t' -v OFS='\t' '$3 == "x-20-INS->1>2>3>4>5>6>7>8>9-63" { $10 = "1/1" } { print }' lift.snp3.vcf > lift.snp3.conflict.vcf
vg lift lift.snp3.gbz lift.snp3.conflict.vcf > lift.snp3.conflict.lifted.vcf 2> lift.snp3.conflict.err
is "$(grep '^gref' lift.snp3.conflict.lifted.vcf | cut -f 10 | sort -u) $(grep -c 'CONFLICT=s;' lift.snp3.conflict.lifted.vcf) $(grep -c '3 nested genotypes are missing' lift.snp3.conflict.err)" "./. 3 1" "a sample with more walks along a contig than haplotypes gets missing genotypes"

# Contig ends: t#2's insertion is the alt contig, s#1 lacks its first base, and s#2 leaves
# it after that base
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTACGTACGTACGT\nS\t21\tT\nS\t22\tCAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGAT\n' > lift.ends.gfa
printf 'S\t9\tCATGCATGCATGCATGCATG\n' >> lift.ends.gfa
for edge in 1,21 1,22 1,9 21,22 21,9 22,9; do
    printf 'L\t%s\t+\t%s\t+\t0M\n' ${edge%,*} ${edge#*,} >> lift.ends.gfa
done
printf 'P\tx\t1+,9+\t*\nW\ts\t1\tc\t1000\t1100\t>1>22>9\nW\ts\t2\tc\t1000\t1041\t>1>21>9\n' >> lift.ends.gfa
printf 'W\tt\t1\tc\t1000\t1040\t>1>9\nW\tt\t2\tc\t1000\t1101\t>1>21>22>9\n' >> lift.ends.gfa
vg paths -x lift.ends.gfa -Q x --compute-gref > lift.ends.gref.vg
vg convert -f lift.ends.gref.vg > lift.ends.gref.gfa
vg gbwt -G lift.ends.gref.gfa --gbz-format -g lift.ends.gbz
vg deconstruct -P x -a lift.ends.gbz > lift.ends.vcf 2> /dev/null
vg lift lift.ends.gbz lift.ends.vcf 2> /dev/null | grep '^gref' | awk -F '\t' '{ print $2, $4, length($4), $5, $10, $11 }' > lift.ends.txt
is "$(cut -d ' ' -f 1,3- lift.ends.txt | tr '\n' ';') $(head -1 lift.ends.txt | cut -d ' ' -f 2)" "1 2 C 1|0 .|0;1 61 T 0|1 .|0; TC" "variants at the start and the end of an alt contig are broken out, padded after at the start"

# A walk on the reference that shares a contig's last node doesn't read the contig: node 4 is a
# SNP allele that s#1 carries and that ends t#2's insertion, the alt contig
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTACGTACGTACGT\nS\t2\tGATTACAGATTACAGATTAC\nS\t3\tA\nS\t4\tG\nS\t9\tCATGCATGCATGCATGCATG\n' > lift.touch.gfa
printf 'S\t5\tTTTCCTCATGCAATTCAAAACCATGTCCGT\nS\t6\tA\nS\t7\tCC\nS\t8\tAATGTAGGCGAAATAGTAAACCATTTTAC\n' >> lift.touch.gfa
for edge in 1,2 2,3 3,9 2,4 4,9 1,5 5,6 5,7 6,8 7,8 8,4; do
    printf 'L\t%s\t+\t%s\t+\t0M\n' ${edge%,*} ${edge#*,} >> lift.touch.gfa
done
printf 'P\tx\t1+,2+,3+,9+\t*\nW\ts\t1\tc\t1000\t1061\t>1>2>4>9\nW\ts\t2\tc\t1000\t1061\t>1>2>3>9\n' >> lift.touch.gfa
printf 'W\tt\t1\tc\t1000\t1101\t>1>5>6>8>4>9\nW\tt\t2\tc\t1000\t1102\t>1>5>7>8>4>9\n' >> lift.touch.gfa
vg paths -x lift.touch.gfa -Q x --compute-gref > lift.touch.gref.vg
vg convert -f lift.touch.gref.vg > lift.touch.gref.gfa
vg gbwt -G lift.touch.gref.gfa --gbz-format -g lift.touch.gbz
vg deconstruct -P x -a lift.touch.gbz > lift.touch.vcf 2> /dev/null
is "$(vg lift lift.touch.gbz lift.touch.vcf 2> /dev/null | grep '^gref' | cut -f 2,4,5,10,11 | tr '\t' ' ' | tr '\n' ';')" "31 CC A ./. 1|0;" "a walk on the reference that shares a contig's last node does not read the contig"

# A contig inside a contig: s#2 and t#1 detour from gref_x_1_alt onto gref_x_2_alt, which
# has a SNP between them
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTACGTACGTACGT\nS\t9\tCATGCATGCATGCATGCATG\n' > lift.nest.gfa
printf 'S\t2\tTTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACG\n' >> lift.nest.gfa
printf 'S\t3\tTGGCTAGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCAT\n' >> lift.nest.gfa
printf 'S\t5\tGGCCCCCCACGATCAGCAGTTCGGCTTGTGAGGTCTTCGCCGGGTGGTCTCCCGCATTTA\n' >> lift.nest.gfa
printf 'S\t4\tACGGGATGTTTAGCGGGGCC\nS\t6\tA\nS\t7\tC\nS\t8\tCAACCAACGCAGTGGTGGCCGGCGTCTTT\n' >> lift.nest.gfa
for edge in 1,2 1,9 2,3 2,4 4,6 4,7 6,8 7,8 8,3 3,5 3,9 5,9; do
    printf 'L\t%s\t+\t%s\t+\t0M\n' ${edge%,*} ${edge#*,} >> lift.nest.gfa
done
printf 'P\tx\t1+,9+\t*\nW\ts\t1\tc\t1000\t1220\t>1>2>3>5>9\nW\ts\t2\tc\t1000\t1210\t>1>2>4>6>8>3>9\n' >> lift.nest.gfa
printf 'W\tt\t1\tc\t1000\t1210\t>1>2>4>7>8>3>9\nW\tt\t2\tc\t1000\t1040\t>1>9\n' >> lift.nest.gfa
vg paths -x lift.nest.gfa -Q x --compute-gref > lift.nest.gref.vg
vg convert -f lift.nest.gref.vg > lift.nest.gref.gfa
vg gbwt -G lift.nest.gref.gfa --gbz-format -g lift.nest.gbz
vg deconstruct -P x -a lift.nest.gbz > lift.nest.vcf 2> /dev/null
vg lift lift.nest.gbz lift.nest.vcf 2> /dev/null | grep '^gref' > lift.nest.lifted.txt
is "$(awk -F '\t' '$1 == "gref_x_2_alt" { print $2, $4, $5, $10, $11 }' lift.nest.lifted.txt) $(grep '^gref_x_2_alt' lift.nest.lifted.txt | grep -o 'CH=[0-9]*')" "21 A C .|0 1|. CH=2" "a SNP on a contig inside a contig is broken out on the inner contig"

is "$(awk -F '\t' '$1 == "gref_x_1_alt" { print $2, length($5) }' lift.nest.lifted.txt | tr '\n' ';')" "60 51;60 51;120 1;" "the detour onto the inner contig is an insertion on the outer one"

# A haplotype that goes through the insertion twice counts once
vg paths -x nesting/nested_snp_in_ins_cycle.gfa -Q x --compute-gref -l 1 > lift.cycle.gref.vg 2> /dev/null
vg convert -f lift.cycle.gref.vg > lift.cycle.gref.gfa
vg gbwt -G lift.cycle.gref.gfa --gbz-format -g lift.cycle.gbz 2> /dev/null
vg deconstruct -P x -a lift.cycle.gbz > lift.cycle.vcf 2> /dev/null
is "$(vg lift lift.cycle.gbz lift.cycle.vcf 2> /dev/null | grep '^gref' | cut -f 1,2,4,5,10,11 | tr '\t' ' ')" "gref_x_1_alt 2 A T 0|1 0|1" "a walk that passes along a contig twice counts once"

rm -f lift.vcf.gfa lift.vcf.gref.vg lift.vcf.gref.gfa lift.vcf.gbz lift.decon.vcf lift.decon.bed lift.decon.lifted.vcf
rm -f lift.pg.vcf lift.pg.vcf.gz lift.pg.lifted.vcf lift.pg.err lift.foreign.vcf
rm -f lift.snp3.gfa lift.snp3.gref.vg lift.snp3.gref.gfa lift.snp3.gbz lift.snp3.vcf lift.snp3.lifted.vcf lift.snp3.rev.vcf
rm -f lift.snp3.conflict.vcf lift.snp3.conflict.lifted.vcf lift.snp3.conflict.err
rm -f lift.ends.gfa lift.ends.gref.vg lift.ends.gref.gfa lift.ends.gbz lift.ends.vcf lift.ends.txt
rm -f lift.touch.gfa lift.touch.gref.vg lift.touch.gref.gfa lift.touch.gbz lift.touch.vcf
rm -f lift.nest.gfa lift.nest.gref.vg lift.nest.gref.gfa lift.nest.gbz lift.nest.vcf lift.nest.lifted.txt
rm -f lift.cycle.gref.vg lift.cycle.gref.gfa lift.cycle.gbz lift.cycle.vcf

rm -f lift.paths.bed lift.paths.exact.bed lift.paths.lift.bed lift.paths.blocks.bed lift.outside.bed
rm -f lift.ins.gbz lift.ins.bed lift.ins.blocks.bed lift.ins.lift.bed lift.ins.tigs.gaf
rm -f lift.names.input.bed lift.names.bed lift.names.gaf lift.names.err lift.tigs.noheader.gaf lift.shared.bed

rm -f lift.gfa lift.gref.vg lift.gref.gfa lift.gbz lift.tigs.gaf lift.tigs.fa lift.bed lift.gaf lift.asm.bed
rm -f lift.place.txt lift.feature.txt lift.support.txt lift.s1.bed lift.t1.bed lift.t4.bed lift.t1.gaf lift.t4.gaf
rm -f lift.unknown.bed lift.long.bed lift.plain.gbz lift.plain.gaf lift.plain.fa lift.plain.bed lift.plain.out.gaf lift.plain.asm.bed

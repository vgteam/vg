#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 55

vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 11 x.vg

is  "$(vg map -s GCTGTGAAGATTAAATTAGGTGAT -x x.xg -g x.gcsa -j - | jq -r '.path.mapping[0].position.offset')" "3" "offset counts unused bases from the start of the node on the forward strand"

is  "$(vg map -s GCTGTGAAGATTAAATTAGGTGAT -x x.xg -g x.gcsa --xdrop-alignment -j - | jq -r '.path.mapping[0].position.offset')" "3" "xdrop alignment obtains the expected result"

is  "$(vg map -s ATCACCTAATTTAATCTTCACAGC -x x.xg -g x.gcsa --xdrop-alignment -j - | jq -r '.path.mapping[0].position.offset')" "5" "xdrop alignment obtains the expected result for the reverse complement"

is  "$(vg map --score-matrix default.mat -s GCTGTGAAGATTAAATTAGGTGAT -x x.xg -g x.gcsa -j - | jq -r '.path.mapping[0].position.offset')" "3" "score-matrix defaults match 1 mismatch -4 should produce same results: with matrixoffset counts unused bases from the start of the node on the forward strand"

is  "$(vg map -s ATCACCTAATTTAATCTTCACAGC -x x.xg -g x.gcsa -j - | jq -r '.path.mapping[0].position.offset')" "5" "offset counts unused bases from the start of the node on the reverse strand"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -j | tr ',' '\n' | grep node_id | grep "72\|73\|76\|77" | wc -l) 4 "global alignment traverses the correct path"

is "$(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -j | jq -r '.score')" 58 "alignment score is as expected"

is "$(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 -x x.xg -g x.gcsa -j | jq -r '.score')" 106 "scoring parameters are respected"

is "$(vg map --score-matrix 2_2.mat -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --gap-open 3 --gap-extend 1 -x x.xg -g x.gcsa -j | jq -r '.score')" 106 "score-matrix file should give same results as --match 2 --mismatch 2: scoring parameters are respected"

is "$(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --full-l-bonus 5 -x x.xg -g x.gcsa -j | jq -r '.score')" 58 "full length bonus always be included"

is "$(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 --full-l-bonus 0 -x x.xg -g x.gcsa -j | jq -r '.score')" 96 "full length bonus can be set to 0"

vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -d x >/dev/null
is $? 0 "vg map takes -d as input without a variant graph"

is $(vg map -s TCAGATTCTCATCCCTCCTCAAGGGCGTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -x x.xg -g x.gcsa -j | jq -r . | grep '"sequence": "G"' | wc -l) 1 "vg map can align across a SNP"

is $(vg map --reads <(vg sim -n 1000 -l 100 -x x.xg) -x x.xg -g x.gcsa  | vg view -a - | jq -r -c '.score == 110 // [.score, .sequence]' | grep true | wc -l) 1000 "alignment works on a small graph"

seq=TCAGATTCTCATCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG
is $(vg map -s $seq -x x.xg -g x.gcsa | vg view -a - | jq -r -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   $(vg map -s $seq -j -x x.xg -g x.gcsa | jq -r -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   "binary alignment format is equivalent to json version"

is $(vg map -b small/x.bam -x x.xg -g x.gcsa -j | jq -r .quality | grep null | wc -l) 0 "alignment from BAM correctly handles qualities"

is $(vg map -s $seq -w 30 -x x.xg -g x.gcsa -j | wc -l) 1 "chunky-banded alignment works"

scores=$(vg map -s GCACCAGGACCCAGAGAGTTGGAATGCCAGGCATTTCCTCTGTTTTCTTTCACCG -x x.xg -g x.gcsa -j -M 2 | jq -r '.score' | tr '\n' ',')
is "${scores}" $(printf ${scores} | tr ',' '\n' | sort -nr | tr '\n' ',')  "multiple alignments are returned in descending score order"

is "$(vg map -s GCACCAGGACCCAGAGAGTTGGAATGCCAGGCATTTCCTCTGTTTTCTTTCACCG -x x.xg -g x.gcsa -j -M 2 | jq -r -c 'select(.is_secondary | not)' | wc -l)" "1" "only a single primary alignment is returned"

vg map -d x -f small/r.fa > /dev/null
is $? 0 "vg can map multi-line FASTA sequences" 

rm -f x.vg x.xg x.gcsa x.gcsa.lcp
rm -rf x.vg.index

vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 >giab.vg
vg index -x giab.xg -g giab.gcsa -k 16 giab.vg

is $(vg map -K -b minigiab/NA12878.chr22.tiny.bam -x giab.xg -g giab.gcsa | vg view -a - | wc -l) $(samtools view minigiab/NA12878.chr22.tiny.bam | wc -l) "mapping of BAM file produces expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa | vg view -a - | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | grep ^@ | wc -l) "mapping from a fastq produces the expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -M 2 -j | jq -r -c 'select(.is_secondary | not)' | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -M 1 -j | wc -l) "allowing secondary alignments with MEM mapping does not change number of primary alignments"

count_prev=$(samtools sort -n minigiab/NA12878.chr22.tiny.bam -T sorting.tmp | samtools bam2fq - 2>/dev/null | vg map -if - -x giab.xg -g giab.gcsa | vg view -a - | jq -r .fragment_prev.name | grep null | wc -l)
count_next=$(samtools sort -n minigiab/NA12878.chr22.tiny.bam -T sorting.tmp | samtools bam2fq - 2>/dev/null | vg map -if - -x giab.xg -g giab.gcsa | vg view -a - | jq -r .fragment_next.name | grep null | wc -l)

is $count_prev $count_next "vg connects paired-end reads in gam output"

rm -f giab.vg giab.xg giab.gcsa giab.gcsa.lcp sorting.tmp*

#vg index -s -k 27 -e 7 graphs/199754000:199755000.vg

#a=$(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 1000 -j | jq -r '.path.mapping[0].position.offset' -c)
#b=$(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 500 -j | jq -r '.path.mapping[0].position.offset' -c)
#is $a $b "banded alignment works correctly even with varied band size"

#c=$(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 1000 -j | jq -r -c '.path.mapping[0].position.offset')
#is $c 29 "unitig mapping produces the correct position"

#is $(for i in $(seq 500 50 2000); do vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B $i -j | jq -r '.path.mapping[0].position.offset' -c; done | sort | uniq | wc -l) 1 "varying the bandwidth does not change the mapping start position"

#rm -rf graphs/199754000:199755000.vg.index

# I was having a problem when updating an edge due to a flipped end node made it
# identical to an already existing edge that hadn't yet been updated. This makes
# sure that that isn't happening.
vg paths -d -v cyclic/orient_must_swap_edges.vg > e.vg
vg index -x e.xg -g e.gcsa -k 10 e.vg
vg map -s "ACACCTCCCTCCCGGACGGGGCGGCTGGCC" -d e >/dev/null
is $? 0 "mapping to graphs that can't be oriented without swapping edges works correctly"

rm -Rf e.xg e.gcsa e.vg

vg index -k10 -g g.idx.gcsa -x g.idx.xg graphs/multimap.vg
is $(vg map -M 2 -s "GCTAAGAGTAGGCCGGGGGTGTAGACCTTTGGGGTTGAATAAATCTATTGTACTAATCGG" -d g.idx -j | jq -r -c 'select(.is_secondary == true)' | wc -l) 1 "reads multi-map to multiple possible locations"

rm -f g.idx.gcsa g.idx.xg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 16 x.vg
is $(vg map -T small/x-s1337-n1000.reads -d x -j | jq -r -c '.path.mapping[0].position.node_id' | wc -l) 1000 "vg map works based on gcsa and xg indexes"

is $(vg map -T <(head -1 small/x-s1337-n1000.reads) -d x -j -t 1 -Q 30 | jq -r .mapping_quality) 30 "the mapping quality may be capped"

vg index -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -k 16 graphs/refonly-lrc_kir.vg

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -u 4 -j  > temp_paired_alignment.json
vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -u 4 -j  > temp_independent_alignment.json
paired_score=$(jq -r ".score" < temp_paired_alignment.json | awk '{ sum+=$1} END {print sum}')
independent_score=$(jq -r ".score" < temp_independent_alignment.json | awk '{ sum+=$1} END {print sum}')
is $(printf "%s\t%s\n" $paired_score $independent_score | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent have lower score than unrestricted alignments"

paired_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_paired_alignment.json| sort | rs -T | awk '{print ($2 - $1)}')
independent_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_independent_alignment.json| sort | rs -T | awk '{print ($2 - $1)}')
is $(printf "%s\t%s\n" $paired_range $independent_range | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent are closer together in node id space than unrestricted alignments"
is $(vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -j -M 4 -UI 832:332.192:50.0602:0:1 | jq -r ".mapping_quality" | grep -v null | wc -l) 2 "only primary alignments have mapping quality scores"
is $(vg map -T small/x-s1337-n1000.reads -x x.xg -g x.gcsa -k 22 -j | jq -r ".mapping_quality" | wc -l) 1000 "unpaired reads produce mapping quality scores"

rm temp_paired_alignment.json temp_independent_alignment.json

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -j -t 1 --surject-to sam > temp_paired_alignment_xg.sam
vg map -x graphs/refonly-lrc_kir.vg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -j -t 1 --surject-to sam > temp_paired_alignment_vg.sam
diff temp_paired_alignment_xg.sam temp_paired_alignment_vg.sam
is "$?" 0 "map returns same output on xg and vg (single threaded)"

rm -f temp_paired_alignment_xg.sam temp_paired_alignment_vg.sam

vg map -f alignment/mismatch_full_qual.fq -x x.xg -g x.gcsa -k 22 -j -A | jq -r -c '.score' > temp_scores_full_qual.txt
vg map -f alignment/mismatch_reduced_qual.fq -x x.xg -g x.gcsa -k 22 -j -A | jq -r -c '.score' > temp_scores_reduced_qual.txt
is $(paste temp_scores_full_qual.txt temp_scores_reduced_qual.txt | column -s $'\t' -t | awk '{if ($1 < $2) count++} END{print count}') 10 "base quality adjusted alignment produces higher scores if mismatches have low quality"
rm temp_scores_full_qual.txt temp_scores_reduced_qual.txt

# Make sure map checks input files to be sure they're there
vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/NONEXISTENT -u 4 -j
is $? 1 "error on vg map -f <nonexistent-file> (unpaired)"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/NONEXISTENT -i -u 4 -j
is $? 1 "error on vg map -f <nonexistent-file> (interleaved)"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -f reads/NONEXISTENT -u 4 -j
is $? 1 "error on vg map -f <nonexistent-file> (paired, LHS)"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/NONEXISTENT -f reads/grch38_lrc_kir_paired.fq -u 4 -j
is $? 1 "error on vg map -f <nonexistent-file> (paired, RHS)"

# Now do the GBWT
vg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -v small/x.vcf.gz --gbwt-name x.gbwt -k 16 x.vg

# This read is all ref which matches no haplotype in x.vcf.gz and visits some unused nodes
is "$(vg map -x x.xg -g x.gcsa --gbwt-name x.gbwt --hap-exp 1 --full-l-bonus 0 -f reads/x.unvisited.fq -j | jq -r '.score')" "36" "mapping a read that touches unused nodes gets the base score"
# This read is all alt which does match a haplotype
is "$(vg map -x x.xg -g x.gcsa --gbwt-name x.gbwt --hap-exp 1 --full-l-bonus 0 -f reads/x.match.fq -j | jq -r '.score')" "35" "mapping a read that matches a haplotype gets a small penalty"
is "$(vg map -x x.xg -g x.gcsa --gbwt-name x.gbwt --hap-exp 0 --full-l-bonus 0 -f reads/x.match.fq -j | jq -r '.score')" "36" "mapping a read that matches a haplotype with exponent 0 gets the base score"
# This read matches no haplotypes but only visits used nodes
is "$(vg map -x x.xg -g x.gcsa --gbwt-name x.gbwt --hap-exp 1 --full-l-bonus 0 -f reads/x.offhap.fq -j | jq -r '.score')" "21" "mapping a read that matches no haplotypes gets a larger penalty"

# Test paired surjected mapping
vg map -d x -iG <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --surject-to SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of paired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '762\n321')" "surjection of paired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "1" "surjection of paired reads to SAM yields properly matched QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '=\n=')" "surjection of paired reads to SAM produces correct pair partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '163\n83')" "surjection of paired reads to SAM produces correct flags"

# And unpaired surjected mapping
vg map -d x -G <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --surject-to SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of unpaired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '0\n0')" "surjection of unpaired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "2" "surjection of unpaired reads to SAM yields distinct QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '*\n*')" "surjection of unpaired reads to SAM produces absent partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '0\n16')" "surjection of unpaired reads to SAM produces correct flags"

rm -f x.vg.idx x.vg.gcsa x.vg.gcsa.lcp x.vg x.xg x.gcsa x.gcsa.lcp x.gbwt graphs/refonly-lrc_kir.vg.xg graphs/refonly-lrc_kir.vg.gcsa graphs/refonly-lrc_kir.vg.gcsa.lcp surjected.sam

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz >tiny.vg
vg index -k 16 -x tiny.xg -g tiny.gcsa tiny.vg
is $(vg map -d tiny -f tiny/tiny.fa -j | jq -r .identity | head -c1) 1 "mapper can read FASTA input"
cat <<EOF >t.fa
>x
CAAATAAGGCTTGGAAATTTTCTGGA
GTTCTATTATATTCCAACTCTCTG
EOF
is $(vg map -d tiny -F t.fa -j | jq -r .identity | head -c1) 1 "mapper can read multiline FASTA input"

rm -f tiny.vg tiny.xg tiny.gcsa tiny.gcsa.lcp t.fa t.fa.fai

vg view -Fv graphs/revdrop.gfa >x.vg
vg index -x x.xg -g x.gcsa x.vg
is $(vg map -s GGTAGGGAACATTAAGGGTATGGAATTGGCAGGACAAGGCACCTGACTGGATTGGGAGAGATAAAGAGGAAAAGCGTCGAGAATGAGCTTGGTGCACTTTGGGCACAGGTGAGTATGCAGAGCGCAACAGGAGGCCTTGGGAACTCATAA -d x -j --xdrop | jq .score) 130 "xdrop works on a reversed read"

rm -f x.vg x.xg x.gcsa x.gcsa.lcp

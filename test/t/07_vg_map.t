#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 33

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 11 x.vg
vg index -sk 11 -d x.idx x.vg

is  "$(vg map -s GCTGTGAAGATTAAATTAGGTGAT -x x.xg -g x.gcsa -J - | jq '.path.mapping[0].position.offset')" "3" "offset counts unused bases from the start of the node on the forward strand"

is  "$(vg map -s ATCACCTAATTTAATCTTCACAGC -x x.xg -g x.gcsa -J - | jq '.path.mapping[0].position.offset')" "5" "offset counts unused bases from the start of the node on the reverse strand"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -J | tr ',' '\n' | grep node_id | grep "72\|73\|76\|77" | wc -l) 4 "global alignment traverses the correct path"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -J | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 58 "alignment score is as expected"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 -x x.xg -g x.gcsa -J | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 106 "scoring parameters are respected"

vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa >/dev/null
is $? 0 "vg map takes -d as input without a variant graph"

is $(vg map -s TCAGATTCTCATCCCTCCTCAAGGGCGTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -x x.xg -g x.gcsa -J | jq . | grep '"sequence": "G"' | wc -l) 1 "vg map can align across a SNP"

is $(vg map -r <(vg sim -s 69 -n 1000 -l 100 x.vg) -x x.xg -g x.gcsa  | vg view -a - | jq -c '.score == 100 // [.score, .sequence]' | grep -v true | wc -l) 0 "alignment works on a small graph"

seq=TCAGATTCTCATCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG
is $(vg map -s $seq -x x.xg -g x.gcsa | vg view -a - | jq -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   $(vg map -s $seq -J -x x.xg -g x.gcsa | jq -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   "binary alignment format is equivalent to json version"

is $(vg map -b small/x.bam -x x.xg -g x.gcsa -J | jq .quality | grep null | wc -l) 0 "alignment from BAM correctly handles qualities"

is $(vg map -s $seq -B 30 -x x.xg -g x.gcsa | vg surject -x x.xg -s - | wc -l) 4 "banded alignment produces a correct alignment"

scores=$(vg map -s GCACCAGGACCCAGAGAGTTGGAATGCCAGGCATTTCCTCTGTTTTCTTTCACCG -x x.xg -g x.gcsa -J -M 2 | jq -r '.score' | tr '\n' ',')
is "${scores}" $(printf ${scores} | tr ',' '\n' | sort -nr | tr '\n' ',')  "multiple alignments are returned in descending score order"

is "$(vg map -s GCACCAGGACCCAGAGAGTTGGAATGCCAGGCATTTCCTCTGTTTTCTTTCACCG -x x.xg -g x.gcsa -J -M 2 | jq -c 'select(.is_secondary | not)' | wc -l)" "1" "only a single primary alignment is returned"

rm -f x.vg x.xg x.gcsa x.gcsa.lcp
rm -rf x.vg.index

vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 >giab.vg
vg index -x giab.xg -g giab.gcsa -k 16 giab.vg

is $(vg map -K -b minigiab/NA12878.chr22.tiny.bam -x giab.xg -g giab.gcsa | vg view -a - | wc -l) $(samtools view minigiab/NA12878.chr22.tiny.bam | wc -l) "mapping of BAM file produces expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa | vg view -a - | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | grep ^@ | wc -l) "mapping from a fastq produces the expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -M 2 -J | jq -c 'select(.is_secondary | not)' | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -M 1 -J | wc -l) "allowing secondary alignments with MEM mapping does not change number of primary alignments"

count_prev=$(samtools sort -no minigiab/NA12878.chr22.tiny.bam x | samtools bam2fq - 2>/dev/null | vg map -if - -x giab.xg -g giab.gcsa | vg view -a - | jq .fragment_prev.name | grep null | wc -l)
count_next=$(samtools sort -no minigiab/NA12878.chr22.tiny.bam x | samtools bam2fq - 2>/dev/null | vg map -if - -x giab.xg -g giab.gcsa | vg view -a - | jq .fragment_next.name | grep null | wc -l)

is $count_prev $count_next "vg connects paired-end reads in gam output"

rm -f giab.vg giab.xg giab.gcsa giab.gcsa.lcp

#vg index -s -k 27 -e 7 graphs/199754000:199755000.vg

#a=$(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 1000 -J | jq '.path.mapping[0].position.offset' -c)
#b=$(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 500 -J | jq '.path.mapping[0].position.offset' -c)
#is $a $b "banded alignment works correctly even with varied band size"

#c=$(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 1000 -J | jq -c '.path.mapping[0].position.offset')
#is $c 29 "unitig mapping produces the correct position"

#is $(for i in $(seq 500 50 2000); do vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B $i -J | jq '.path.mapping[0].position.offset' -c; done | sort | uniq | wc -l) 1 "varying the bandwidth does not change the mapping start position"

#rm -rf graphs/199754000:199755000.vg.index

# I was having a problem when updating an edge due to a flipped end node made it
# identical to an already existing edge that hadn't yet been updated. This makes
# sure that that isn't happening.
vg index -s -k10 -d e.idx cyclic/orient_must_swap_edges.vg
vg map -s "ACACCTCCCTCCCGGACGGGGCGGCTGGCC" -d e.idx >/dev/null
is $? 0 "mapping to graphs that can't be oriented without swapping edges works correctly"

rm -Rf e.idx

vg index -s -k10 -d g.idx graphs/multimap.vg
is $(vg map -M 2 -s "GCTAAGAGTAGGCCGGGGGTGTAGACCTTTGGGGTTGAATAAATCTATTGTACTAATCGG" -d g.idx -J | jq -c 'select(.is_secondary == true)' | wc -l) 1 "reads multi-map to multiple possible locations"

rm -Rf g.idx

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.vg.idx -g x.vg.gcsa -k 16 -X 2 x.vg
vg index -x x.xg -g x.gcsa -k 11 x.vg
vg sim -s 1337 -n 1000 -x x.vg.idx >x.reads
is $(vg map -r x.reads -x x.vg.idx -g x.vg.gcsa -J -t 1 -J | jq -c '.path.mapping[0].position.node_id' | wc -l) 1000 "vg map works based on gcsa and xg indexes"

is $(vg map -r x.reads -x x.vg.idx -g x.vg.gcsa -J -n 5 -t 1 -J | jq -c '.path.mapping[0].position.node_id' | wc -l) 1000 "mem mapping works"

is $(vg map -r x.reads -x x.vg.idx -g x.vg.gcsa -J -t 1 -J -a | jq -c '.path.mapping[0].position.node_id' | wc -l) 1000 "id-cluster mapping works"

vg index -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -k 16 graphs/refonly-lrc_kir.vg

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -W 300 -u 4 -J  > temp_paired_alignment.json
vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -u 4 -J  > temp_independent_alignment.json
paired_score=$(jq -r ".score" < temp_paired_alignment.json | awk '{ sum+=$1} END {print sum}')
independent_score=$(jq -r ".score" < temp_independent_alignment.json | awk '{ sum+=$1} END {print sum}')
is $(printf "%s\t%s\n" $paired_score $independent_score | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent have lower score than unrestricted alignments"

paired_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_paired_alignment.json| sort | rs -T | awk '{print ($2 - $1)}')
independent_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_independent_alignment.json| sort | rs -T | awk '{print ($2 - $1)}')
is $(printf "%s\t%s\n" $paired_range $independent_range | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent are closer together in node id space than unrestricted alignments"
is $(vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -W 300 -u 4 -J -M 10 -W 750 | jq -r ".mapping_quality" | grep -v null | wc -l) 2 "only primary alignments have mapping quality scores"
is $(vg map -r x.reads -x x.xg -g x.gcsa -k 22 -J | jq -r ".mapping_quality" | wc -l) 1000 "unpaired reads produce mapping quality scores"

rm temp_paired_alignment.json temp_independent_alignment.json

is $(vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -W 300 -u 4 -v 1 -J | jq -r ".mapping_quality") $(vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -U -u 4 -v 2 -J | jq -r ".mapping_quality") "mapping quality approximation is equal to exact calculation in a clear cut case"

vg map -f alignment/mismatch_full_qual.fq -x x.xg -g x.gcsa -k 22 -J -1 | jq -c '.score' > temp_scores_full_qual.txt
vg map -f alignment/mismatch_reduced_qual.fq -x x.xg -g x.gcsa -k 22 -J -1 | jq -c '.score' > temp_scores_reduced_qual.txt
is $(paste temp_scores_full_qual.txt temp_scores_reduced_qual.txt | column -s $'\t' -t | awk '{if ($1 < $2) count++} END{print count}') 10 "base quality adjusted alignment produces higher scores if mismatches have low quality"
rm temp_scores_full_qual.txt temp_scores_reduced_qual.txt

is $(vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -W 300 -u 4 -W 750 -J | jq -r 'select(.name == "ERR194147.679985061/1") | .path.mapping[0].position.node_id') 8121 "paired-end reads are pulled to consistent locations at the cost of non-optimal individual alignments"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/NONEXISTENT -W 300 -u 4 -J
is $? 1 "error on vg map -f <nonexistent-file> (unpaired)"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/NONEXISTENT -i -W 300 -u 4 -J
is $? 1 "error on vg map -f <nonexistent-file> (interleaved)"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -f reads/NONEXISTENT -W 300 -u 4 -J
is $? 1 "error on vg map -f <nonexistent-file> (paired, LHS)"

vg map -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/NONEXISTENT -f reads/grch38_lrc_kir_paired.fq -W 300 -u 4 -J
is $? 1 "error on vg map -f <nonexistent-file> (paired, RHS)"

rm -f x.vg.idx x.vg.gcsa x.vg.gcsa.lcp x.vg x.reads x.xg x.gcsa graphs/refonly-lrc_kir.vg.xg graphs/refonly-lrc_kir.vg.gcsa graphs/refonly-lrc_kir.vg.gcsa.lcp

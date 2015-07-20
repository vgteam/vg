#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 14

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s -k 11 x.vg

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg -J | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "global alignment traverses the correct path"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg -J | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 96 "alignment score is as expected"

vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -d x.vg.index >/dev/null
is $? 0 "vg map takes -d as input without a variant graph"

is $(vg map -s TCAGATTCTCATCCCTCCTCAAGGGCGTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG x.vg -J | jq . | grep '"sequence": "G"' | wc -l) 1 "vg map can align across a SNP"

is $(vg map -r <(vg sim -s 69 -n 1000 -l 100 x.vg) x.vg | vg view -a - | jq -c '.score == 200 // [.score, .sequence]' | grep -v true | wc -l) 0 "alignment works on a small graph"

seq=TCAGATTCTCATCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG
is $(vg map -s $seq x.vg | vg view -a - | jq -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   $(vg map -s $seq -J x.vg | jq -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   "binary alignment format is equivalent to json version"

is $(vg map -b small/x.bam x.vg -J | jq .quality | grep null | wc -l) 0 "alignment from BAM correctly handles qualities"

is $(vg map -s $seq -B 30 x.vg | vg surject -d x.vg.index -s - | wc -l) 4 "banded alignment produces a correct alignment"

rm x.vg
rm -rf x.vg.index

vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz >giab.vg
vg index -s -k 27 -e 7 giab.vg                                                   
is $(vg map -b minigiab/NA12878.chr22.tiny.bam giab.vg | vg view -a - | wc -l) $(samtools view minigiab/NA12878.chr22.tiny.bam | wc -l) "mapping of BAM file produces expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - giab.vg | vg view -a - | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | grep ^@ | wc -l) "mapping from a fastq produces the expected number of alignments"

count_prev=$(samtools sort -n minigiab/NA12878.chr22.tiny.bam -o x | samtools bam2fq - 2>/dev/null | vg map -if - giab.vg | vg view -a - | jq .fragment_prev.name | grep null | wc -l)
count_next=$(samtools sort -n minigiab/NA12878.chr22.tiny.bam -o x | samtools bam2fq - 2>/dev/null | vg map -if - giab.vg | vg view -a - | jq .fragment_next.name | grep null | wc -l)

is $count_prev $count_next "vg connects paired-end reads in gam output"

rm giab.vg
rm -rf giab.vg.index

vg index -s -k 27 -e 7 graphs/199754000:199755000.vg

is $(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 1000 -J | jq '.path.mapping[0].position' -c -S) $(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 500 -J | jq '.path.mapping[0].position' -c -S) "banded alignment works correctly even with varied band size"

is $(vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B 1000 -J | grep '"offset": 29' | wc -l) 1 \
   "unitig mapping position is as expected"

is $(for i in $(seq 500 50 2000); do vg map -f graphs/2086553952_1469228759.mag -d graphs/199754000:199755000.vg.index -B $i -J | jq '.path.mapping[0].position.offset' -c; done | sort | uniq | wc -l) 1 "varying the bandwidth does not change the mapping start position"

rm -rf graphs/199754000:199755000.vg.index

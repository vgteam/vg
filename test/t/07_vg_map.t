#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 19

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -k 11 x.vg
vg index -sk 11 -d x.idx x.vg

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -k 22 -J | tr ',' '\n' | grep node_id | grep "72\|74\|75\|77" | wc -l) 4 "global alignment traverses the correct path"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -k 22 -J | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 48 "alignment score is as expected"

is $(vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG --match 2 --mismatch 2 --gap-open 3 --gap-extend 1 -x x.xg -g x.gcsa -k 22 -J | tr ',' '\n' | grep score | sed "s/}//g" | awk '{ print $2 }') 96 "scoring parameters are respected"

vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa -k 22  >/dev/null
is $? 0 "vg map takes -d as input without a variant graph"

is $(vg map -s TCAGATTCTCATCCCTCCTCAAGGGCGTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -x x.xg -g x.gcsa -k 22 -J | jq . | grep '"sequence": "G"' | wc -l) 1 "vg map can align across a SNP"

is $(vg map -r <(vg sim -s 69 -n 1000 -l 100 x.vg) -x x.xg -g x.gcsa -k 22  | vg view -a - | jq -c '.score == 100 // [.score, .sequence]' | grep -v true | wc -l) 0 "alignment works on a small graph"

seq=TCAGATTCTCATCCCTCCTCAAGGGCTTCTAACTACTCCACATCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG
is $(vg map -s $seq -x x.xg -g x.gcsa -k 22 | vg view -a - | jq -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   $(vg map -s $seq -J -x x.xg -g x.gcsa -k 22 | jq -c '[.score, .sequence, .path.node_id]' | md5sum | awk '{print $1}') \
   "binary alignment format is equivalent to json version"

is $(vg map -b small/x.bam -x x.xg -g x.gcsa -k 22 -J | jq .quality | grep null | wc -l) 0 "alignment from BAM correctly handles qualities"

is $(vg map -s $seq -B 30 -x x.xg -g x.gcsa -k 22 | vg surject -d x.idx -s - | wc -l) 4 "banded alignment produces a correct alignment"

scores=$(vg map -s GCACCAGGACCCAGAGAGTTGGAATGCCAGGCATTTCCTCTGTTTTCTTTCACCG -x x.xg -g x.gcsa -k 22 -J -M 2 | jq -r '.score' | tr '\n' ',')
is "${scores}" $(printf ${scores} | tr ',' '\n' | sort -nr | tr '\n' ',')  "multiple alignments are returned in descending score order"

rm -f x.vg x.xg x.gcsa x.gcsa.lcp
rm -rf x.vg.index

vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 >giab.vg
vg index -x giab.xg -g giab.gcsa -k 11 giab.vg

is $(vg map -K -b minigiab/NA12878.chr22.tiny.bam -k 27 -x giab.xg -g giab.gcsa | vg view -a - | wc -l) $(samtools view minigiab/NA12878.chr22.tiny.bam | wc -l) "mapping of BAM file produces expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -k 27 | vg view -a - | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | grep ^@ | wc -l) "mapping from a fastq produces the expected number of alignments"

is $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -M 2 -J | jq -c 'select(.is_secondary | not)' | wc -l) $(samtools bam2fq minigiab/NA12878.chr22.tiny.bam 2>/dev/null | vg map -f - -x giab.xg -g giab.gcsa -M 1 -J | wc -l) "allowing secondary alignments with MEM mapping does not change number of primary alignments"

count_prev=$(samtools sort -no minigiab/NA12878.chr22.tiny.bam x | samtools bam2fq - 2>/dev/null | vg map -if - -x giab.xg -g giab.gcsa -k 27 | vg view -a - | jq .fragment_prev.name | grep null | wc -l)
count_next=$(samtools sort -no minigiab/NA12878.chr22.tiny.bam x | samtools bam2fq - 2>/dev/null | vg map -if - -x giab.xg -g giab.gcsa -k 27 | vg view -a - | jq .fragment_next.name | grep null | wc -l)

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
vg sim -s 1337 -n 1000 x.vg >x.reads
is $(vg map -r x.reads -x x.vg.idx -g x.vg.gcsa -J -k 16 -t 1 -J | jq -c '.path.mapping[0].position.node_id' | wc -l) 1000 "vg map works based on gcsa and xg indexes"

is $(vg map -r x.reads -x x.vg.idx -g x.vg.gcsa -J -n 5 -t 1 -J | jq -c '.path.mapping[0].position.node_id' | wc -l) 1000 "mem mapping works"

is $(vg map -r x.reads -V x.vg -k 16 -t 1 -J | jq -c '.path.mapping[0].position.node_id' | wc -l) 1000 "vg map can build its own in-memory indexes"

rm -f x.vg.idx x.vg.gcsa x.vg.gcsa.lcp x.vg x.reads

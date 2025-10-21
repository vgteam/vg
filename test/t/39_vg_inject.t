#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 33

vg construct -r small/x.fa > j.vg
vg index -x j.xg j.vg
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index -k 11 -g x.gcsa -x x.xg x.vg

#samtools view -f 0 small/x.bam | awk '$6 ~ /100M/' > p.bam
#samtools view -f 0 small/x.bam | awk '$6 !~ /100M/' > pm.bam
#samtools view -O bam -f 16 small/x.bam > small/i.bam
#samtools view -f 16 small/x.bam | awk '$6 ~ /100M/' > i.bam
#samtools view -f 16 small/x.bam | awk '$6 !~ /100M/' > im.bam

is $(vg inject -x x.xg small/x.bam | vg view -a - | wc -l) \
    1000 "reads are generated"

is $(vg inject -x x.xg small/x.bam | vg surject -x x.xg -t 1 - | vg view -a - | wc -l) \
    1000 "vg inject works for all reads included in the original bam"

is "$(samtools view small/x.bam | cut -f 1 | sort)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 1 | sort)" \
    "vg inject retains read names"

is "$(samtools view small/x.bam | cut -f 4)" "$(vg inject -x j.xg small/x.bam | vg surject -p x -x j.xg -t 1 -s - | samtools view -S - | cut -f 4 | sort -n)" \
    "vg inject works perfectly for the position of alignment in the simple graph"

is "$(samtools view small/x.bam | cut -f 4)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 4 | sort -n)" \
    "vg inject works perfectly for the position of alignment"

is "$(samtools view small/x.bam | cut -f 5 | grep 60 | wc -l)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | cut -f 5 | grep 60 | wc -l)" \
    "vg inject works perfectly for the mapping quality of alignment"

is "$(samtools view small/x.bam | sort -k 1 | cut -f 6)" "$(vg inject -x j.xg small/x.bam | vg surject -p x -x j.xg -t 1 -s - | samtools view -S - | sort -k 1 | cut -f 6)" \
    "vg inject works perfectly for the cigar of alignment in the simple graph"

is "$(samtools view small/x.bam | sort -k 1 | cut -f 6)" "$(vg inject -x x.xg small/x.bam | vg surject -p x -x x.xg -t 1 -s - | samtools view -S - | sort -k 1 | cut -f 6)" \
    "vg inject works perfectly for the cigar of alignment"

is "$(vg inject -x x.xg small/i.bam | vg view -a - | wc -l)" 470 "vg inject preserves all reads"

is "$(vg inject -x x.xg small/i.bam | vg view -a - | jq .path.mapping[0].position.is_reverse | grep -v true | wc -l)" \
    0 "vg inject works perfectly for the reads flagged as is_reverse"

cat <(samtools view -H small/x.bam) <(printf "name\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGT\tHHHHHHHHHHHH\n") > unmapped.sam
is "$(vg inject -x x.xg unmapped.sam | vg view -aj - | grep "path" | wc -l)" 0 "vg inject does not make an alignment for an umapped read"
is "$(echo $?)" 0 "vg inject does not crash on an unmapped read"

is $(vg inject -x x.xg small/x.bam -o GAF | wc -l) \
    1000 "vg inject supports GAF output"

samtools cat small/x.bam small/i.bam > all.bam
is "$(vg inject -x x.xg all.bam | vg validate -A -c -a - x.xg 2>&1 | grep invalid | wc -l)" 0 "vg inject produces valid alignments of the read and reference"

vg inject -x x.xg small/atstart.sam >/dev/null 2>log.txt
is "${?}" 0 "vg inject can inject a read abutting the start of the contig"
is "$(cat log.txt | wc -l)" "0" "vg inject does not warn about a read abutting the start of the contig"

# TODO: We can't see a 0-position-in-the-SAM past-start mapped read properly,
# because htslib prints a warning and fixes it to unmapped for us.
vg inject -x x.xg small/paststart.sam >/dev/null 2>log.txt
is "$(grep '\[W' log.txt | wc -l)" "1" "vg inject warns about a read mapping extending past the start of the contig"

vg inject -x x.xg small/atend.sam >/dev/null 2>log.txt
is "${?}" 0 "vg inject can inject a read abutting the end of the contig"
is "$(cat log.txt | wc -l)" "0" "vg inject does not warn about a read abutting the end of the contig"

vg inject -x x.xg small/pastend.sam >/dev/null 2>log.txt
is "${?}" 1 "vg inject aborts when given a read extending past the end of the contig"
is "$(grep error log.txt | grep 1001 | grep 1002 | wc -l)" "1" "vg inject reports a useful error message about a read extending past the end of the contig"

vg inject -x x.xg small/unmapped.sam | vg stats -a - >stats.txt
is "$(grep "Total alignments: 4" stats.txt | wc -l)" "1" "Injecting reads flagged as unmapped produces alignment records"
is "$(grep "Total aligned: 0" stats.txt | wc -l)" "1" "Injecting reads flagged as unmapped produces unaligned alignment records"

vg inject -x x.xg small/bad_contig.sam >/dev/null 2>log.txt
is "${?}" 1 "vg inject aborts when given a read on a contig not in the graph"
is "$(grep error log.txt | grep 'not present in graph' | wc -l)" "1" \
"vg inject reports a useful error message about a read on a contig not in the graph"

vg inject -t 2 -x x.xg small/bad_contig.sam >/dev/null 2>log.txt
is "${?}" 1 "vg inject aborts when given a read on a contig not in the graph (parallel)"
is "$(grep error log.txt | grep 'not present in graph' | wc -l)" "1" \
"vg inject reports a useful error message about a read on a contig not in the graph (parallel)"

vg inject -a -x x.xg small/bad_contig.sam > bad_contig.gam
is "${?}" 0 "vg inject -a allows a read on a contig not in the graph"
vg stats -a bad_contig.gam >stats.txt
is "$(grep "Total alignments: 2" stats.txt | wc -l)" "1" "Injecting reads on missing contigs still produces alignment records"
is "$(grep "Total aligned: 1" stats.txt | wc -l)" "1" "Injecting reads on missing contigs produces unmapped alignment records"

vg inject -t 2 -a -x x.xg small/bad_contig.sam > bad_contig.gam
is "${?}" 0 "vg inject -a allows a read on a contig not in the graph (parallel)"
vg stats -a bad_contig.gam >stats.txt
is "$(grep "Total alignments: 2" stats.txt | wc -l)" "1" "Injecting reads on missing contigs still produces alignment records (parallel)"
is "$(grep "Total aligned: 1" stats.txt | wc -l)" "1" "Injecting reads on missing contigs produces unmapped alignment records (parallel)"

rm j.vg j.xg x.vg x.gcsa x.gcsa.lcp x.xg unmapped.sam all.bam log.txt stats.txt bad_contig.gam

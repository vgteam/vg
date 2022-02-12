#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 36

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz -a >x2.vg
vg index -x x.xg x.vg
vg index -G x.gbwt -v small/x.vcf.gz x2.vg

is $(vg sim -l 100 -n 100 -x x.xg | wc -l) 100 \
    "vg sim creates the correct number of reads"

is $(vg sim -l 100 -n 1 -e 0.0 -i 0.0 -J -x x.xg | jq .score) 110 "end bonuses are included"

is $(vg sim -l 100 -n 100 -a -x x.xg | vg view -a - | wc -l) 100 \
   "alignments may be generated rather than read sequences"

is $(vg sim -l 100 -n 100 -J -x x.xg | grep -o is_reverse | uniq | wc -l) 1 "alignments are produced on both strands"

is $(vg sim -l 100 -n 100 -e 0.1 -i 0.1 -J -x x.xg | jq .sequence | wc -c) 10300 "high simulated error rates do not change the number of bases generated"

is $(vg sim -l 100 -n 100 -x x.xg -aJ | jq 'select(.path.mapping[0].is_reverse)' | wc -l) 0 \
   "vg sim creates forward-strand reads when asked"

vg view -j x.vg | jq -r '.path[].mapping[].position.node_id' | sort > path.txt
is $(vg sim -l 100 -n 100 -x x.xg -P "x" -aJ | jq -r '.path.mapping[].position.node_id' | sort | uniq | comm -23 - path.txt | wc -l) "0" "vg sim can simulate from just a path"

is $(vg sim -F reads/grch38_lrc_kir_paired.fq -n 100 -x x.xg -P "x" -aJ | jq -r '.path.mapping[].position.node_id' | sort | uniq | comm -23 - path.txt | wc -l) "0" "vg sim can simulate from just a path in FASTQ mode"

vg sim -n 100 -l 100 -A -x x.xg -r > reads.txt 2> log.txt
is $? 0 "reads can be simulated from all paths"
is $(grep -c "Selecting all 1 paths" log.txt) 1 "simulation was successful"

vg sim -n 100 -l 100 -m 1 -g x.gbwt -x x.xg -r > reads.txt 2> log.txt
is $? 0 "reads can be simulated from GBWT samples"
is $(grep -c "Inserted 2 paths" log.txt) 1 "simulation was successful"

rm -f path.txt reads.txt log.txt

is $(vg sim -n 10 -i 0.005 -l 10 -p 50 -v 50 -x x.xg -J | wc -l) 20 "pairs simulated even when fragments overlap"

cat tiny/tiny.fa | sed s/GCTTGGA/GCNTGGA/ >n.fa
vg construct -r n.fa >n.vg
vg index -x n.xg n.vg
is $(vg sim -n 1000 -l 20 -x n.xg | grep N | wc -l) 0 "sim does not emit reads with Ns"

is $(vg sim -N -n 1000 -l 20 -x n.xg | grep -o N | uniq | wc -l) 1 "sim can emit reads with Ns when asked to"

is $(vg sim -n 1000 -l 2 -p 5 -e 0.1 -x n.xg | grep N | wc -l) 0 "sim doesn't emit Ns even with pair and errors"

rm -f x.vg x2.vg x.xg x.gbwt n.vg n.fa n.xg

vg construct -r small/xy.fa -v small/x.vcf.gz -a >xy.vg
vg index -x xy.xg xy.vg
vg index -G xy.gbwt -v small/x.vcf.gz xy.vg

vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 --any-path >/dev/null
is $? "0" "Sample simulation works along with --any-path"

is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:0 | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "not having any regexes skips all unvisited contigs"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex candyfloss:12345 | jq -c 'select(.refpos[].name == "x")' | wc -l | awk '{print ($0>450&&$0<550)?1:0}')" "1" "assigning a ploidy of 2 by default gives reads in a ~1-1 ratio vs. diploid covered contigs"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:1 | jq -c 'select(.refpos[].name == "x")' | wc -l | awk '{print ($0>620&&$0<720)?1:0}')" "1" "assigning a ploidy of 1 by regex gives reads in a ~1-2 ratio vs. diploid covered contigs"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:0 | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "assigning a ploidy of 0 by regex excludes reads from that contig"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex a:5,b:10,y:0,c:6 | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "multiple regexes are interpreted"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex "[^x]:0" | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "regexes are interpreted as regexes"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:1E20 | jq -c 'select(.refpos[].name == "y")' | wc -l)" "1000" "assigning a very high ploidy by regex starves out the other paths"

vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 --any-path -F minigiab/NA12878.chr22.tiny.fq.gz >/dev/null
is $? "0" "Sample simulation works along with --any-path with training FASTQ"

is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:0 -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "not having any regexes skips all unvisited contigs with training FASTQ"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex candyfloss:12345 -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "x")' | wc -l | awk '{print ($0>450&&$0<550)?1:0}')" "1" "assigning a ploidy of 2 by default gives reads in a ~1-1 ratio vs. diploid covered contigs with training FASTQ"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:1 -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "x")' | wc -l | awk '{print ($0>620&&$0<720)?1:0}')" "1" "assigning a ploidy of 1 by regex gives reads in a ~1-2 ratio vs. diploid covered contigs with training FASTQ"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:0 -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "assigning a ploidy of 0 by regex excludes reads from that contig with training FASTQ"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex a:5,b:10,y:0,c:6 -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "multiple regexes are interpreted with training FASTQ"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex "[^x]:0" -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "x")' | wc -l)" "1000" "regexes are interpreted as regexes with training FASTQ"
is "$(vg sim -s 12345 -n 1000 -l 2 -e 0.1 -x xy.xg -g xy.gbwt --sample-name 1 -aJ --ploidy-regex y:1E20 -F minigiab/NA12878.chr22.tiny.fq.gz | jq -c 'select(.refpos[].name == "y")' | wc -l)" "1000" "assigning a very high ploidy by regex starves out the other paths with training FASTQ"

rm -f xy.vg xy.xg xy.gbwt

vg construct -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg x.vg
is "$(vg sim -s 12345 -n 1000 -l 64 -e 0.1 -x x.xg -aJ | jq -c 'select(.refpos | length == 1)' | wc -l)" "1000" "All reads are annotated with single reference positions by default"
is "$(vg sim -s 12345 -n 1000 -l 64 -e 0.1 -x x.xg -aJ --multi-position | jq -c 'select(.refpos | length > 1)' | wc -l)" "1000" "All reads can be annotated with multiple reference positions when requested"

rm -f x.vg x.xg

vg view -Fv graphs/cactus-BRCA2.gfa >cactus-BRCA2.vg
vg index -x cactus-BRCA2.xg cactus-BRCA2.vg
is $(vg sim -x cactus-BRCA2.xg -n 100 -l 150 -p 1000 -v 100 -e 0.01 -i 0.005 -F minigiab/NA12878.chr22.tiny.fq.gz | wc -l) 100 "ngs trained simulator works"
is $(vg sim -x cactus-BRCA2.xg -n 100 -l 150 -p 1000 -v 100 -e 0.01 -i 0.005 -a -F minigiab/NA12878.chr22.tiny.fq.gz | vg view -a - | wc -l) 200 "ngs trained simulator generates gam"
rm -f cactus-BRCA2.xg cactus-BRCA2.vg

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 23

vg construct -r small/x.fa >j.vg
vg index -x j.xg j.vg
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -k 11 -g x.gcsa -x x.xg x.vg

# Simulate some reads from just j
vg map -G <(vg sim -a -s 1337 -n 100 -x j.xg) -g x.gcsa -x x.xg > j.gam
# And some from all of x
vg map -G <(vg sim -a -s 1337 -n 100 -x x.xg) -g x.gcsa -x x.xg > x.gam

is $(vg view -aj j.gam | wc -l) \
    100 "reads are generated"

is $(vg surject -p x -x x.xg -t 1 j.gam | vg view -a - | jq .score | grep 100 | wc -l) \
    100 "vg surject works perfectly for perfect reads derived from the reference"
    
is $(vg surject -p x -x x.xg -t 1 -s j.gam | grep -v "@" | cut -f3 | grep x | wc -l) \
    100 "vg surject actually places reads on the correct path"

is $(vg surject -x x.xg -t 1 -s j.gam | grep -v "@" | cut -f3 | grep x | wc -l) \
    100 "vg surject doesn't need to be told which path to use"

is $(vg surject -p x -x x.xg -t 1 x.gam | vg view -a - | wc -l) \
    100 "vg surject works for every read simulated from a dense graph"

is $(vg surject -p x -x x.xg -s x.gam | grep -v ^@ | wc -l) \
    100 "vg surject produces valid SAM output"

is $(vg map -G <(vg sim -a -s 1337 -n 100 -x x.xg) -g x.gcsa -x x.xg --surject-to sam | grep -v ^@ | wc -l) \
    100 "vg map may surject reads to produce valid SAM output"

is $(vg map -G <(vg sim -a -s 1337 -n 100 -x x.xg) -g x.gcsa -x x.xg --surject-to bam | samtools view - | grep -v ^@ | wc -l) \
    100 "vg map may surject reads to produce valid BAM output"

is $(vg view -aj j.gam | jq '.name = "Alignment"' | vg view -JGa - | vg surject -p x -x x.xg - | vg view -aj - | jq -c 'select(.name)' | wc -l) \
   100 "vg surject retains read names"

is $(vg map -s GTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCTCCTTGAGTCTCCTATGTCCAGGAATGAACCGCTGCT -d x | vg surject -x x.xg -s - | grep 29S | wc -l) 1 "we respect the original mapping's softclips"

# These sequences have edits in them, so we can test CIGAR reversal as well
SEQ="ACCGTCATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAGCTCTCCTTTGTGTCCAGAGTG"
SEQ_RC="CACTCTGGACACAAAGGAGAGCTGTCCTTGAGCCTCGACTCTGGGTGAGCCCTCTGGGTCTTAGATTTGAGATGCAATTTTCAAACTTGAAGATGACGGT"

is "$(vg map -s $SEQ -g x.gcsa -x x.xg | vg surject -p x -x x.xg - -s | cut -f1,3,4,5,6,7,8,9,10)" "$(vg map -s $SEQ_RC -g x.gcsa -x x.xg | vg surject -p x -x x.xg - -s | cut -f1,3,4,5,6,7,8,9,10)" "forward and reverse orientations of a read produce the same surjected SAM, ignoring flags"

is $(vg map -G <(vg sim -a -s 1337 -n 100 -x x.xg) -g x.gcsa -x x.xg | vg surject -p x -x x.xg -b - | samtools view - | wc -l) \
    100 "vg surject produces valid BAM output"

#is $(vg map -G <(vg sim -a -s 1337 -n 100 x.vg) x.vg | vg surject -p x -g x.gcsa -x x.xg -c - | samtools view - | wc -l) \
#    100 "vg surject produces valid CRAM output"

echo '{"sequence": "GATTACA", "path": {"mapping": [{"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}]}, "mapping_quality": 99}' | vg view -JGa - > read.gam
is "$(vg surject -p x -x x.xg read.gam | vg view -aj - | jq '.mapping_quality')" "99" "mapping quality is preserved through surjection"

echo '{"name": "read/2", "sequence": "GATTACA", "path": {"mapping": [{"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}]}, "fragment_prev": {"name": "read/1"}}' | vg view -JGa - > read.gam
is "$(vg surject -p x -x x.xg -i read.gam | vg view -aj - | jq -r '.fragment_prev.name')" "read/1" "read pairing is preserved through GAM->GAM surjection"

vg map -d x -iG <(vg sim -a -s 13241 -n 1 -p 500 -v 300 -x x.xg | vg view -a - | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) | vg surject -x x.xg -p x -s -i - >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort | cut -f 4)" "$(printf '321\n762')" "surjection of paired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort | cut -f 8)" "$(printf '762\n321')" "surjection of paired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "1" "surjection of paired reads to SAM yields properly matched QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '=\n=')" "surjection of paired reads to SAM produces correct pair partner contigs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 2 | sort -n)" "$(printf '83\n131')" "surjection of paired reads to SAM produces correct flags"

rm -rf j.vg x.vg j.gam x.gam x.idx j.xg x.xg x.gcsa read.gam reads.gam surjected.sam

vg mod -c graphs/fail.vg >f.vg
vg index -k 11 -g f.gcsa -x f.xg f.vg

read=TTCCTGTGTTTATTAGCCATGCCTAGAGTGGGATGCGCCATTGGTCATCTTCTGGCCCCTGTTGTCGGCATGTAACTTAATACCACAACCAGGCATAGGTGAAAGATTGGAGGAAAGATGAGTGACAGCATCAACTTCTCTCACAACCTAG
revcompread=CTAGGTTGTGAGAGAAGTTGATGCTGTCACTCATCTTTCCTCCAATCTTTCACCTATGCCTGGTTGTGGTATTAAGTTACATGCCGACAACAGGGGCCAGAAGATGACCAATGGCGCATCCCACTCTAGGCATGGCTAATAAACACAGGAA
is $(vg map -s $read -g f.gcsa -x f.xg | vg surject -p 6646393ec651ec49 -x f.xg -s - | grep $revcompread | wc -l) 1 "surjection works for a longer (151bp) read"

rm -rf f.xg f.gcsa f.vg

vg mod -c graphs/fail2.vg >f.vg
vg index -k 11 -g f.gcsa -x f.xg f.vg

read=TATTTACGGCGGGGGCCCACCTTTGACCCTTTTTTTTTTTCAAGCAGAAGACGGCATACGAGATCACTTCGAGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCTACCCTAACCCTAACCCCAACCCCTAACCCTAACCCCA
is $(vg map -s $read -g f.gcsa -x f.xg | vg surject -p ad93c27f548fc1ae -x f.xg -s - | grep $read | wc -l) 1 "surjection works for another difficult read"

rm -rf f.xg f.gcsa f.vg

vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz >minigiab.vg
vg index -k 11 -g m.gcsa -x m.xg minigiab.vg
is $(vg map -b minigiab/NA12878.chr22.tiny.bam -x m.xg -g m.gcsa | vg surject -p q -x m.xg -s - | grep chr22.bin8.cram:166:6027 | grep BBBBBFBFI | wc -l) 1 "mapping reproduces qualities from BAM input"
is $(vg map -f minigiab/NA12878.chr22.tiny.fq.gz -x m.xg -g m.gcsa | vg surject -p q -x m.xg -s - | grep chr22.bin8.cram:166:6027 | grep BBBBBFBFI | wc -l) 1 "mapping reproduces qualities from fastq input"

rm -rf minigiab.vg* m.xg m.gcsa




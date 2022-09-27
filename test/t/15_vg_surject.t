#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 44

vg construct -r small/x.fa >j.vg
vg index -x j.xg j.vg
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -k 11 -g x.gcsa -x x.xg x.vg

# We have already simulated some reads from just j
vg map -G small/x-allref-nohptrouble.gam -g x.gcsa -x x.xg > j.gam
# Simulate some from all of x
vg map -G <(vg sim -a -n 100 -x x.xg) -g x.gcsa -x x.xg > x.gam

is $(vg view -aj j.gam | wc -l) \
    100 "reads are generated"

# Surjection uses path anchored surject which keeps aligned stuff aligned even if there's a better alignment that shifts it.
# This means arbitrarily chosen homopolymer indel alignment that arbitrarily chose wrong won't be fixed.
# We generate GAMs that don't have that problem.

is $(vg surject -p x -x x.xg -t 1 j.gam | vg view -a - | jq .score | grep 110 | wc -l) \
   100 "vg surject works perfectly for perfect reads without misaligned homopolymer indels derived from the reference"

is $(vg convert x.xg -G j.gam | vg surject -p x -x x.xg -t 1 -G - | vg view -a - | jq .score | grep 110 | wc -l) \
    100 "vg surject works perfectly for perfect reads without misaligned homopolymer indels derived from the reference"
    
is $(vg surject -p x -x x.xg -t 1 -s j.gam | grep -v "@" | cut -f3 | grep x | wc -l) \
    100 "vg surject actually places reads on the correct path"

is $(vg surject -x x.xg -t 1 -s j.gam | grep -v "@" | cut -f3 | grep x | wc -l) \
    100 "vg surject doesn't need to be told which path to use"
    
vg paths -X -x x.vg | vg view -aj - | jq '.name = "sample#0#x#0"' | vg view -JGa - > paths.gam
vg paths -X -x x.vg | vg view -aj - | jq '.name = "ref#0#x[55]"' | vg view -JGa - >> paths.gam
vg augment x.vg -i paths.gam > x.aug.vg
vg index -x x.aug.xg x.aug.vg

is $(vg surject -x x.aug.xg -t 1 -s j.gam | grep -v "@" | cut -f3 | grep "ref#0#x" | wc -l) \
    100 "vg surject picks a reference-sense path if it is present"

rm x.aug.vg x.aug.xg paths.gam

is $(vg surject -p x -x x.xg -t 1 x.gam | vg view -a - | wc -l) \
    100 "vg surject works for every read simulated from a dense graph"

is $(vg surject -S -p x -x x.xg -t 1 x.gam | vg view -a - | wc -l) \
    100 "vg surject spliced algorithm works for every read simulated from a dense graph"

is $(vg surject -p x -x x.xg -s x.gam | grep -v ^@ | wc -l) \
    100 "vg surject produces valid SAM output"

is $(vg map -G <(vg sim -a -n 100 -x x.xg) -g x.gcsa -x x.xg --surject-to sam | grep -v ^@ | wc -l) \
    100 "vg map may surject reads to produce valid SAM output"

is $(vg map -G <(vg sim -a -n 100 -x x.xg) -g x.gcsa -x x.xg --surject-to bam | samtools view - | grep -v ^@ | wc -l) \
    100 "vg map may surject reads to produce valid BAM output"

is $(vg view -aj j.gam | jq '.name = "Alignment"' | vg view -JGa - | vg surject -p x -x x.xg - | vg view -aj - | jq -c 'select(.name)' | wc -l) \
   100 "vg surject retains read names"
   
is $(vg surject -p x -x x.xg j.gam --sample "NA12345" --read-group "RG1" | vg view -aj - | jq -c 'select(.sample_name == "NA12345" and .read_group == "RG1")' | wc -l) \
   100 "vg surject can set sample and read group"

is $(vg map -s GTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCTCCTTGAGTCTCCTATGTCCAGGAATGAACCGCTGCT -d x | vg surject -x x.xg -s - | grep 29S | wc -l) 1 "we respect the original mapping's softclips"

# These sequences have edits in them, so we can test CIGAR reversal as well
SEQ="ACCGTCATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAGCTCTCCTTTGTGTCCAGAGTG"
SEQ_RC="CACTCTGGACACAAAGGAGAGCTGTCCTTGAGCCTCGACTCTGGGTGAGCCCTCTGGGTCTTAGATTTGAGATGCAATTTTCAAACTTGAAGATGACGGT"
QUAL="CCCFFFFFHHHHHJJJJJHFDDDD&((((+>(26:&)()(+((+3((8A(280<32(+(&+(38>B&&)&&)2(+(&)&))8((28()0&09&05&05<&"
QUAL_R="&<50&50&90&0)(82((8))&)&(+(2)&&)&&B>83(+&(+(23<082(A8((3+((+()()&:62(>+((((&DDDDFHJJJJJHHHHHFFFFFCCC"

printf "@read\n${SEQ}\n+\n${QUAL}\n" > fwd.fq
printf "@read\n${SEQ_RC}\n+\n${QUAL_R}\n" > rev.fq

vg map -f fwd.fq -g x.gcsa -x x.xg > mapped.fwd.gam
vg map -f rev.fq -g x.gcsa -x x.xg > mapped.rev.gam

is "$(vg view -aj mapped.rev.gam | jq -r '.quality' | base64 -d | xxd -p -c1 | tac | xxd -p -r | xxd)" "$(vg view -aj mapped.fwd.gam | jq -r '.quality' | base64 -d | xxd)" "quality strings we will use for testing are oriented correctly"

is "$(vg surject -p x -x x.xg mapped.fwd.gam -s | cut -f1,3,4,5,6,7,8,9,10,11)" "$(vg surject -p x -x x.xg mapped.rev.gam -s | cut -f1,3,4,5,6,7,8,9,10,11)" "forward and reverse orientations of a read produce the same surjected SAM, ignoring flags"

rm -f fwd.fq rev.fq mapped.fwd.gam mapped.rev.gam

is $(vg map -G <(vg sim -a -n 100 -x x.xg) -g x.gcsa -x x.xg | vg surject -p x -x x.xg -b - | samtools view - | wc -l) \
    100 "vg surject produces valid BAM output"

#is $(vg map -G <(vg sim -a -n 100 x.vg) x.vg | vg surject -p x -g x.gcsa -x x.xg -c - | samtools view - | wc -l) \
#    100 "vg surject produces valid CRAM output"

echo '{"sequence": "CAAATAA", "path": {"mapping": [{"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}]}, "mapping_quality": 99}' | vg view -JGa - > read.gam
is "$(vg surject -p x -x x.xg read.gam | vg view -aj - | jq '.mapping_quality')" "99" "mapping quality is preserved through surjection"

echo '{"name": "read/2", "sequence": "CAAATAA", "path": {"mapping": [{"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}]}, "fragment_prev": {"name": "read/1"}}{"name": "read/1", "sequence": "CTTATTT", "path": {"mapping": [{"position": {"node_id": 1, "is_reverse": true}, "edit": [{"from_length": 7, "to_length": 7}]}]}, "fragment_next": {"name": "read/2"}}' | vg view -JGa - > read.gam
is "$(vg surject -p x -x x.xg -i read.gam | vg view -aj - | jq -r 'select(.name == "read/2") | .fragment_prev.name')" "read/1" "read pairing is preserved through GAM->GAM surjection"

vg surject -p x -x x.xg -i read.gam -s > read.gam.surject.sam
vg convert x.xg -G read.gam -t 1 | vg surject -p x -x x.xg -i -G - -s > read.gaf.surject.sam
diff read.gam.surject.sam read.gaf.surject.sam
is $? 0 "interleaved surjection produces same SAM when using GAF and GAM inputs"
rm -f read.gam.surject.sam read.gaf.surject.sam

vg map -d x -iG <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) | vg surject -x x.xg -p x -s -i -N Sample1 -R RG1 - >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort | cut -f 4)" "$(printf '321\n762')" "surjection of paired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort | cut -f 8)" "$(printf '762\n321')" "surjection of paired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "1" "surjection of paired reads to SAM yields properly matched QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '=\n=')" "surjection of paired reads to SAM produces correct pair partner contigs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 2 | sort -n)" "$(printf '83\n163')" "surjection of paired reads to SAM produces correct flags"
is "$(cat surjected.sam | grep -v '^@' | grep 'RG1' | wc -l)" "2" "surjection of paired reads to SAM tags both reads with a read group"
is "$(cat surjected.sam | grep '@RG' | grep 'RG1' | grep 'Sample1' | wc -l)" "1" "surjection of paired reads to SAM creates RG header"

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
is $(vg map -f minigiab/NA12878.chr22.tiny.fq.gz -x m.xg -g m.gcsa --gaf | vg surject -p q -x m.xg -s - -G | grep chr22.bin8.cram:166:6027 | grep BBBBBFBFI | wc -l) 1 "mapping reproduces qualities from GAF input"

is "$(zcat < minigiab/NA12878.chr22.tiny.fq.gz | head -n 4000 | vg mpmap -B -p -x m.xg -g m.gcsa -M 1 -f - | vg surject -m -x m.xg -p q -s - | samtools view | wc -l)" 1000 "surject works on GAMP input"

is "$(vg sim -x m.xg -n 500 -l 150 -a -s 768594 -i 0.01 -e 0.01 -p 250 -v 50 | vg view -aX - | vg mpmap -B -p -b 200 -x m.xg -g m.gcsa -i -M 1 -f - | vg surject -m -x m.xg -i -p q -s - | samtools view | wc -l)" 1000 "surject works on paired GAMP input"

rm -rf minigiab.vg* m.xg m.gcsa

vg construct -r small/x.fa >j.vg
vg index j.vg -g j.gcsa
vg map -x j.vg -g j.gcsa -s TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA > r.gam
vg surject -x j.vg r.gam -s > r.sam

cat small/x.fa | sed -e 's/x/x[500]/g' > x.sub.fa
vg construct -r x.sub.fa >j.sub.vg
vg index j.sub.vg -g j.sub.gcsa
vg map -x j.sub.vg -g j.sub.gcsa -s TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA > r.sub.gam
vg surject -x j.sub.vg r.sub.gam -s > r.sub.sam

cat r.sam | sed -e 's/LN:1001/LN:1501/g' -e 's/161/661/g' > r.manual.sam
diff r.manual.sam r.sub.sam
is "$?" 0 "vg surject correctly handles subpath suffix in path name"

printf "x\t2000\n" > path_info.tsv
rm -f r.sub.sam r.manual.sam
vg surject -x j.sub.vg r.sub.gam -s --ref-paths path_info.tsv > r.sub.sam
cat r.sam | sed -e 's/LN:1001/LN:2000/g' -e 's/161/661/g' > r.manual.sam
diff r.manual.sam r.sub.sam
is "$?" 0 "vg surject correctly fetches base path length from input file"

rm -f h.vg h.gcsa r.gam r.sam x.sub.fa j.sub.vg j.sub.gcsa r.sub.gam r.sub.sam r.sub.sam

vg surject -s -x surject/perpendicular.vg surject/perpendicular.gam > perpendicular.sam
is "$?" 0 "vg surject does not crash when surjecting a read that grazes the reference with a deletion"
is "$(cat perpendicular.sam | grep -v "^@" | cut -f2)" "4" "vg surject leaves a read that grazes the reference with a deletion unmapped"

vg surject -s --prune-low-cplx -x surject/perpendicular.vg surject/perpendicular.gam > perpendicular.sam
is "$?" 0 "vg surject does not crash when surjecting a read that grazes the reference with a deletion and pruning low complexity anchors"
is "$(cat perpendicular.sam | grep -v "^@" | cut -f2)" "4" "vg surject leaves a read that grazes the reference with a deletion unmapped when pruning low complexity anchors"

rm -f perpendicular.sam

vg construct -r small/x.fa > x.vg
cat <(vg view x.vg) <(vg view x.vg | grep P | sed 's/P\tx/P\ty/') | vg convert -g - > x.pathdup.vg
vg index -x x.xg -g x.gcsa x.pathdup.vg
vg sim -x x.xg -n 20 -l 40 -p 60 -v 10 -a > x.gam
vg mpmap -x x.xg -g x.gcsa -n dna --suppress-mismapping -B -G x.gam -i -F GAM -I 60 -D 10 > mapped.gam
vg mpmap -x x.xg -g x.gcsa -n dna --suppress-mismapping -B -G x.gam -i -F GAMP -I 60 -D 10 > mapped.gamp

is "$(vg surject -x x.xg -s -t 1 mapped.gam | grep -v '@' | wc -l)" 40 "GAM surject can return only primaries"
is "$(vg surject -x x.xg -M -s -t 1 mapped.gam | grep -v '@' | wc -l)" 80 "GAM surject can return multimappings"
is "$(vg surject -x x.xg -s -m -t 1 mapped.gamp | grep -v '@' | wc -l)" 40 "GAMP surject can return only primaries"
is "$(vg surject -x x.xg -M -m -s -t 1 mapped.gamp | grep -v '@' | wc -l)" 80 "GAMP surject can return multimappings"

rm x.vg rm x.pathdup.vg x.xg x.gcsa x.gcsa.lcp x.gam mapped.gam mapped.gamp


#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 117

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


head -c -10 j.gam >j-truncated.gam
vg surject -p x -x x.xg -t 1 j-truncated.gam >/dev/null 2>log.txt
is "${?}" "1" "vg surject stops when the input read file is truncated"
is "$(grep "truncated input" log.txt | wc -l)" "1" "vg surject reports that files are truncated"
rm -f j-truncated.gam log.txt

is $(vg surject -x x.xg -t 1 -s x.gam | grep AS | wc -l) 100 "vg surject reports alignment scores"

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

is "$(vg filter --tsv-out quality mapped.rev.gam | tac -rs 'x\|[^x]' | head -2 | tail -1)" "$(vg filter --tsv-out quality mapped.fwd.gam | tail -1)" "quality strings we will use for testing are oriented correctly"

is "$(vg surject -p x -x x.xg mapped.fwd.gam -s | cut -f1,3,4,5,6,7,8,9,10,11)" "$(vg surject -p x -x x.xg mapped.rev.gam -s | cut -f1,3,4,5,6,7,8,9,10,11)" "forward and reverse orientations of a read produce the same surjected SAM, ignoring flags"

rm -f fwd.fq rev.fq mapped.fwd.gam mapped.rev.gam

is "$(vg map -G <(vg sim -a -n 100 -x x.xg) -g x.gcsa -x x.xg | vg surject -p x -x x.xg -b - | samtools view - | wc -l)" \
    "100" "vg surject produces valid BAM output"

is "$(vg map -G <(vg sim -a -n 100 -x x.vg) -g x.gcsa -x x.vg | vg surject -p x -x x.xg -c - | samtools view - | wc -l)" \
    "100" "vg surject produces valid CRAM output"

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

# a uniform random sequence
printf "@read TG:Z:val\nGGCGACGTACTAGGGACTACAGTCCTTCGTCTTTCTCTCTCGACTCCGAA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n" > x.fq
vg map -d x -t 1 -f x.fq --comments-as-tags >x.gam
is $(vg surject -p x -x x.xg --bam-output x.gam | samtools view -f 4 | grep "TG:Z:val" | wc -l | sed 's/^[[:space:]]*//') 1 "Tags are preserved on unmapped reads"
vg map -d x -t 1 -f x.fq --comments-as-tags --gaf >x.gaf
is $(vg surject -p x -x x.xg --bam-output --gaf-input x.gaf | samtools view -f 4 | grep "TG:Z:val" | wc -l | sed 's/^[[:space:]]*//') 1 "Tags are preserved on unmapped reads in GAF"

# A sequence that should map
printf "@read TG:Z:val\nACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n" > x2.fq
vg map -d x -t 1 -f x2.fq --comments-as-tags --gaf >x2.gaf
is $(vg surject -p x -x x.xg --bam-output --gaf-input x2.gaf | samtools view | grep "TG:Z:val" | wc -l | sed 's/^[[:space:]]*//') 1 "Tags are preserved on mapped reads in GAF"

rm -rf j.vg x.vg j.gam x.gam x.gaf x.idx j.xg x.xg x.gcsa read.gam reads.gam surjected.sam x.fq x2.fq x2.gaf

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

cat r.sam | sed -e 's/LN:1001/LN:1501/g' -e 's/161/661/g' -e 's/.M5:[a-zA-Z0-9]*//g' > r.manual.sam
diff r.manual.sam r.sub.sam
is "$?" 0 "vg surject correctly handles subpath suffix in path name"

printf "x\t2000\n" > path_info.tsv
rm -f r.sub.sam r.manual.sam
vg surject -x j.sub.vg r.sub.gam -s --ref-paths path_info.tsv > r.sub.sam
cat r.sam | sed -e 's/LN:1001/LN:2000/g' -e 's/161/661/g' -e 's/.M5:[a-zA-Z0-9]*//g' > r.manual.sam
diff r.manual.sam r.sub.sam
is "$?" 0 "vg surject correctly fetches base path length from input file"

is "$(vg surject -x j.vg -b --graph-aln r.gam | samtools view | grep 'GR:Z:' | wc -l | sed 's/^[[:space:]]*//')" "1" "BAMs can be annotated with the graph-space alignment"

rm -f h.vg h.gcsa r.gam r.sam x.sub.fa j.vg j.gcsa j.gcsa.lcp j.sub.vg j.sub.gcsa j.sub.gcsa.lcp r.sub.gam r.sub.sam r.sub.sam path_info.tsv r.manual.sam

vg surject -U -s -x surject/perpendicular.vg surject/perpendicular.gam > perpendicular.sam
is "$?" 0 "vg surject does not crash when surjecting a read that grazes the reference with a deletion"
is "$(cat perpendicular.sam | grep -v "^@" | cut -f2)" "4" "vg surject leaves a read that grazes the reference with a deletion unmapped"

vg surject -U -s --prune-low-cplx -x surject/perpendicular.vg surject/perpendicular.gam > perpendicular.sam
is "$?" 0 "vg surject does not crash when surjecting a read that grazes the reference with a deletion and pruning low complexity anchors"
is "$(cat perpendicular.sam | grep -v "^@" | cut -f2)" "4" "vg surject leaves a read that grazes the reference with a deletion unmapped when pruning low complexity anchors"

rm -f perpendicular.sam

vg construct -r small/x.fa > x.vg
cat <(vg view x.vg) <(vg view x.vg | grep P | sed 's/P\tx/P\ty/') | vg convert -g - > x.pathdup.vg
vg index -x x.xg -g x.gcsa x.pathdup.vg
vg sim -x x.xg -n 20 -l 40 -p 60 -v 10 -a --random-seed 123 > x.gam
vg mpmap -x x.xg -g x.gcsa -n dna --suppress-mismapping -B -G x.gam -i -F GAM -I 60 -D 10 -t 1 > mapped.gam
vg mpmap -x x.xg -g x.gcsa -n dna --suppress-mismapping -B -G x.gam -i -F GAMP -I 60 -D 10 -t 1 > mapped.gamp

vg surject -x x.xg -s -t 1 mapped.gam >/dev/null 2>/dev/null
is "$?" "1" "GAM surject produces an error when given paired reads when not in interleaved mode"

vg surject -x x.xg -m -s -t 1 mapped.gamp >/dev/null 2>/dev/null
is "$?" "1" "GAMP surject produces an error when given paired reads when not in interleaved mode"

vg convert x.xg -G mapped.gam -t 1 >mapped.gaf
vg surject -x x.xg -G -s -t 1 mapped.gaf >/dev/null 2>/dev/null
is "$?" 1 "GAF surject produces an error when given paired reads when not in interleaved mode"

is "$(vg surject -x x.xg -U -s -t 1 mapped.gam | grep -v '@' | wc -l)" 40 "GAM surject can return only primaries"
is "$(vg surject -x x.xg -M -U -s -t 1 mapped.gam | grep -v '@' | wc -l)" 80 "GAM surject can return multimappings"
is "$(vg surject -x x.xg -M -i -s -t 1 mapped.gam | grep -v '@' | wc -l)" 80 "GAM surject can return paired multimappings"
is "$(vg surject -x x.xg -U -s -m -t 1 mapped.gamp | grep -v '@' | wc -l)" 40 "GAMP surject can return only primaries"
is "$(vg surject -x x.xg -M -U -m -s -t 1 mapped.gamp | grep -v '@' | wc -l)" 80 "GAMP surject can return multimappings"
is "$(vg surject -x x.xg -M -i -m -s -i -t 1 mapped.gamp | grep -v '@' | wc -l)" 80 "GAMP surject can return paired multimappings"

vg construct -r tiny/tiny.fa > tiny.vg
vg surject -x tiny.vg -s -t 1 mapped.gam >/dev/null 2>err.txt
is "${?}" "1" "Surjection fails when using the wrong graph for GAM"
is "$(cat err.txt | grep 'cannot be interpreted' | wc -l)" "1" "Surjection of GAM to the wrong graph reports the problem"
vg surject -x tiny.vg -s -t 1 -m mapped.gamp >/dev/null 2>err.txt
is "${?}" "1" "Surjection fails when using the wrong graph for GAMP"
cat err.txt 1>&2
is "$(cat err.txt | grep 'cannot be interpreted' | wc -l)" "1" "Surjection of GAMP to the wrong graph reports the problem"

rm x.vg x.pathdup.vg x.xg x.gcsa x.gcsa.lcp x.gam mapped.gam mapped.gamp tiny.vg err.txt

is "$(vg surject -p CHM13#0#chr8 -x surject/opposite_strands.gfa --prune-low-cplx --sam-output --gaf-input surject/opposite_strands.gaf | grep -v "^@" | cut -f3-12 | sort | uniq | wc -l)" 1 "vg surject low compelxity pruning gets the same alignment regardless of read orientation"

is "$(vg surject -p CHM13#0#chr8 -x surject/opposite_strands.gfa --read-length long --sam-output --gaf-input surject/opposite_strands.gaf)" "$(vg surject -p CHM13#0#chr8 -x surject/opposite_strands.gfa --prune-low-cplx --sam-output --gaf-input surject/opposite_strands.gaf)" "vg surject long read preset uses low-complexity pruning"

vg autoindex -p d -w map -g graphs/long_deletion.gfa
printf "@read\nGGGAGAGAGAGAGA\n+\nHHHHHHHHHHHHHH\n" > d.fq
vg map -d d -f d.fq > d.gam
vg surject -u -b -x d.xg d.gam > d.bam
is "$(samtools view -f 2048 d.bam | wc -l)" "1" "Supplementary alignments can be produced"
is "$(samtools view -F 2048 d.bam | grep "SA:Z:" | wc -l)" "1" "Primary alignments get the SA tag for supplementaries"
vg view -ak d.gam > d.gamp
vg surject -u -b -x d.xg -m d.gamp > d2.bam
is "$(samtools view -f 2048 d2.bam | wc -l)" "1" "Supplementary alignments can be produced with GAMP input"
is "$(samtools view -F 2048 d2.bam | grep "SA:Z:" | wc -l)" "1" "Primary alignments get the SA tag for supplementaries with GAMP input"
printf "@read\nTTTCTCTCTCTCTC\n+\nHHHHHHHHHHHHHH\n" > e.fq
vg map -d d -f d.fq -f e.fq > e.gam
vg surject -u -b -x d.xg -i e.gam > e.bam
is "$(samtools view -f 2048 e.bam | wc -l)" "2" "Paired supplementary alignments can be produced"
is $(samtools view -f 2048 e.bam | awk '{if ($7 == "=" || $7 == x) {print $0}}' | wc -l) "2" "Paired supplementary alignments have correct mate contig"
is $(samtools view -f 2112 e.bam | awk '{print $8}') $(samtools view -F 2048 -f 128 e.bam | awk '{print $4}') "Read 1 of paired supplementary alignments has correct mate position"
is $(samtools view -f 2176 e.bam | awk '{print $8}') $(samtools view -F 2048 -f 64 e.bam | awk '{print $4}') "Read 2 of paired supplementary alignments has correct mate position"
vg view -ak e.gam > e.gamp
vg surject -u -m -b -x d.xg -i e.gamp > e2.bam
is "$(samtools view -f 2048 e2.bam | wc -l)" "2" "Paired supplementary alignments can be produced with GAMP input"
is $(samtools view -f 2048 e2.bam | awk '{if ($7 == "=" || $7 == x) {print $0}}' | wc -l) "2" "Paired supplementary alignments have correct mate contig with GAMP input"
is $(samtools view -f 2112 e2.bam | awk '{print $8}') $(samtools view -F 2048 -f 128 e2.bam | awk '{print $4}') "Read 1 of paired supplementary alignments has correct mate position with GAMP input"
is $(samtools view -f 2176 e2.bam | awk '{print $8}') $(samtools view -F 2048 -f 64 e2.bam | awk '{print $4}') "Read 2 of paired supplementary alignments has correct mate position with GAMP input"

vg autoindex -p f -w map -g graphs/long_inversion.gfa
vg map -d f -f d.fq > f.gam
vg surject -u -b -x f.xg f.gam > f.bam
is "$(samtools view -f 2048 f.bam | wc -l)" "1" "Supplementary alignment can be produced with an inversion"
is "$(samtools view f.bam | cut -f 3 | uniq | wc -l)" "1" "Inversion supplementary is on the same contig"
is "$(samtools view -F 16 f.bam | wc -l)" "1" "One of inverted primary/supplementary pair is on forward strand"
is "$(samtools view -f 16 f.bam | wc -l)" "1" "One of inverted primary/supplementary pair is on reverse strand"

rm d.xg d.gcsa d.gcsa.lcp d.fq d.gam d.gamp d.bam d2.bam e.fq e.gam e.gamp e.bam e2.bam f.xg f.gcsa f.gcsa.lcp f.gam f.bam

vg gbwt -G graphs/haplotypes.gfa -g haplotypes.gbz --gbz-format
vg view -JGa reads/haplotypes_read.json >read.gam
vg surject -x haplotypes.gbz -p 'KOLF2.1J#1#chr1_1#0' --sam-output read.gam >surjected.sam

is "$(cat surjected.sam | tail -n1 | cut -f3)" "KOLF2.1J#1#chr1_1#0" "surjecting explicitly to a haplotype in a GBZ puts a read on that haplotype"

rm haplotypes.gbz read.gam surjected.sam

vg autoindex -p g -w map -g graphs/long_insertion.gfa
vg map -d g -f reads/ts.fq | vg surject -x g.xg -b --off-ref-position - > g.bam
is $(samtools view g.bam | grep "NR:Z:x:8+" | wc -l | sed 's/^[[:space:]]*//') "1" "off reference reads can be annotated with the nearest reference position"

rm g.xg g.gcsa g.gcsa.lcp g.bam

# The reference follows 1 -> 2 -> 3 -> 4 -> 5; the read follows the shortcut 1 -> 6 -> 5.
# Node 1 (12 bp) and node 3 share the same sequence; nodes 4 and 6 share the same 8-base
# sequence; node 2 is an 80-base run of C's separating them; node 5 is the 32-base body
# of the read. The read's left_tail_length=12 marks node 1 as the tail. Pruning that
# anchor lets the surjector place the whole read through 3 -> 4 -> 5 starting at
# SAM position 93. Without pruning, node 1 (ref pos 1) and node 5 (ref pos 113) are
# too far apart to emit as one alignment, so -u splits the tail into a supplementary record.
vg surject -x surject/tail-pruning.gfa -p ref -t 1 -u -s surject/tail-pruning.gam > tail-baseline.sam || exit 1
is "$(grep -v '^@' tail-baseline.sam | cut -f2-4,6 | sort -n)" "$(printf '0\tref\t105\t12S40M\n2048\tref\t1\t12M40S')" \
    "Without tail pruning, the misplaced tail is emitted separately"
vg surject -x surject/tail-pruning.gfa -p ref -t 1 -u -s --prune-tail-region surject/tail-pruning.gam > tail-pruned.sam || exit 1
is "$(grep -v '^@' tail-pruned.sam | cut -f2-4,6)" "$(printf '0\tref\t93\t52M')" \
    "Tail pruning places the entire read at reference position 93"

vg surject -x surject/tail-pruning.gfa -p ref -t 1 -u -s -G --prune-tail-region surject/tail-pruning.gaf > tail-gaf.sam || exit 1
is "$(grep -v '^@' tail-gaf.sam | cut -f2-4,6)" "$(printf '0\tref\t93\t52M')" \
    "GAF tail annotations produce the expected pruned placement"

rm tail-baseline.sam tail-pruned.sam tail-gaf.sam

# The haplotypes share flanks but differ by a substitution and six inserted bases.
# Each matching read should prefer its own haplotype over the indel alignment.
vg surject -x surject/diploid-map.gfa --diploid-map sample -s -t 1 surject/diploid-map.gam > diploid-selection.sam || exit 1
# Compare flags, path, position, MAPQ, CIGAR, and alignment score for each candidate.
awk 'BEGIN {OFS="\t"} !/^@/ {score=""; for(i=12;i<=NF;i++) if($i ~ /^AS:i:/) score=substr($i,6); print $1,$2,$3,$4,$5,$6,score}' diploid-selection.sam > diploid-selection.tsv
is "$(awk '$1 == "hap1_forward"' diploid-selection.tsv | cut -f2-)" \
    "$(printf '0\tsample#1#chr1\t1\t37\t52M\t62\n256\tsample#2#chr1\t1\t37\t32M6D20M\t46')" \
    "Haplotype 1 read prefers the exact match over a substitution and deletion"
is "$(awk '$1 == "hap2_forward"' diploid-selection.tsv | cut -f2-)" \
    "$(printf '0\tsample#2#chr1\t1\t37\t58M\t68\n256\tsample#1#chr1\t1\t37\t32M6I20M\t46')" \
    "Haplotype 2 read prefers the exact match over a substitution and insertion"
is "$(awk '$1 == "hap1_reverse"' diploid-selection.tsv | cut -f2-)" \
    "$(printf '16\tsample#1#chr1\t1\t37\t52M\t62\n272\tsample#2#chr1\t1\t37\t32M6D20M\t46')" \
    "Reverse-complement input selects the same haplotype with reverse-strand flags"
# This group's input primary is on node 5, which belongs to neither target path.
is "$(awk '$1 == "secondary_wins"' diploid-selection.tsv | cut -f2-)" \
    "$(printf '0\tsample#1#chr1\t1\t37\t52M\t62\n256\tsample#2#chr1\t1\t37\t32M6D20M\t46')" \
    "A secondary input supplies the winning placement when the primary is off target"

# Repeat these cases with distinct read names to exercise parallel grouped output.
vg surject -x surject/diploid-map.gfa --diploid-map sample -b -t 4 surject/diploid-map-parallel.gam > diploid-output.bam || exit 1
samtools view diploid-output.bam > diploid-output.sam
is "$(wc -l < diploid-output.sam)" 2050 "Every diploid read emits both target-path candidates"
awk 'NR % 2 == 1 {name=$1; flag=$2; if ((flag != 0 && flag != 16) || seen[name]++) exit 1} NR % 2 == 0 {if ($1 != name || $2 != flag+256) exit 1}' diploid-output.sam
is "$?" 0 "Multithreaded output keeps each primary and secondary together"
# Haplotype quality reaches 60 for these score differences; input MAPQ caps output at 37.
# The stale aq:i:1 is replaced by the input primary's MAPQ; ZZ:Z:keep is preserved.
awk '$5 != 37 {exit 1} {hp=0; hq=0; aq=0; zz=0; expected=(NR % 2 ? "hp:Z:pri_hap" : "hp:Z:sec_hap"); for(i=12;i<=NF;i++){if($i==expected)hp++; if($i=="hq:i:60")hq++; if($i=="aq:i:37")aq++; if($i=="ZZ:Z:keep")zz++} if(hp!=1 || hq!=1 || aq!=1 || zz!=1)exit 1}' diploid-output.sam
is "$?" 0 "Diploid output has capped MAPQ and preferred/alternative haplotype tags"
# Inputs with no diploid candidate are emitted as unmapped records.
vg surject -x surject/diploid-map.gfa -d sample -p 'sample#1#chr1' -b -t 1 surject/diploid-map-unmapped.gam > diploid-unmapped.bam
is "$?" 0 "Diploid reads with no candidate can be written to BAM"
samtools quickcheck diploid-unmapped.bam
is "$?" 0 "Unmapped diploid BAM is complete and readable"
is "$(samtools view diploid-unmapped.bam | cut -f1-6)" "$(printf 'empty_path\t4\t*\t0\t0\t*\noff_target\t4\t*\t0\t0\t*')" \
    "Empty paths and off-target placements produce unmapped BAM records"
rm diploid-selection.sam diploid-selection.tsv diploid-output.bam diploid-output.sam diploid-unmapped.bam

# Paired diploid selection reuses the unit-test layout: two separate paths,
# equal read scores, but fragment lengths 100 and 200. The model should select 200.
paired_dir=$(mktemp -d) || exit 1
paired_sequence=$(printf 'ACGT%.0s' {1..100})
printf 'H\tVN:Z:1.1\tRS:Z:sample\nS\t1\t%s\nS\t2\t%s\nW\tsample\t1\tchr1\t0\t400\t>1\nW\tsample\t2\tchr1\t0\t400\t>2\n' "$paired_sequence" "$paired_sequence" > "$paired_dir/graph.gfa"
jq -cn 'range(0;1025) as $i | (1,2) as $node | (1,2) as $mate |
    ("pair_" + ($i|tostring)) as $name |
    {name:($name + "/" + ($mate|tostring)), sequence:"ACGTACGTACGTACGTACGT",
     mapping_quality:(if $mate == 1 then 90 else 17 end), is_secondary:($node == 2),
     path:{mapping:[{position:{node_id:$node,is_reverse:($mate == 2),
          offset:(if $mate == 1 then 0 elif $node == 1 then 300 else 200 end)},
          edit:[{from_length:20,to_length:20}]}]}}
    + (if $mate == 1 then {fragment_next:{name:($name + "/2")}}
       else {fragment_prev:{name:($name + "/1")}} end)' > "$paired_dir/reads.json"
vg view -JGa "$paired_dir/reads.json" > "$paired_dir/reads.gam"
vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-mean 200 --fragment-stdev 10 -b -t 4 "$paired_dir/reads.gam" > "$paired_dir/four.bam"
is "$?" 0 "Paired diploid GAM runs with a fixed fragment model"
samtools quickcheck "$paired_dir/four.bam"
is "$?" 0 "Grouped paired BAM is complete"
samtools view "$paired_dir/four.bam" > "$paired_dir/four.records"
is "$(wc -l < "$paired_dir/four.records")" 4100 "Every fragment emits both complete pair alternatives"
awk 'NR%4==1 {name=$1; if(seen[name]++ || $2!=99 || $3!="sample#2#chr1" || $8!=181 || $9!=200)exit 1}
     NR%4==2 {if($1!=name || $2!=147 || $8!=1 || $9!=-200)exit 1}
     NR%4==3 {if($1!=name || $2!=355 || $3!="sample#1#chr1" || $9!=100)exit 1}
     NR%4==0 {if($1!=name || $2!=403 || $9!=-100)exit 1}' "$paired_dir/four.records"
is "$?" 0 "Fragment scoring selects the longer pair and keeps each group contiguous with correct flags and TLEN"
awk 'NR%2==1 {if($5!=60)exit 1} NR%2==0 {if($5!=17)exit 1}
     {expected=(NR%2==1 ? "aq:i:90" : "aq:i:17"); found=0; for(i=12;i<=NF;i++)if($i==expected)found++; if(found!=1)exit 1}' "$paired_dir/four.records"
is "$?" 0 "Each mate retains its own original MAPQ cap and aq tag"
vg convert "$paired_dir/graph.gfa" -G "$paired_dir/reads.gam" -t 1 > "$paired_dir/reads.gaf"
vg surject -x "$paired_dir/graph.gfa" -d sample -i -G --fragment-mean 200 --fragment-stdev 10 -s -t 1 "$paired_dir/reads.gaf" > "$paired_dir/one.sam"
is "$?" 0 "Paired grouped GAF accepts the same input contract"
cmp -s <(grep -v '^@' "$paired_dir/one.sam" | sort) <(sort "$paired_dir/four.records")
is "$?" 0 "Paired GAM/BAM and GAF/SAM agree across thread counts"
vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-mean 200 --fragment-stdev 10 -t 4 "$paired_dir/reads.gam" > "$paired_dir/output.gam"
vg view -aj "$paired_dir/output.gam" | jq -se 'length==4100 and all(.[]; (.refpos|length)==1 and (.fragment_next.name // .fragment_prev.name)!=null)' > /dev/null
is "$?" 0 "Paired GAM output retains reference positions and reciprocal mate links"
head -n 2 "$paired_dir/reads.json" | vg view -JGa - > "$paired_dir/single.gam"
vg surject -x "$paired_dir/graph.gfa" -d sample -i -f 50 --fragment-mean 200 --fragment-stdev 10 -s "$paired_dir/single.gam" > "$paired_dir/fallback.sam"
is "$(grep -v '^@' "$paired_dir/fallback.sam" | cut -f2,5)" "$(printf '97\t0\n145\t0')" "No compatible pair uses the documented MAPQ-zero improper fallback"
# Learning uses confident, unambiguous pairs, and replays the buffered prefix.
jq -c 'select(.is_secondary==false) | .mapping_quality=60' "$paired_dir/reads.json" | vg view -JGa - > "$paired_dir/learn.gam"
vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-sample-size 2 -b -t 4 "$paired_dir/learn.gam" > "$paired_dir/learn.bam" 2> "$paired_dir/learn.log"
is "$?" 0 "Paired diploid learns a model before parallel output"
grep -q 'Learned diploid fragment mean 100, stdev 1 from 2 pairs' "$paired_dir/learn.log"
is "$?" 0 "Constant observed distances use a finite fragment standard deviation"
is "$(samtools view -c "$paired_dir/learn.bam")" 2050 "Learning replays buffered pairs without losing or duplicating reads"
vg surject -x "$paired_dir/graph.gfa" -d sample -i -b "$paired_dir/single.gam" > "$paired_dir/short.bam" 2> "$paired_dir/short.log"
grep -q 'using alignment scores only' "$paired_dir/short.log"
is "$?" 0 "EOF with insufficient training data reports score-only fallback"
is "$(samtools view -c "$paired_dir/short.bam")" 2 "EOF flushes a partial learning buffer"
vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-buffer-size 2 -b -t 4 "$paired_dir/reads.gam" > "$paired_dir/ambiguous.bam" 2> "$paired_dir/ambiguous.log"
grep -q 'using alignment scores only' "$paired_dir/ambiguous.log"
is "$?" 0 "Ambiguous pairs do not train the model and the buffer is bounded"
is "$(samtools view -c "$paired_dir/ambiguous.bam")" 4100 "Buffer-limit fallback preserves all pairs"
# Both empty and one-empty pairs must survive BAM serialization.
head -n 2 "$paired_dir/reads.json" | jq -c 'del(.path)' | vg view -JGa - > "$paired_dir/empty.gam"
vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-mean 200 --fragment-stdev 10 -b "$paired_dir/empty.gam" > "$paired_dir/empty.bam"
is "$(samtools view -f 12 -c "$paired_dir/empty.bam")" 2 "Fully unmapped pairs retain paired and mate-unmapped flags"
head -n 2 "$paired_dir/reads.json" | jq -c 'if .fragment_prev then del(.path) else . end' | vg view -JGa - > "$paired_dir/half.gam"
vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-mean 200 --fragment-stdev 10 -s "$paired_dir/half.gam" > "$paired_dir/half.sam"
is "$(grep -v '^@' "$paired_dir/half.sam" | cut -f2,9)" "$(printf '73\t0\n133\t0')" "One unmapped mate has correct mate flags and zero TLEN"
vg surject -d sample -i --fragment-mean 200 "$paired_dir/single.gam" > /dev/null 2> "$paired_dir/error"
is "$?" 1 "A fixed fragment model requires both parameters"
head -n 1 "$paired_dir/reads.json" | vg view -JGa - > "$paired_dir/odd.gam"
(ulimit -c 0; vg surject -x "$paired_dir/graph.gfa" -d sample -i --fragment-mean 200 --fragment-stdev 10 "$paired_dir/odd.gam" > /dev/null 2> "$paired_dir/error")
test "$?" -ne 0 && grep -q 'incomplete final pair' "$paired_dir/error"
is "$?" 0 "An incomplete interleaved pair is rejected instead of silently dropped"
rm -rf -- "$paired_dir"

# Reuse the split-anchor unit-test layout to exercise candidate-local mate links
# and atomic emission when each pair alternative also has a supplementary piece.
paired_dir=$(mktemp -d) || exit 1
jq -cn '{node:[{id:1,sequence:("A"*60)},{id:2,sequence:("G"*200)},{id:3,sequence:("C"*40)},
                   {id:4,sequence:("A"*60)},{id:5,sequence:("G"*200)},{id:6,sequence:("C"*40)}],
         edge:[{from:1,to:2},{from:2,to:3},{from:1,to:3},{from:4,to:5},{from:5,to:6},{from:4,to:6}],
         path:[{name:"sample#1#chr1",mapping:[{position:{node_id:1},rank:1},{position:{node_id:2},rank:2},{position:{node_id:3},rank:3}]},
               {name:"sample#2#chr1",mapping:[{position:{node_id:4},rank:1},{position:{node_id:5},rank:2},{position:{node_id:6},rank:3}]}]}' | vg view -Jv - > "$paired_dir/graph.vg"
jq -cn 'range(0;1025) as $i | (1,4) as $node | ("split_"+($i|tostring)) as $name |
    {name:($name+"/1"),sequence:(("A"*60)+("C"*40)),mapping_quality:60,is_secondary:($node==1),fragment_next:{name:($name+"/2")},
     path:{mapping:[{position:{node_id:$node},rank:1,edit:[{from_length:60,to_length:60}]},
                    {position:{node_id:($node+2)},rank:2,edit:[{from_length:40,to_length:40}]}]}},
    {name:($name+"/2"),sequence:("G"*40),mapping_quality:60,is_secondary:($node==1),fragment_prev:{name:($name+"/1")},
     path:{mapping:[{position:{node_id:($node+2),is_reverse:true},rank:1,edit:[{from_length:40,to_length:40}]}]}}' | vg view -JGa - > "$paired_dir/reads.gam"
vg surject -x "$paired_dir/graph.vg" -d sample -i -u --no-prune-low-cplx --fragment-mean 300 --fragment-stdev 10 -b -t 4 "$paired_dir/reads.gam" > "$paired_dir/split.bam"
is "$?" 0 "Paired diploid supplementary output can be written to BAM"
samtools view "$paired_dir/split.bam" > "$paired_dir/split.records"
is "$(wc -l < "$paired_dir/split.records")" 6150 "All paired alternatives retain their supplementary pieces"
awk 'NR%6==1 {name=$1; if(seen[name]++)exit 1}
     {part=(NR-1)%6; path=(part<3 ? "sample#1#chr1" : "sample#2#chr1");
      if($1!=name || $3!=path || $7!="=" || $8!=(part%3==1 ? 1 : 261))exit 1;
      if(part%3==2 && $9!=0)exit 1;
      if(part%3!=1){found=0; for(i=12;i<=NF;i++)if(index($i,"SA:Z:"path",")==1)found++; if(found!=1)exit 1}}' "$paired_dir/split.records"
is "$?" 0 "Supplementaries stay with their own pair and have candidate-local mate positions and SA tags"
rm -rf -- "$paired_dir"

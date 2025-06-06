#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 65

vg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg
vg gbwt -o x-haps.gbwt -x x.vg -v small/x.vcf.gz
vg gbwt -o x-paths.gbwt -x x.vg --index-paths
vg gbwt -o x-merged.gbwt -m x-haps.gbwt x-paths.gbwt
vg gbwt -o x.gbwt --augment-gbwt -x x.vg x-merged.gbwt
vg index -j x.dist x.vg
vg minimizer -k 29 -w 11 -d x.dist -g x.gbwt -o x.shortread.withzip.min -z x.shortread.zipcodes x.xg


# For later tests we expect this to make x.giraffe.gbz so we can't have an x.gbz around.
rm -f x.gbz x.giraffe.gbz
vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with xg + gbwt + min + zips + dist specified without crashing"
rm -f x-haps.gbwt x-paths.gbwt x-merged.gbwt

vg giraffe -Z x.giraffe.gbz -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with gbz + min + zips + dist specified without crashing"

rm -f x.giraffe.gbz
vg gbwt -x x.xg -g x.gg x.gbwt
vg giraffe -g x.gg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with gg + gbwt + min + zips + dist specified without crashing"

rm -f x.shortread.withzip.min
rm -f x.shortread.zipcodes
vg giraffe -Z x.giraffe.gbz -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with the minimizer index and zipcodes being regenerated"

rm -f x.shortread.withzip.min
vg giraffe -Z x.giraffe.gbz -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with the minimizer index regenerated"

rm -f x.shortread.zipcodes
vg giraffe -Z x.giraffe.gbz -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with the zipcodes being regenerated"

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with the indexes being inferred by name"

is "$(vg view -aj mapped1.gam | grep 'time_used' | wc -l | sed 's/^[[:space:]]*//')" "1" "Mapping logs runtime per read"

is "$(vg view -aj mapped1.gam | jq '.score')" "73" "Mapping produces the correct score"

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq -b fast >/dev/null
is "${?}" "0" "a read can be mapped with the fast preset"

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq -b default >/dev/null
is "${?}" "0" "a read can be mapped with the default preset"

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq -b chaining-sr >/dev/null
is "${?}" "0" "a read can be mapped with the short read chaining preset"

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq -f reads/small.middle.ref.fq -b chaining-sr >/dev/null 2>log.txt
is "${?}" "1" "a read pair cannot be mapped with the short read chaining preset"
is "$(cat log.txt | grep "not yet implemented" | wc -l)" "1" "trying to map paired-end data with chaining produces an informative error"
rm -f log.txt

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.mismatched.fq -b chaining-sr >/dev/null
is "${?}" "0" "a read with a mismatch can be mapped with the short read chaining preset"

rm -Rf grid-out
mkdir grid-out
vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq --output-basename grid-out/file --hard-hit-cap 5:6
is "$(ls grid-out/*.gam | wc -l)" "2" "Grid search works end-inclusive"
rm -Rf grid-out

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq --full-l-bonus 0 > mapped-nobonus.gam
is "$(vg view -aj  mapped-nobonus.gam | jq '.score')" "63" "Mapping without a full length bonus produces the correct score"
rm -f mapped-nobonus.gam

is "$(vg giraffe -Z x.giraffe.gbz -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f reads/small.middle.ref.fq -o BAM --add-graph-aln | samtools view | grep "GR:Z:" | wc -l | sed 's/^[[:space:]]*//')" "1" "BAMs can be annotated with graph alignment"

vg minimizer -k 29 -b -s 18 -d x.dist -g x.gbwt -o x.sync x.xg

vg giraffe -x x.xg -H x.gbwt -m x.sync -d x.dist -f reads/small.middle.ref.fq > mapped.sync.gam
is "${?}" "0" "a read can be mapped with syncmer indexes without crashing"

chmod 400 x.dist
vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f reads/small.middle.ref.fq >/dev/null 
is "${?}" "0" "a read can be mapped when the distance index is not writable"

echo "@read" >read.fq
echo "GATTACATTAGGAGATAGCCATACGACGTAGCATCTAGCTCAGCCACA$(cat small/x.fa | head -n2 | tail -n1)" >>read.fq
echo "+" >>read.fq
echo "GATTACATTAGGAGATAGCCATACGACGTAGCATCTAGCTCAGCCACA$(cat small/x.fa | head -n2 | tail -n1)" | tr 'ACGTN' '(((((' >>read.fq

vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f read.fq > read.gam
LOOP_LINES="$(vg view -aj read.gam | jq -c 'select(.path.mapping[0].position.node_id == .path.mapping[1].position.node_id)' | wc -l | sed 's/^[[:space:]]*//')"
is "${LOOP_LINES}" "0" "a read which softclips does not appear to loop"


printf "@read1\tT1:A:t T2:i:1\t T3:f:3.5e-7\nCACCGTGATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAG\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n" > tagged1.fq
printf "@read2 T4:Z:str T5:H:FF00\tT6:B:S,0,10 T7:B:f,8.0,5.0\nCACCGTGATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAG\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n" > tagged2.fq
vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f tagged1.fq --comments-as-tags -o BAM > t1.bam
vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f tagged2.fq --comments-as-tags -o BAM > t2.bam
vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f tagged1.fq -f tagged2.fq --comments-as-tags -o BAM > t3.bam
vg giraffe -x x.xg -H x.gbwt -m x.shortread.withzip.min -z x.shortread.zipcodes -d x.dist -f tagged1.fq --comments-as-tags -o GAF > t1.gaf


is "$(samtools view t1.bam | grep T1 | grep T2 | grep T3 | wc -l | sed 's/^[[:space:]]*//')" "1" "BAM tags are preserved on read 1"
is "$(samtools view t2.bam | grep T4 | grep T5 | grep T6 | wc -l | sed 's/^[[:space:]]*//')" "1" "BAM tags are preserved on read 2"
is "$(samtools view t3.bam | grep T1 | grep T2 | grep T3 | grep read1 | wc -l | sed 's/^[[:space:]]*//')" "1" "BAM tags are preserved on paired read 1"
is "$(samtools view t3.bam | grep T4 | grep T5 | grep T6 | grep read2 | wc -l | sed 's/^[[:space:]]*//')" "1" "BAM tags are preserved on paired read 2"
is "$(cat t1.gaf | grep T1 | grep T2 | grep T3 | wc -l | sed 's/^[[:space:]]*//')" "1" "GAF tags are preserved on read 1"


rm t1.bam t2.bam t3.bam t1.gaf tagged1.fq tagged2.fq
rm -f read.fq read.gam
rm -f x.vg x.xg x.gbwt x.shortread.zipcodes x.shortread.withzip.min x.sync x.dist x.gg
rm -f x.giraffe.gbz


cp small/x.fa .
cp small/x.vcf.gz .
cp small/x.vcf.gz.tbi .

vg giraffe x.fa x.vcf.gz -f reads/small.middle.ref.fq > mapped2.gam
is "${?}" "0" "a read can be mapped with just FASTA and VCF without crashing"

# These files can differ as serialized and still represent the same data, due to protobuf field order not being specified.
# Tripping through JSON will sort all the keys.
# Make sure to also remove time info
vg view -aj mapped1.gam | jq -c '.time_used = null' >mapped1.json
vg view -aj mapped2.gam | jq -c '.time_used = null' >mapped2.json
vg view -aj mapped.sync.gam | jq -c '.time_used = null' >mapped.sync.json

# Make sure at least one file converted successfully
SIZE="$(wc -c mapped2.json | cut -f1 -d' ')"
EMPTY=0
if [ "${SIZE}" == "0" ] ; then
    EMPTY=1
fi
is "${EMPTY}" "0" "mapping with just a FASTA and a VCF produced JSON-able alignments"

diff mapped1.json mapped2.json
is "${?}" "0" "mapping to manually-generated indexes and automatically-generated indexes is the same"

is "$(jq '.path' mapped1.json)" "$(jq '.path' mapped.sync.json)" "mapping with syncmers produces the same alignment as mapping with minimizers"

rm -rf mapped1.gam mapped1.json mapped2.gam mapped2.json mapped.sync.gam mapped.sync.json

vg giraffe x.fa x.vcf.gz -f small/x.fa_1.fastq > single.gam
is "$(vg view -aj single.gam | jq -c 'select((.fragment_next | not) and (.fragment_prev | not))' | wc -l | sed 's/^[[:space:]]*//')" "1000" "unpaired reads lack cross-references"

vg giraffe x.fa x.vcf.gz -f small/x.fa_1.fastq -f small/x.fa_1.fastq --fragment-mean 300 --fragment-stdev 100 > paired.gam
is "$(vg view -aj paired.gam | jq -c 'select((.fragment_next | not) and (.fragment_prev | not))' | wc -l | sed 's/^[[:space:]]*//')" "0" "paired reads have cross-references"

# Test paired surjected mapping
vg giraffe x.fa x.vcf.gz --show-work -iG <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --output-format SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of paired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '762\n321')" "surjection of paired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l | sed 's/^[[:space:]]*//')" "1" "surjection of paired reads to SAM yields properly matched QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '=\n=')" "surjection of paired reads to SAM produces correct pair partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '163\n83')" "surjection of paired reads to SAM produces correct flags"

# And unpaired surjected mapping
vg giraffe x.fa x.vcf.gz -G <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --output-format SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of unpaired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '0\n0')" "surjection of unpaired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l | sed 's/^[[:space:]]*//')" "2" "surjection of unpaired reads to SAM yields distinct QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '*\n*')" "surjection of unpaired reads to SAM produces absent partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '0\n16')" "surjection of unpaired reads to SAM produces correct flags"

rm -f x.vg x.gbwt x.xg x.min x.shortread.withzip.min x.shortread.zipcodes x.dist x.gg x.fa x.fa.fai x.vcf.gz x.vcf.gz.tbi single.gam paired.gam surjected.sam
rm -f x.giraffe.gbz

rm -f xy.vg xy.gbwt xy.xg xy.shortread.zipcodes xy.shortread.withzip.min xy.dist xy.gg xy.fa xy.fa.fai xy.vcf.gz xy.vcf.gz.tbi
cp small/xy.fa .
cp small/xy.vcf.gz .
cp small/xy.vcf.gz.tbi .
vg giraffe xy.fa xy.vcf.gz -f small/x.fa_1.fastq -o SAM --ref-paths small/yx.dict | sed 's/.M5:[a-zA-Z0-9]*//g' | grep -E "^@(SQ|HD)" > surjected-yx.dict
vg giraffe xy.fa xy.vcf.gz -f small/x.fa_1.fastq -o SAM --ref-paths small/xy.dict | sed 's/.M5:[a-zA-Z0-9]*//g' | grep -E "^@(SQ|HD)" > surjected-xy.dict

diff surjected-yx.dict small/yx.dict
is "${?}" "0" "surjecting with a sequence dictionary in non-sorted order produces headers in non-sorted order"

diff surjected-xy.dict small/xy.dict
is "${?}" "0" "surjecting with a sequence dictionary in sorted order produces headers in sorted order"

rm -f surjected-yx.dict surjected-xy.dict

vg construct -r small/x.fa > x.vg
vg sim -a -p 200 -v 10 -l 50 -n 1000 -s 12345 -x x.vg >x.gam

vg giraffe xy.fa xy.vcf.gz -G x.gam -o SAM -i --fragment-mean 200 --fragment-stdev 10 --distance-limit 50 >xy.sam
X_HITS="$(cat xy.sam | grep -v "^@" | cut -f3 | grep x | wc -l)"
if [ "${X_HITS}" -lt 1200 ] && [ "${X_HITS}" -gt 800 ] ; then
    IN_RANGE="1"
else
    IN_RANGE="0"
fi
is "${IN_RANGE}" "1" "paired reads are evenly split between equivalent mappings"

vg giraffe xy.fa xy.vcf.gz -G x.gam -o SAM >xy.sam
X_HITS="$(cat xy.sam | grep -v "^@" | cut -f3 | grep x | wc -l)"
if [ "${X_HITS}" -lt 1200 ] && [ "${X_HITS}" -gt 800 ] ; then
    IN_RANGE="1"
else
    IN_RANGE="0"
fi
is "${IN_RANGE}" "1" "unpaired reads are evenly split between equivalent mappings"

vg giraffe xy.fa xy.vcf.gz -G x.gam --track-provenance --discard
is $? "0" "provenance tracking succeeds for unpaired reads"

vg giraffe xy.fa xy.vcf.gz -G x.gam --track-provenance --track-correctness -o json >xy.json
is $? "0" "correctness tracking succeeds for unpaired reads"

is "$(cat xy.json | grep "correct-minimizer-coverage" | wc -l)" "2000" "unpaired reads are annotated with minimizer coverage"

vg giraffe xy.fa xy.vcf.gz -G x.gam -i --fragment-mean 200 --fragment-stdev 10 --distance-limit 50 --track-provenance --discard
is $? "0" "provenance tracking succeeds for paired reads"

vg giraffe xy.fa xy.vcf.gz -G x.gam -i --fragment-mean 200 --fragment-stdev 10 --distance-limit 50 --track-provenance --track-correctness -o json >xy.json
is $? "0" "correctness tracking succeeds for paired reads"

is "$(cat xy.json | grep "correct-minimizer-coverage" | wc -l | sed 's/^[[:space:]]*//')" "2000" "paired reads are annotated with minimizer coverage"
is "$(cat xy.json | jq '.annotation["correct-minimizer-coverage"] | select(. > 0)' | wc -l | sed 's/^[[:space:]]*//')" "2000" "paired reads all have nonzero correct minimizer coverage"

rm -f x.vg xy.sam xy.json
rm -f xy.vg xy.gbwt xy.xg xy.shortread.zipcodes xy.shortread.withzip.min xy.dist xy.gg xy.fa xy.fa.fai xy.vcf.gz xy.vcf.gz.tbi

vg giraffe -Z xy.giraffe.gbz -G x.gam -o BAM >xy.bam
is $? "0" "vg giraffe can map to BAM using a GBZ alone"
is "$(samtools view xy.bam | wc -l | sed 's/^[[:space:]]*//')" "2000" "GBZ-based BAM contains the expected number of reads"

rm -f x.gam xy.bam
rm -f xy.giraffe.gbz
rm -f xy.vg xy.gbwt xy.xg xy.shortread.zipcodes xy.shortread.withzip.min xy.dist xy.gg xy.fa xy.fa.fai xy.vcf.gz xy.vcf.gz.tbi

vg autoindex -p brca -w giraffe -g graphs/cactus-BRCA2.gfa 
vg sim -s 100 -x brca.giraffe.gbz -n 200 -a > reads.gam
vg giraffe -Z brca.giraffe.gbz -m brca.shortread.withzip.min -z brca.shortread.zipcodes -d brca.dist -G reads.gam --named-coordinates > mapped.gam
is "$?" "0" "Mapping reads to named coordinates on a nontrivial graph without walks succeeds"
is "$(vg view -aj mapped.gam | jq -r '.score' | grep -v "^0" | grep -v "null" | wc -l | sed 's/^[[:space:]]*//')" "200" "Reads are all mapped"
is "$(vg view -aj mapped.gam | jq -r '.path.mapping[].position.name' | grep null | wc -l | sed 's/^[[:space:]]*//')" "0" "GFA segment names are all set"

vg giraffe -Z brca.giraffe.gbz -m brca.shortread.withzip.min -z brca.shortread.zipcodes -d brca.dist -G reads.gam --named-coordinates -o gaf > mapped.gaf
is "$?" "0" "Mapping reads as GAF to named coordinates on a nontrivial graph without walks succeeds"

vg view -aj mapped.gam | jq -r '.path.mapping[].position.name' | sort | uniq > gam_names.txt
cat mapped.gaf | cut -f6 | tr '><' '\n\n' | grep "." | sort | uniq > gaf_names.txt

is "$(md5sum gaf_names.txt | cut -f1 -d' ')" "$(md5sum gam_names.txt | cut -f1 -d' ')" "Mapping reads as named GAF uses the same names as named GAM"

rm -f reads.gam mapped.gam mapped.gaf brca.* gam_names.txt gaf_names.txt

# Try long read alignment with Distance Index 2
vg construct -S -a -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >1mb1kgp.vg 2>/dev/null
vg index -j 1mb1kgp.dist  1mb1kgp.vg
vg autoindex -p 1mb1kgp -w giraffe -P "VG w/ Variant Paths:1mb1kgp.vg" -P "Giraffe Distance Index:1mb1kgp.dist" -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz
vg giraffe -Z 1mb1kgp.giraffe.gbz -f reads/1mb1kgp_longread.fq >longread.gam -U 300 --progress --track-provenance --align-from-chains --set-refpos
# This is an 8001 bp read with 1 insert and 1 substitution
# 7999 * 1 + 1 * -4 + -6 + 5 + 5 = 7999
is "$(vg view -aj longread.gam | jq -r '.score')" "7999" "A long read can be correctly aligned"
is "$(vg view -aj longread.gam | jq -c '.path.mapping[].edit[] | select(.sequence)' | wc -l | sed 's/^[[:space:]]*//')" "2" "A long read has the correct edits found"
is "$(vg view -aj longread.gam | jq -c '. | select(.annotation["filter_3_cluster-coverage_cluster_passed_size_total"] <= 300)' | wc -l | sed 's/^[[:space:]]*//')" "1" "Long read minimizer set is correctly restricted"
is "$(vg view -aj longread.gam | jq -c '.refpos[]' | wc -l)" "$(vg view -aj longread.gam | jq -c '.path.mapping[]' | wc -l)" "Giraffe sets refpos for each reference node"
is "$(vg view --extract-tag PARAMS_JSON longread.gam | jq '.["track-provenance"]')" "true" "Giraffe embeds parameters in GAM"

rm -f longread.gam 1mb1kgp.dist 1mb1kgp.giraffe.gbz 1mb1kgp.shortread.withzip.min 1mb1kgp.shortread.zipcodes log.txt


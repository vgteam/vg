#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 45

vg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg
vg gbwt -o x-haps.gbwt -x x.vg -v small/x.vcf.gz
vg gbwt -o x-paths.gbwt -x x.vg --index-paths
vg gbwt -o x-merged.gbwt -m x-haps.gbwt x-paths.gbwt
vg gbwt -o x.gbwt --augment-gbwt -x x.vg x-merged.gbwt
vg snarls --include-trivial x.vg > x.snarls
vg index -j x.dist x.vg
vg minimizer -k 29 -w 11 -g x.gbwt -o x.min x.xg

vg giraffe -x x.xg -H x.gbwt -m x.min -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with xg + gbwt + min + dist specified without crashing"
rm -f x-haps.gbwt x-paths.gbwt x-merged.gbwt

vg giraffe -Z x.giraffe.gbz -m x.min -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with gbz + min + dist specified without crashing"

rm -f x.giraffe.gbz
vg gbwt -x x.xg -g x.gg x.gbwt
vg giraffe -g x.gg -H x.gbwt -m x.min -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with gg + gbwt + min + dist specified without crashing"

rm -f x.min
vg giraffe -Z x.giraffe.gbz -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with the minimizer index being regenerated"

vg giraffe -Z x.giraffe.gbz -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with the indexes being inferred by name"

is "$(vg view -aj mapped1.gam | grep 'time_used' | wc -l)" "1" "Mapping logs runtime per read"

vg minimizer -k 29 -b -s 18 -g x.gbwt -o x.sync x.xg

vg giraffe -x x.xg -H x.gbwt -m x.sync -d x.dist -f reads/small.middle.ref.fq > mapped.sync.gam
is "${?}" "0" "a read can be mapped with syncmer indexes without crashing"

chmod 400 x.dist
vg giraffe -x x.xg -H x.gbwt -m x.min -d x.dist -f reads/small.middle.ref.fq >/dev/null 
is "${?}" "0" "a read can be mapped when the distance index is not writable"

rm -f x.vg x.xg x.gbwt x.snarls x.min x.sync x.dist x.gg
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
is "$(vg view -aj single.gam | jq -c 'select((.fragment_next | not) and (.fragment_prev | not))' | wc -l)" "1000" "unpaired reads lack cross-references"

vg giraffe x.fa x.vcf.gz -f small/x.fa_1.fastq -f small/x.fa_1.fastq --fragment-mean 300 --fragment-stdev 100 > paired.gam
is "$(vg view -aj paired.gam | jq -c 'select((.fragment_next | not) and (.fragment_prev | not))' | wc -l)" "0" "paired reads have cross-references"

# Test paired surjected mapping
vg giraffe x.fa x.vcf.gz -iG <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --output-format SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of paired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '762\n321')" "surjection of paired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "1" "surjection of paired reads to SAM yields properly matched QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '=\n=')" "surjection of paired reads to SAM produces correct pair partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '163\n83')" "surjection of paired reads to SAM produces correct flags"

# And unpaired surjected mapping
vg giraffe x.fa x.vcf.gz -G <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --output-format SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of unpaired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '0\n0')" "surjection of unpaired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "2" "surjection of unpaired reads to SAM yields distinct QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '*\n*')" "surjection of unpaired reads to SAM produces absent partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '0\n16')" "surjection of unpaired reads to SAM produces correct flags"

rm -f x.vg x.gbwt x.xg x.snarls x.min x.dist x.gg x.fa x.fa.fai x.vcf.gz x.vcf.gz.tbi single.gam paired.gam surjected.sam
rm -f x.giraffe.gbz

rm -f xy.vg xy.gbwt xy.xg xy.snarls xy.min xy.dist xy.gg xy.fa xy.fa.fai xy.vcf.gz xy.vcf.gz.tbi
cp small/xy.fa .
cp small/xy.vcf.gz .
cp small/xy.vcf.gz.tbi .
vg giraffe xy.fa xy.vcf.gz -f small/x.fa_1.fastq -o SAM --ref-paths small/yx.dict | grep -E "^@(SQ|HD)" > surjected-yx.dict
vg giraffe xy.fa xy.vcf.gz -f small/x.fa_1.fastq -o SAM --ref-paths small/xy.dict | grep -E "^@(SQ|HD)" > surjected-xy.dict

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

is "$(cat xy.json | grep "correct-minimizer-coverage" | wc -l)" "2000" "paired reads are annotated with minimizer coverage"
is "$(cat xy.json | jq '.annotation["correct-minimizer-coverage"] | select(. > 0)' | wc -l)" "2000" "paired reads all have nonzero correct minimizer coverage"

rm -f x.vg xy.sam xy.json
rm -f xy.vg xy.gbwt xy.xg xy.snarls xy.min xy.dist xy.gg xy.fa xy.fa.fai xy.vcf.gz xy.vcf.gz.tbi

vg giraffe -Z xy.giraffe.gbz -G x.gam -o BAM >xy.bam
is $? "0" "vg giraffe can map to BAM using a GBZ alone"
is "$(samtools view xy.bam | wc -l)" "2000" "GBZ-based BAM contains the expected number of reads"

rm -f x.gam xy.bam
rm -f xy.giraffe.gbz

vg autoindex -p brca -w giraffe -g graphs/cactus-BRCA2.gfa 
vg sim -s 100 -x brca.giraffe.gbz -n 200 -a > reads.gam
vg giraffe -Z brca.giraffe.gbz -m brca.min -d brca.dist -G reads.gam --named-coordinates > mapped.gam
is "$?" "0" "Mapping reads to named coordinates on a nontrivial graph without walks succeeds"
is "$(vg view -aj mapped.gam | jq -r '.score' | grep -v "^0" | wc -l)" "200" "Reads are all mapped"
is "$(vg view -aj mapped.gam | jq -r '.path.mapping[].position.name' | grep null | wc -l)" "0" "GFA segment names are all set"

vg giraffe -Z brca.giraffe.gbz -m brca.min -d brca.dist -G reads.gam --named-coordinates -o gaf > mapped.gaf
is "$?" "0" "Mapping reads as GAF to named coordinates on a nontrivial graph without walks succeeds"

vg view -aj mapped.gam | jq -r '.path.mapping[].position.name' | sort | uniq > gam_names.txt
cat mapped.gaf | cut -f6 | tr '><' '\n\n' | grep "." | sort | uniq > gaf_names.txt

is "$(md5sum gaf_names.txt | cut -f1 -d' ')" "$(md5sum gam_names.txt | cut -f1 -d' ')" "Mapping reads as named GAF uses the same names as named GAM"

rm -f reads.gam mapped.gam mapped.gaf brca.* gam_names.txt gaf_names.txt

# Try long read alignment with Distance Index 2
vg construct -S -a -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >1mb1kgp.vg 2>/dev/null
vg index -j 1mb1kgp.dist  1mb1kgp.vg
vg autoindex -p 1mb1kgp -w giraffe -P "VG w/ Variant Paths:1mb1kgp.vg" -P "Giraffe Distance Index:1mb1kgp.dist" -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz
vg giraffe -Z 1mb1kgp.giraffe.gbz -f reads/1mb1kgp_longread.fq >longread.gam -U 300 --progress --track-provenance --align-from-chains
# This is an 8001 bp read with 1 insert and 1 substitution
is "$(vg view -aj longread.gam | jq -r '.score')" "7989" "A long read can be correctly aligned"
is "$(vg view -aj longread.gam | jq -c '.path.mapping[].edit[] | select(.sequence)' | wc -l)" "2" "A long read hasd the correct edits found"
is "$(vg view -aj longread.gam | jq -c '. | select(.annotation["filter_3_cluster-coverage_cluster_passed_size_total"] <= 300)' | wc -l)" "1" "Long read minimizer set is correctly restricted"

rm -f longread.gam 1mb1kgp.dist 1mb1kgp.giraffe.gbz 1mb1kgp.min log.txt


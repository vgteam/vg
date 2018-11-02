#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 5

# Toy example of hand-made pileup (and hand inspected truth) to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg view -J -v call/tiny.json > tiny.vg

# With an empty pileup and loci mode we should assert the primery path.
true > empty.gam
vg augment tiny.vg empty.gam -a pileup -Z empty.aug.trans -S empty.aug.support > empty.aug.vg
vg call empty.aug.vg -z empty.aug.trans -s empty.aug.support -b tiny.vg --no-vcf > calls.loci
rm -f empty.gam

LOCUS_COUNT="$(vg view --locus-in -j calls.loci | wc -l)"
REF_LOCUS_COUNT="$(vg view --locus-in -j calls.loci | jq -c 'select(.genotype[0].allele == [0, 0])' | wc -l)"

is "${REF_LOCUS_COUNT}" "${LOCUS_COUNT}" "all loci on an empty pileup are called reference"

vg mod --sample-graph calls.loci empty.aug.vg > sample.vg

is "$(vg stats -l sample.vg)" "$(vg mod -k x tiny.vg | vg stats -l -)"  "called loci describe the primary path"

vg call empty.aug.vg -z empty.aug.trans -s empty.aug.support -b tiny.vg --no-vcf --call-nodes-by-coverage > calls.loci

LOCUS_COUNT="$(vg view --locus-in -j calls.loci | wc -l)"
EMPTY_LOCUS_COUNT="$(vg view --locus-in -j calls.loci | jq -c 'select(.genotype[0].allele == null)' | wc -l)"

is "${EMPTY_LOCUS_COUNT}" "${LOCUS_COUNT}" "all loci on an empty pileup in coverage-calling mode are called deleted"

rm -f tiny.vg empty.vgpu calls.loci sample.vg empty.aug.trans empty.aug.support empty.aug.vg

echo '{"node": [{"id": 1, "sequence": "CGTAGCGTGGTCGCATAAGTACAGTAGATCCTCCCCGCGCATCCTATTTATTAAGTTAAT"}]}' | vg view -Jv - > test.vg
vg index -x test.xg -g test.gcsa -k 16 test.vg
true >reads.txt
for REP in seq 1 5; do
    echo 'CGTAGCGTGGTCGCATAAGTACAGTANATCCTCCCCGCGCATCCTATTTATTAAGTTAAT' >>reads.txt
done
vg map -x test.xg -g test.gcsa --reads reads.txt > test.gam
cat test.gam | vg augment test.vg - -a pileup -P test.vgpu -Z test.trans -S test.support > test.aug.vg
vg call test.aug.vg -s test.support -z test.trans -b test.vg > /dev/null

N_COUNT=$(vg view -j test.aug.vg | grep "N" | wc -l)

is "${N_COUNT}" "0" "N bases are not augmented into the graph"

rm -rf reads.txt test.vg test.xg test.gcsa test.gcsa.lcp test.gam  test.aug.vg test.trans test.support empty.gam test.vgpu

vg construct -r inverting/miniFasta.fa -v inverting/miniFasta_VCFinversion.vcf.gz -S > miniFastaGraph.vg
vg index -x miniFastaGraph.xg -g miniFastaGraph.gcsa miniFastaGraph.vg
vg sim -x miniFastaGraph.xg -n 1000 -l 30 -a > miniFasta.gam
vg map -G miniFasta.gam -g miniFastaGraph.gcsa -x miniFastaGraph.xg > miniFastaGraph.gam
vg augment -a pileup -Z mappedminitest.trans -S mappedminitest.support --pileup mappedminitest.pileup miniFastaGraph.vg miniFastaGraph.gam > mappedminitest.aug.vg
vg call -z mappedminitest.trans -s mappedminitest.support -b miniFastaGraph.vg mappedminitest.aug.vg > calledminitest.vcf


L_COUNT=$(cat calledminitest.vcf | grep "#" -v | wc -l)
is "${L_COUNT}" "1" "Called microinversion"
 
rm -f miniFastaGraph.vg miniFasta.gam miniFastaGraph.gam mappedminitest.aug.vg calledminitest.vcf mappedminitest.trans mappedminitest.support mappedminitest.pileup miniFastaGraph.xg miniFastaGraph.gcsa





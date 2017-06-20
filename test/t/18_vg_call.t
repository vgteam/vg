#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 5

# Toy example of hand-made pileup (and hand inspected truth) to make sure some
# obvious (and only obvious) SNPs are detected by vg call
vg view -J -v call/tiny.json > tiny.vg

vg view -J -l call/pileup.json -L > tiny.vgpu
vg call tiny.vg tiny.vgpu -A calls_l.vg   > /dev/null 2> /dev/null
is $? "0" "vg call doesn't crash"

rm -f calls_l.json calls_l.vg tiny.vgpu

# With an empty pileup and loci mode we should assert the primery path.
true > empty.gam
vg pileup tiny.vg empty.gam > empty.vgpu
vg call tiny.vg empty.vgpu --aug-graph augmented.vg --no-vcf > calls.loci

LOCUS_COUNT="$(vg view --locus-in -j calls.loci | wc -l)"
REF_LOCUS_COUNT="$(vg view --locus-in -j calls.loci | jq -c 'select(.genotype[0].allele == [0, 0])' | wc -l)"

is "${REF_LOCUS_COUNT}" "${LOCUS_COUNT}" "all loci on an empty pileup are called reference"

vg mod --sample-graph calls.loci augmented.vg > sample.vg

is "$(vg stats -l sample.vg)" "$(vg mod -k x tiny.vg | vg stats -l -)"  "called loci describe the primary path"

vg call tiny.vg empty.vgpu --aug-graph augmented.vg --no-vcf --call-nodes-by-coverage > calls.loci

LOCUS_COUNT="$(vg view --locus-in -j calls.loci | wc -l)"
EMPTY_LOCUS_COUNT="$(vg view --locus-in -j calls.loci | jq -c 'select(.genotype[0].allele == null)' | wc -l)"

is "${EMPTY_LOCUS_COUNT}" "${LOCUS_COUNT}" "all loci on an empty pileup in coverage-calling mode are called deleted"

rm -f tiny.vg empty.vgpu augmented.vg calls.loci sample.vg

echo '{"node": [{"id": 1, "sequence": "CGTAGCGTGGTCGCATAAGTACAGTAGATCCTCCCCGCGCATCCTATTTATTAAGTTAAT"}]}' | vg view -Jv - > test.vg
vg index -x test.xg -g test.gcsa -k 16 test.vg
true >reads.txt
for REP in seq 1 5; do
    echo 'CGTAGCGTGGTCGCATAAGTACAGTANATCCTCCCCGCGCATCCTATTTATTAAGTTAAT' >>reads.txt
done
vg map -x test.xg -g test.gcsa --reads reads.txt > test.gam
vg pileup test.vg test.gam > test.vgpu
vg call test.vg test.vgpu >/dev/null --aug-graph aug.vg

N_COUNT=$(vg view -j aug.vg | grep "N" | wc -l)

is "${N_COUNT}" "0" "N bases are not augmented into the graph"

rm -r reads.txt test.vg test.xg test.gcsa test.gcsa.lcp test.gam test.vgpu aug.vg





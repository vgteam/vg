#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 4

vg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg
vg snarls --include-trivial x.vg > x.snarls
vg index -s x.snarls -j x.dist x.vg
vg minimizer -k 29 -w 11 -g x.gbwt -i x.min x.xg

vg gaffe -x x.xg -H x.gbwt -m x.min -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with all indexes specified without crashing"

rm -f x.vg x.xg x.gbwt x.snarls x.min x.dist x.gg

cp small/x.fa .
cp small/x.vcf.gz .
cp small/x.vcf.gz.tbi .

vg gaffe x.fa x.vcf.gz -f reads/small.middle.ref.fq > mapped2.gam
is "${?}" "0" "a read can be mapped with just FASTA and VCF without crashing"

# These files can differ as serialized and still represent the same data, dur to protobuf field order not being specified.
# Tripping through JSON will sort all the keys.
vg view -aj mapped1.gam >mapped1.json
vg view -aj mapped2.gam >mapped2.json

# Make sure at least one file converted successfully
SIZE="$(wc -c mapped2.json | cut -f1 -d' ')"
echo "Size: ${SIZE}"
EMPTY=0
if [ "${SIZE}" == "0" ] ; then
    EMPTY=1
fi
is "${EMPTY}" "0" "mapping with just a FASTA and a VCF produced JSON-able alignments"

diff mapped1.json mapped2.json
is "${?}" "0" "mapping to manually-generated indexes and automatically-generated indexes is the same"

rm -f x.vg x.gbwt x.gg x.snarls x.min x.dist x.gg x.fa x.fa.fai x.vcf.gz x.vcf.gz.tbi mapped1.gam mapped1.json mapped2.gam mapped2.json



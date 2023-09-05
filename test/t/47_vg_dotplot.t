#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 > hla.vg
vg index -x hla.xg hla.vg
vg dotplot -x hla.xg >/dev/null

is "${?}" "0" "vg dotplot runs successfully"

rm -f hla.vg hla.xg



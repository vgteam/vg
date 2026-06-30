#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg index -x hla.xg msgas/hla.vg
vg dotplot -x hla.xg >/dev/null

is "${?}" "0" "vg dotplot runs successfully"

rm -f hla.xg



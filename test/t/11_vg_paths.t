#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 1

is "$(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg paths --list -v -)" "x" "path listing works"

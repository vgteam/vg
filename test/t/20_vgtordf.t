#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2
vgtordf.sh <(vg view -V  graphs/199754000\:199755000.vg  -j) | jq . 
is $? 0 "Basic syntax check passes"

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -t -r 'http://example.org' - | rapper -i turtle -I 'http://example.org/' -o turtle - 2>&1 | grep 'returned 65 triples' | wc -l) 1 "vg turtle output works"

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 2

vg sort -r ref -g -i sort/test_sort.gfa > sorted.vg

is $(vg view -g sorted.vg | grep "^S" | wc -l) 13 "sorted graph contains expected number of nodes"
is $(vg view -d sorted.vg | grep label | tr '\n' ' ' | grep -E '1:C(.*)6:T(.*)7:C(.*)8:A(.*)2:G(.*)10:G(.*)9:C(.*)3:C(.*)4:G(.*)5:C(.*)11:A(.*)12:T(.*)13:C' | wc -l) 1 "nodes appear in the expected order"

rm -f sorted.vg

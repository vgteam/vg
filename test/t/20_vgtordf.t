#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1
vgtordf.sh <(vg view -V  graphs/199754000\:199755000.vg  -j) | jq . 
is $? 0 "Basic syntax check passes"

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

vg view -J -v snarls/snarls.json > snarls.vg
is $(vg snarls snarls.vg -r st.pb | vg view -R - | wc -l) 3 "vg snarls made right number of protobuf Snarls"
is $(vg view -E st.pb | wc -l) 6 "vg snarls made right number of protobuf SnarlTraversals"

rm -f snarls.vg st.pb 
 


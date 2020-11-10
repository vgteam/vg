#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 5

# Make a big graph that is probably not in sorted order
# Chopping it accomplishes that and makes it have several chunks
vg construct -r minigiab/q.fa -v minigiab/NA12878.chr22.tiny.giab.vcf.gz -m 64 | vg mod -X 1 - >giab.vg

vg view -j giab.vg | jq -cr '.node[].id' > ids.txt
sort -n ids.txt > ids.sorted.txt

diff ids.txt ids.sorted.txt >/dev/null
is "${?}" "1" "ids in our test graph are not initially in sorted order"

# Sort by ID and index
vg sort -a id -I giab.sorted.vg.vgi giab.vg > giab.sorted.vg

is "${?}" "0" "a vg graph can be sorted and indexed by ID without crashing"

vg sort -a topo giab.vg >/dev/null
is "${?}" "0" "a vg graph can be sorted topologically"

vg sort -a eades -r q giab.vg >/dev/null
is "${?}" "0" "a vg graph can be sorted with Eades algorithm without crashing"

vg sort -a max-flow -r q giab.vg >/dev/null
is "${?}" "0" "a vg graph can be sorted with the max-flow algorithm without crashing"

rm -f giab.vg giab.sorted.vg giab.sorted.vg.vgi ids.txt ids.sorted.txt


 


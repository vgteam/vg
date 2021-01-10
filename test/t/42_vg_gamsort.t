#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 2

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg  x.vg
vg sim -n 1000 -l 100 -e 0.01 -i 0.005 -x x.xg -a >x.gam

vg gamsort x.gam >x.sorted.gam

vg view -aj x.sorted.gam | jq -r '.path.mapping | ([.[] | .position.node_id | tonumber] | min)' >min_ids.gamsorted.txt
vg view -aj x.gam | jq -r '.path.mapping | ([.[] | .position.node_id | tonumber] | min)' | sort -n >min_ids.sorted.txt

is "$(md5sum <min_ids.gamsorted.txt)" "$(md5sum <min_ids.sorted.txt)" "Sorting a GAM orders the alignments by min node ID"

vg gamsort x.gam -i x.sorted.gam.gai >x.sorted.gam
is "$?" "0" "sorted GAMs can be indexed during the sort"


rm -f x.vg x.xg x.gam x.sorted.gam x.sorted.2.gam min_ids.gamsorted.txt min_ids.sorted.txt x.sorted.gam.gai x.sorted.2.gam.gai

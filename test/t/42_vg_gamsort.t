#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg


plan tests 15

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg  x.vg
vg sim -n 1000 -l 100 -e 0.01 -i 0.005 -x x.xg -a >x.gam


# GAM
vg gamsort x.gam >x.sorted.gam

vg view -aj x.sorted.gam | jq -r '.path.mapping | ([.[] | .position.node_id | tonumber] | min)' >min_ids.gamsorted.txt
vg view -aj x.gam | jq -r '.path.mapping | ([.[] | .position.node_id | tonumber] | min)' | sort -n >min_ids.sorted.txt

is "$(md5sum <min_ids.gamsorted.txt)" "$(md5sum <min_ids.sorted.txt)" "Sorting a GAM orders the alignments by min node ID"

vg gamsort x.gam -i x.sorted.gam.gai >x.sorted.gam
is "$?" "0" "sorted GAMs can be indexed during the sort"

vg gamsort --shuffle x.sorted.gam >x.shuffled.gam
is "$?" "0" "GAMs can be shuffled"
is "$(vg stats -a x.shuffled.gam)" "$(vg stats -a x.sorted.gam)" "Shuffling preserves read data"

rm -f x.sorted.gam x.sorted.2.gam x.shuffled.gam min_ids.gamsorted.txt min_ids.sorted.txt x.sorted.gam.gai x.sorted.2.gam.gai


# GAF. Correctness is tested in unit tests, so we just check that the commands work.
vg convert -G x.gam x.xg > x.gaf
sort x.gaf > x.gaf.lexicographic
rm -f x.gbwt

vg gamsort -G x.gaf > x.sorted.gaf
is "$?" "0" "GAFs can be sorted"
is "$(ls x.gbwt 2> /dev/null)" "" "GBWT index is not created by default"
sort x.sorted.gaf > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF preserves read data"

vg gamsort -G x.gaf -g x.gbwt > x.sorted.gaf
is "$?" "0" "GAFs can be sorted while creating a GBWT index"
is "$(ls x.gbwt 2> /dev/null)" "x.gbwt" "GBWT index can be created"
sort x.sorted.gaf > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF preserves read data with GBWT index"

vg gamsort -G x.gaf -g x.gbwt --bidirectional > x.sorted.gaf
is "$?" "0" "GAFs can be sorted while creating a bidirectional GBWT index"
is "$(ls x.gbwt 2> /dev/null)" "x.gbwt" "Bidirectional GBWT index can be created"
sort x.sorted.gaf > x.sorted.gaf.lexicographic
cmp x.gaf.lexicographic x.sorted.gaf.lexicographic
is "$?" "0" "Sorting a GAF preserves read data with bidirectional GBWT index"

vg gamsort -G --shuffle x.gaf > x.shuffled.gaf
is "$?" "0" "GAFs can be shuffled"
sort x.shuffled.gaf > x.shuffled.gaf.lexicographic
cmp x.gaf.lexicographic x.shuffled.gaf.lexicographic
is "$?" "0" "Shuffling a GAF preserves read data"

rm -f x.gaf x.gaf.lexicographic x.sorted.gaf x.sorted.gaf.lexicographic x.shuffled.gaf x.shuffled.gaf.lexicographic
rm -f x.gbwt


# Cleanup
rm -f x.vg x.xg x.gam

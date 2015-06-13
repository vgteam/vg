#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 6

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^P | wc -l) \
    $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -k x - | vg view - | grep ^S | wc -l) \
    "vg mod yields a graph with only a particular path"

is $(vg mod -o graphs/orphans.vg | vg view - | wc -l) 8 "orphan edge removal works"

vg construct -r tiny/tiny.fa >t.vg
vg index -s -k 11 t.vg

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 1 "path inclusion does not modify the graph when alignment is a perfect match"

is $(vg map -s CAAATAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 5 "path inclusion with a complex variant introduces the right number of nodes"

is $(vg map -s CAAAAAGGCTTGGAAAGGGTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | wc -l) 24 "path inclusion works for deletions"

is $(vg map -s CAAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 1 "soft clips in alignments don't affect the graph when we introduce the alignment paths into the graph"

vg map -r <(vg sim -l 50 -n 10 t.vg -e 0.05 -i 0.005 -s 9669) t.vg | vg mod -i - -c t.vg >y.vg
vg index -s -k 11 y.vg

score_sum_before=$(vg map -r <(vg sim -l 50 -n 10 t.vg -e 0.05 -i 0.005 -s 9669) t.vg | vg view -a - | jq .score | awk '{ sum += $1 } END { print sum }')
score_sum_after=$(vg map -r <(vg sim -l 50 -n 10 t.vg -e 0.05 -i 0.005 -s 9669) y.vg | vg view -a - | jq .score | awk '{ sum += $1 } END { print sum }')

#echo $score_sum_before
#echo $score_sum_after

#is $(echo $score_sum_after / $score_sum_before | bc) 1 "modifying the graph to include alignments improves mapping of original reads"

rm t.vg y.vg
rm -rf t.vg.index y.vg.index

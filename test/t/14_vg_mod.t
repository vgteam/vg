#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 9

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

is $(vg map -s CAAATAAGGCTTGGAAATTTTCTGCAGTTCTATTATATTCCAACTCTCTG t.vg | vg mod -i - t.vg | vg view - | grep ^S | wc -l) 4 "SNPs can be included in the graph"

rm t.vg
rm -rf t.vg.index

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 8 -e 4 - | vg kmers -g -k 8 - | sort | md5sum | awk '{ print $1 }') 7caff8f7497ad9cf2762e2d95d54a241 "graph complexity reduction works as expected"

is $( vg construct -r small/x.fa -v small/x.vcf.gz | vg mod -pl 8 -e 4 -t 16 - | vg mod -S -l 200 - | vg view - | grep ^S | wc -l) 152 "short subgraph pruning works"

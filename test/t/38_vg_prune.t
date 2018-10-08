#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 15


# Build a graph with one path and two threads
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg

# Basic pruning: 5 components, 31 nodes, 31 edges
vg prune -e 1 x.vg > y.vg
is $(vg stats -s y.vg | wc -l) 5 "pruning produces the correct number of components"
is $(vg stats -N y.vg) 31 "pruning leaves the correct number of nodes"
is $(vg stats -E y.vg) 31 "pruning leaves the correct number of edges"
rm -f y.vg

# Restore paths: 1 component, 44 nodes, 48 edges
vg prune -r -e 1 x.vg > y.vg
is $(vg stats -s y.vg | wc -l) 1 "pruning with path restoring produces the correct number of components"
is $(vg stats -N y.vg) 44 "pruning with path restoring leaves the correct number of nodes"
is $(vg stats -E y.vg) 48 "pruning with path restoring leaves the correct number of edges"
rm -f y.vg

# Unfold paths and threads: 1 component, 60 nodes, 72 edges
vg prune -u -g x.gbwt -e 1 x.vg > y.vg
is $(vg stats -s y.vg | wc -l) 1 "pruning with path/thread unfolding produces the correct number of components"
is $(vg stats -N y.vg) 60 "pruning with path/thread unfolding produces the correct number of nodes"
is $(vg stats -E y.vg) 72 "pruning with path/thread unfolding produces the correct number of edges"
rm -f y.vg

# Unfold only paths: 1 component, 44 nodes, 48 edges
vg prune -u -e 1 x.vg > y.vg
is $(vg stats -s y.vg | wc -l) 1 "pruning with path unfolding produces the correct number of components"
is $(vg stats -N y.vg) 44 "pruning with path unfolding produces the correct number of nodes"
is $(vg stats -E y.vg) 48 "pruning with path unfolding produces the correct number of edges"
rm -f y.vg

rm -f x.vg x.gbwt


# Build vg graphs for two chromosomes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j -m xy.mapping x.vg y.vg
vg index -G x.gbwt -v small/xy2.vcf.gz x.vg
vg index -G y.gbwt -v small/xy2.vcf.gz y.vg

# Check the range in the node mapping at each point
is $(head -c 16 xy.mapping | hexdump -v -e '4/4 "%u_"') 99_0_99_0_ "empty mapping starts from the right node id"
vg prune -u -g x.gbwt -e 1 -a -m xy.mapping x.vg > /dev/null
is $(head -c 16 xy.mapping | hexdump -v -e '4/4 "%u_"') 99_0_128_0_ "the first unfolded graph adds the correct number of nodes to the mapping"
vg prune -u -g y.gbwt -e 1 -a -m xy.mapping y.vg > /dev/null
is $(head -c 16 xy.mapping | hexdump -v -e '4/4 "%u_"') 99_0_156_0_ "the second unfolded graph adds the correct number of nodes to the mapping"

rm -f x.vg y.vg x.gbwt y.gbwt xy.mapping

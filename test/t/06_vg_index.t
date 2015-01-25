#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 9

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -s x.vg
is $? 0 "indexing nodes and edges of graph"

vg index -k 11 x.vg
is $? 0 "indexing 11mers"

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
is $? 0 "construction of 1mb graph succeeds"

vg index -s z.vg
is $? 0 "indexing nodes, and edges of 1mb graph"

is $(tail -1 z.vg.index/LOG | grep ERROR | wc -l) 0 "indexing does not close on error"

# clean up
rm -rf z.vg.index z.vg

num_records=$(vg index -D x.vg | wc -l)
is $? 0 "dumping graph index"
is $num_records 7933 "correct number of records in graph index"

rm -rf x.vg.index
rm -f x.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg ids -j x.vg y.vg z.vg
vg index -s -d q.vg.index x.vg y.vg z.vg
vg index -s x.vg

single=$(vg index -D -d x.vg.index | wc -l)
triple=$(vg index -D -d q.vg.index | wc -l)

is $triple $(echo "$single * 3" | bc) "storage of multiple graphs in an index succeeds"

rm x.vg y.vg z.vg
rm -rf x.vg.index q.vg.index

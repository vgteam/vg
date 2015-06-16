#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 14

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -s x.vg
is $? 0 "indexing nodes and edges of graph"

vg index -k 11 x.vg
is $? 0 "indexing 11mers"

vg index -C x.vg
is $? 0 "index compaction"

#vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
#is $? 0 "construction of 1mb graph succeeds"

#vg index -s z.vg
#is $? 0 "indexing nodes, and edges of 1mb graph"

#is $(tail -1 z.vg.index/LOG | grep ERROR | wc -l) 0 "indexing does not close on error"

# clean up
#rm -rf z.vg.index z.vg

num_records=$(vg index -D x.vg | wc -l)
is $? 0 "dumping graph index"
is $num_records 3207 "correct number of records in graph index"

vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg index -a - -d x.vg.aln
is $(vg index -D -d x.vg.aln | wc -l) 100 "index can store alignments"
is $(vg index -A -d x.vg.aln | vg view -a - | wc -l) 100 "index can dump alignments"

vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg index -m - -d x.vg.map
is $(vg index -D -d x.vg.map | wc -l) $(vg map -r <(vg sim -s 1337 -n 100 x.vg) x.vg | vg view -a - | jq -c '.path.mapping[]' | sort | uniq | wc -l) "index stores all unique mappings"

rm -rf x.vg.index x.vg.map x.vg.aln
rm -f x.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg concat x.vg y.vg z.vg >q.vg

vg index -s q.vg
vg index -s x.vg

single=$(vg index -D -d x.vg.index | wc -l)
triple=$(vg index -D -d q.vg.index | wc -l)

# subtract two for metadata lines about paths that aren't duplicated in the merged
is $triple $(echo "$single * 3 - 2" | bc) "storage of multiple graphs in an index succeeds"

rm x.vg y.vg z.vg q.vg
rm -rf x.vg.index q.vg.index

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -s x.vg
is $(vg index -D x.vg | grep +g | grep +p | wc -l) $(vg view x.vg | grep ^P | wc -l) "correct number of elements in path index"
is $(vg index -D x.vg | grep +path_id | wc -l) 1 "path id recorded"
is $(vg index -D x.vg | grep +path_name | wc -l) 1 "path name recorded"
rm -rf x.vg.index x.vg

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -v small/x.vcf.gz -r small/x.fa | vg view - | sed s/x/y/ | vg view -v - >y.vg
vg ids -j x.vg y.vg

vg index -s -d q.idx x.vg y.vg
is $(vg index -L -d q.idx | tail -1 | awk '{ print $3 }' ) 420 "end of the second path found correctly"

rm -rf q.idx x.vg y.vg


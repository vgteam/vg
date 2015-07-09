#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 19

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg stats -z - | grep nodes | cut -f 2) 210 "construction produces the right number of nodes"

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg stats -z - | grep edges | cut -f 2) 291 "construction produces the right number of edges"

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

nodes=$(vg stats -z z.vg | head -1 | cut -f 2)
is $nodes 84553 "the 1mb graph has the expected number of nodes"

edges=$(vg stats -z z.vg | tail -1 | cut -f 2)
is $edges 115357 "the 1mb graph has the expected number of edges"

rm -f z.vg

vg construct -r complex/c.fa -v complex/c.vcf.gz >c.vg
is $? 0 "construction of a very complex region succeeds"

nodes=$(vg stats -z c.vg | head -1 | cut -f 2)
is $nodes 71 "the complex graph has the expected number of nodes"

edges=$(vg stats -z c.vg | tail -1 | cut -f 2)
is $edges 116 "the complex graph has the expected number of edges"

rm -f c.vg

order_a=$(vg construct -r order/n.fa -v order/x.vcf.gz | md5sum | cut -f 1 -d\ )
order_b=$(vg construct -r order/n.fa -v order/y.vcf.gz | md5sum | cut -f 1 -d\ )

is $order_a $order_b "the ordering of variants at the same position has no effect on the resulting graph"

vg construct -r order/n.fa -v order/z.vcf.gz -R n:47-73 >/dev/null
is $? 0 "construction does not fail when the first position in the VCF is repeated and has an indel"

x1=$(for i in $(seq 100); do size=$(shuf -i 1-100 -n 1); threads=1; vg construct -r small/x.fa -v small/x.vcf.gz -z $size -t $threads | vg view -g - | sort -n -k 2 | md5sum; done | sort | uniq | wc -l)

is $x1 1 "the size of the regions used in construction has no effect on the graph"

x2=$(for i in $(seq 100); do
    size=10;
    threads=$(shuf -i 1-100 -n 1);
    vg construct -r small/x.fa -v small/x.vcf.gz -z $size -t $threads | vg view -g - | sort -n -k 2 | md5sum;
    done | sort | uniq | wc -l)

is $x2 1 "the number of threads used in construction has no effect on the graph"

x3=$(for i in $(seq 100); do
    size=$(shuf -i 1-100 -n 1);
    threads=$(shuf -i 1-100 -n 1);
    vg construct -r small/x.fa -v small/x.vcf.gz -z $size -t $threads | vg view -g - | sort -n -k 2 | md5sum;
    done | sort | uniq | wc -l)

is $x3 1 "the number of threads and regions used in construction has no effect on the graph"

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz -R z:10-20 >/dev/null
is $? 0 "construction of a graph with two head nodes succeeds"

# in case there were failures in topological sort
rm -f fail.vg

# check that we produce a full graph

# build deps
cd ../vcflib && make vcf2tsv && cd -
cd ../fastahack && make && cd -

refbp=$(../fastahack/fastahack -r x small/x.fa | tr '\n' ' ' | sed 's/ //' | wc -c)
variantbp=$(zcat small/x.vcf.gz | ../vcflib/bin/vcf2tsv \
    | cut -f 5,4 | tail -n+2 \
    | awk '{ x=length($2)-length($1); if (x > 0) { print x; } else if (x == 0) { print length($2); } }' \
        | awk '{ sum += $1 } END { print sum }')
rm ../vcflib/bin/vcf2tsv ../fastahack/fastahack

graphbp=$(vg construct -r small/x.fa -v small/x.vcf.gz | vg stats -l - | cut -f 2)

is $graphbp $(echo "$refbp + $variantbp" | bc) "the graph contains all the sequence in the reference and VCF"

is $(for i in $(seq 100); do vg construct -r small/x.fa -v small/x.vcf.gz -m $i | vg stats -l - | cut -f 2; done | sort | uniq | wc -l) 1 "varying the max node size does not affect graph length"

max_node_size=$(vg construct -r small/x.fa -v small/x.vcf.gz -m 12 | vg view -g - | grep ^S | cut -f 3 | awk '{ print length($1) }' | sort -n | tail -1)

is $max_node_size 12 "nodes are correctly capped in size"

is $(vg construct -R z:10000-20000 -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz | vg view - | awk '{ print length($3); }' | sort -n | tail -1) 241 "-R --region flag is respected" 

is $(vg construct -r small/x.fa -m 50 | vg view - | wc -l) 63 "vg construct does not require a VCF and respects node size limit"

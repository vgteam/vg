#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 30

is $(vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | vg stats -z - | grep nodes | cut -f 2) 210 "construction produces the right number of nodes"

is $(vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | vg stats -z - | grep edges | cut -f 2) 291 "construction produces the right number of edges"

is $(vg construct -r small/x.fa --rename chrX=x -R chrX:1-2 | vg stats -l - | cut -f 2 ) 2 "construction obeys rename and region options"

vg construct -m 1000 -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz >z.vg
is $? 0 "construction of a 1 megabase graph from the 1000 Genomes succeeds"

nodes=$(vg stats -z z.vg | head -1 | cut -f 2)
is $nodes 84559 "the 1mb graph has the expected number of nodes"

edges=$(vg stats -z z.vg | tail -1 | cut -f 2)
is $edges 115375 "the 1mb graph has the expected number of edges"

# Our 1000 Genomes subset file has changed variant start positions but not stored end positions.
vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz -S >/dev/null
is "${?}" "0" "vg construct can skip self-inconsistent structural variants"

rm -f z.vg

is $(vg construct -r 1mb1kgp/z.fa | vg view -j - | jq -c '.node[] | select((.sequence | length) >= 1024)' | wc -l) 0 "node size is manageable by default"

vg construct -m 1000 -r complex/c.fa -v complex/c.vcf.gz >c.vg
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

refbp=$(../deps/fastahack/fastahack -r x small/x.fa | tr -d '\n' | wc -c)
variantbp=$(zcat < small/x.vcf.gz | ../deps/vcflib/bin/vcf2tsv \
    | cut -f 5,4 | tail -n+2 \
    | awk '{ x=length($2)-length($1); if (x > 0) { print x; } else if (x == 0) { print length($2); } }' \
        | awk '{ sum += $1 } END { print sum }')

graphbp=$(vg construct -r small/x.fa -v small/x.vcf.gz | vg stats -l - | cut -f 2)

is $graphbp $(echo "$refbp + $variantbp" | bc) "the graph contains all the sequence in the reference and VCF"

is $(for i in $(seq 100); do vg construct -r small/x.fa -v small/x.vcf.gz -m $i | vg stats -l - | cut -f 2; done | sort | uniq | wc -l) 1 "varying the max node size does not affect graph length"

max_node_size=$(vg construct -r small/x.fa -v small/x.vcf.gz -m 12 | vg view -g - | grep ^S | cut -f 3 | awk '{ print length($1) }' | sort -n | tail -1)

is $max_node_size 12 "nodes are correctly capped in size"

## Check the length of the longest node
is $(vg construct -m 1000 -R z:10000-20000 -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz | vg view - | grep "^S" | awk '{ print length($3); }' | sort -n | tail -1) 241 "-R --region flag is respected"

vg construct -r small/x.fa >/dev/null
is $? 0 "vg construct does not require a vcf"

(( short_enough=$(vg construct -r small/x.fa -m 50 | vg view - | grep "^S" | cut -f3 | wc -l) < 51 ? 1 : 0 ))
is $short_enough 1 "vg construct respects node size limit"

is $(vg construct -CR 'gi|568815592:29791752-29792749' -r GRCh38_alts/FASTA/HLA/V-352962.fa | vg view - | grep TCTAGAAGAGTCCACGGGGACAGGTAAG | wc -l) 1 "--region can be interpreted to be a reference sequence (and not parsed as a region spec)"

is "$(vg construct -r sv/x.fa -v sv/x.inv.vcf -S | vg view - | sort | md5sum | cut -f 1 -d\ )" "$(cat sv/x.inv.gfa | sort | md5sum | cut -f1 -d\ )" "vg constructs the correct graph for inversions"

is "$(vg construct -r tiny/tiny.fa -v tiny/dots.vcf | vg mod -u - | vg stats -z - | grep nodes | cut -f2)" "1" "vg construct skips variants with . ALTs"

vg construct -r small/x.fa -v small/x.v4.3.vcf.gz >x.vg 2>log.txt
is $? 0  "VCF with * alleles is accepted"
grep "skipping variant" log.txt >/dev/null
is $? 0  "VCF record with * allele is rejected"
rm -f x.vg log.txt

vg construct -r tiny/ambiguous.fa >tiny.vg
is $? 0 "Reference with ambiguity codes has them coerced to Ns"
is "$(vg view -j tiny.vg | jq -r '.node[].sequence' | tr -d 'ACGT\n' | wc -c)" "10" "Expected number of Ns are created"
rm -f tiny.vg


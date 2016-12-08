#!/bin/bash

if [ ! $# -eq 6 ];
then
    echo "usage: " $0 " [prefix] '[input spec]' [linear_xg] [linear_gcsa] [pan_xg] [pan_gcsa]"
    exit
fi

prefix=$1
input=$2
linear_xg=$3
linear_gcsa=$4
pan_xg=$5
pan_gcsa=$6

# take the pan and ref indexes
# map a set of reads against both

mkdir -p $prefix

echo mapping $prefix against linear graph
vg map -x $linear_xg -g $linear_gcsa $input -u 1 -L 10 -a >$prefix/linear.gam

echo mapping $prefix against pan graph
vg map -x $pan_xg -g $pan_gcsa $input -u 1 -L 10 -a >$prefix/pan.gam

# join the results (identity, score, mapping quality) on read name
echo comparing results
( echo group name ref.identity ref.score ref.mapq pan.identity pan.score pan.mapq
  join -j 1 \
       <(vg view -a $prefix/linear.gam | jq '[.name, .identity, .score, .mapping_quality] | @sh' | sed 's/"//g' | sed "s/'//g" | sed 's/ null/ 0/g' | sort) \
       <(vg view -a $prefix/pan.gam | jq '[.name, .identity, .score, .mapping_quality] | @sh' | sed 's/"//g' | sed "s/'//g" | sed 's/ null/ 0/g' | sort) \
      | while read line; do echo $prefix $line; done ) | tr ' ' '\t' | gzip >$prefix/vs.tsv.gz

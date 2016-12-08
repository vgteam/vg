#!/bin/bash

if [ ! $# -eq 7 ];
then
    echo "usage: " $0 " [prefix] [fq1] [fq2] [linear_xg] [linear_gcsa] [pan_xg] [pan_gcsa]"
    exit
fi

prefix=$1
fq1=$2
fq2=$3
linear_xg=$4
linear_gcsa=$5
pan_xg=$6
pan_gcsa=$7

# take the pan and ref indexes
# map a set of reads against both

mkdir -p $prefix

echo mapping against linear graph
vg map -x $linear_xg -g $linear_gcsa -f $fq1 -f $fq2 -u 1 -L 10 >$prefix/linear.gam

echo mapping against pan graph
vg map -x $pan_xg -g $pan_gcsa -f $fq1 -f $fq2 -u 1 -L 10 >$prefix/pan.gam

# join the results (identity, score, mapping quality) on read name
echo comparing results
( echo 'name\tref.identity\tref.score\tref.mapq\tpan.identity\tpan.score\tpan.mapq'
  join -j 1 \
       <(vg view -a $prefix/linear.gam | jq '[.name, .identity, .score, .mapping_quality] | @sh' | sed 's/"//g' | sed s/null/0/ | sort) \
       <(vg view -a $prefix/pan.gam | jq '[.name, .identity, .score, .mapping_quality] | @sh' | sed 's/"//g' | sed s/null/0/ | sort)
) | gzip >$prefix/vs.tsv.gz

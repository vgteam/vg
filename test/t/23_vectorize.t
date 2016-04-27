#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 0

#vg construct -r ../tiny/tiny.fa -v ../tiny/tiny.vcf.gz > tiny.vg
#vg index -x tiny.xg -g tiny.gcsa -k 16 tiny.vg
#vg sim -l 10 -s 3 -x tiny.xg tiny.vg > tiny.reads

## Can we handle streaming alignments?
#is $(vg map -x tiny.xg -g tiny.gcsa -r tiny.reads | vg vectorize -x tiny.xg - | wc -l) 10 "Streaming produces the correct number of vectors."

## Can we read from a file?
#is $(vg vectorize -x tiny.xg tiny.reads | wc -l) 10 "Reading from file produces the expected number of vectors."


## Test that vectorize correctly vectorizes a small alignment in a_hot format
#is $(vg vectorize -a -x tiny.xg tiny.gam | md5sum) hh345jj "Alignments can be vectorized to a-hot format."

## Test if we get the expected one-hot format
#is $(vg vectorize -x tiny.xg tiny.gam | md5sum) hh345jj "Alignments can be vectorized to one-hot format."

## Test if the alignment can be output in identity_hot
#is $(vg vectorize -x tiny.xg tiny.gam | md5sum) hh345jj "Alignments can be vectorized to identity-hot format."

## Test if names can be relabeled using -l
#is $(vg vectorize -l test -w -x tiny.xg tiny.gam | head -n 1 | cut -f 1) "test" "Vectorized alignments can be renamed using -l."


## Test if vectorize correctly produces formatted output (tab separated with sequence name)

## Check if vectorize can produce expected wabbit output.
#is $(vg vectorize -w -x tiny.xg tiny.gam | head -n 1 | md5sum) hh345jj "Vectorize -w produces the expected vowpal-wabbit compatible output."


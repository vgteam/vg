#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 3

vg index -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -k 16 graphs/refonly-lrc_kir.vg

vg mpmap -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -I 10 -D 50 -S | vg view -aj -  > temp_paired_alignment.json
vg mpmap -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -I 100000 -D 5 -S | vg view -aj -  > temp_distant_alignment.json
vg mpmap -x graphs/refonly-lrc_kir.vg.xg -g graphs/refonly-lrc_kir.vg.gcsa -f reads/grch38_lrc_kir_paired.fq -i -S | vg view -aj - > temp_independent_alignment.json
paired_score=$(jq -r ".score" < temp_paired_alignment.json | awk '{ sum+=$1} END {print sum}')
independent_score=$(jq -r ".score" < temp_independent_alignment.json | awk '{ sum+=$1} END {print sum}')
is $(printf "%s\t%s\n" $paired_score $independent_score | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent have lower score than unrestricted alignments"

paired_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_paired_alignment.json | sort | rs -T | awk '{print ($2 - $1)}')
distant_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_distant_alignment.json | sort | rs -T | awk '{print ($2 - $1)}')
independent_range=$(jq -r ".path.mapping[0].position.node_id" <  temp_independent_alignment.json | sort | rs -T | awk '{print ($2 - $1)}')
is $(printf "%s\t%s\n" $paired_range $independent_range | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be consistent are closer together in node id space than unrestricted alignments"
is $(printf "%s\t%s\n" $paired_range $distant_range | awk '{if ($1 < $2) print 1; else print 0}') 1 "paired read alignments forced to be near each other are closer together in node id space than those forced to be far apart"

rm temp_paired_alignment.json temp_independent_alignment.json

rm -f graphs/refonly-lrc_kir.vg.xg graphs/refonly-lrc_kir.vg.gcsa graphs/refonly-lrc_kir.vg.gcsa.lcp

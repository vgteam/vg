#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 5

vg autoindex -p auto -w map -r tiny/tiny.fa -v tiny/tiny.vcf.gz -t 1
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with basic input"
is $(ls auto* | wc -l) 3 "autoindexing makes 3 outputs for vg map" 
is $(ls auto.xg | wc -l) 1 "autoindexing makes an XG for vg map"
is $(ls auto.gcsa* | wc -l) 2 "autoindexing makes a GCSA2/LCP pair for vg map"

rm auto.xg auto.gcsa auto.gcsa.lcp

vg autoindex -p auto -w map -r small/x.fa -v small/x.vcf.gz -r small/y.fa -v small/y.vcf.gz -t 1
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with chunked input"

rm auto.xg auto.gcsa auto.gcsa.lcp
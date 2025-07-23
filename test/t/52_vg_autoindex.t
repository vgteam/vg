#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 56

rm auto.*

vg autoindex -p auto -w map -r tiny/tiny.fa -v tiny/tiny.vcf.gz --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with basic input"
is $(ls auto.* | wc -l) 3 "autoindexing makes 3 outputs for vg map" 
is $(ls auto.xg | wc -l) 1 "autoindexing makes an XG for vg map"
is $(ls auto.gcsa* | wc -l) 2 "autoindexing makes a GCSA2/LCP pair for vg map"
vg sim -x auto.xg -n 20 -a -l 10 | vg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg map"

old_xg_creation_time=`ls -l --time-style=full-is auto.xg | tr -s ' ' | cut -d ' ' -f7`
vg autoindex -p auto -w map -r tiny/tiny.fa -v tiny/tiny.vcf.gz --force-unphased --no-guessing
cur_xg_creation_time=`ls -l --time-style=full-is auto.xg | tr -s ' ' | cut -d ' ' -f7`
same_time="no"
if [ "$old_xg_creation_time" == "$cur_xg_creation_time" ]; then
    same_time="yes"
fi
is "$same_time" "no" "autoindexing with --no-guessing overwrites existing files"

rm auto.*

vg autoindex -p auto -w map -r small/x.fa -v small/x.vcf.gz -r small/y.fa -v small/y.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with chunked input"
vg sim -x auto.xg -n 20 -a -l 10 | vg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "chunked autoindexing results can be used by vg map"

rm auto.*

vg autoindex -p auto -w map -r small/xy.fa -v small/xy2.vcf.gz
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with phased input"
vg sim -x auto.xg -n 20 -a -l 10 | vg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "phased autoindexing results can be used by vg map"

rm auto.*

vg autoindex -p auto -w mpmap -w rpvg -r tiny/tiny.fa -v tiny/tiny.vcf.gz -x tiny/tiny.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with unchunked input"
is $(ls auto.* | wc -l) 6 "autoindexing creates 6 files for mpmap/rpvg"
vg sim -x auto.spliced.xg -n 20 -a -l 10 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg mpmap"
is $(vg paths -g auto.haplotx.gbwt -L | wc -l) 6 "haplotype transcript GBWT made by autoindex is valid"
is $(cat auto.txorigin.tsv | wc -l) 7 "transcript origin table has expected number of rows" 

rm auto.*

vg autoindex -p auto -w mpmap  -r tiny/tiny.fa -v tiny/tiny.vcf.gz -x tiny/tiny.gtf --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with unchunked, unphased input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 files for mpmap/rpvg"
vg sim -x auto.spliced.xg -n 20 -a -l 10 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "basic unphased autoindexing results can be used by vg mpmap"

rm auto.*

vg autoindex -p auto -w mpmap -r tiny/tiny.fa -x tiny/tiny.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap without variants"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 files for mpmap without variants"
vg sim -x auto.spliced.xg -n 20 -a -l 10 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with no variants can be used by vg mpmap"

rm auto.*

vg autoindex -p auto -w mpmap -w rpvg -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/x.gtf -x small/y.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with chunked input"
is $(ls auto.* | wc -l) 6 "autoindexing creates 6 files for mpmap/rpvg with chunked input"

rm auto.*

vg autoindex -p auto -w mpmap -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/x.gtf -x small/y.gtf --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with unphased chunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 files for mpmap/rpvg with chunked input"
vg sim -x auto.spliced.xg -n 20 -a -l 10 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with chunked unphased input can be used by vg mpmap"

rm auto.*

vg autoindex -p auto -w mpmap -w rpvg -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/xy.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with partially chunked input"

rm auto.*

vg autoindex -p auto -w mpmap -w rpvg -r small/xy.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/xy.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with another partially chunked input"

rm auto.*

vg autoindex -p auto -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with GFA input"
vg sim -x auto.spliced.xg -n 20 -a -l 5 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with GFA input can be used by vg mpmap"

rm auto.*

vg autoindex -p auto -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf -H graphs/gfa_with_w_lines_haplo.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with GFA input and a haplotype GTF"
vg sim -x auto.spliced.xg -n 20 -a -l 5 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "autoindexing results with GFA input and a haplotype GTF can be used by vg mpmap"

rm auto.*

vg autoindex -p auto -w rpvg -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap and rpvg with GFA input"

rm auto.*

vg autoindex -p auto -w rpvg -w mpmap -g graphs/gfa_with_w_lines.gfa -x graphs/gfa_with_w_lines.gtf -H graphs/gfa_with_w_lines_haplo.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap and rpvg with GFA input and a haplotype GTF"

rm auto.*

vg autoindex -p auto -w giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for vg giraffe"
vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
vg index -x t.xg t.vg
vg sim -x t.xg -n 20 -a -l 10 | vg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -z auto.shortread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg giraffe"

rm auto.*
rm t.*

vg autoindex -p auto -w sr-giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for short read vg giraffe"
vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
vg index -x t.xg t.vg
vg sim -x t.xg -n 20 -a -l 10 | vg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -z auto.shortread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg giraffe"

rm auto.*
rm t.*

vg autoindex -p auto -w lr-giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for long read vg giraffe"
vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
vg index -x t.xg t.vg
vg sim -x t.xg -n 20 -a -l 10 | vg giraffe -Z auto.giraffe.gbz -m auto.longread.withzip.min -z auto.longread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg giraffe"

rm auto.*
rm t.*

vg autoindex -p auto -w giraffe -g graphs/gfa_with_w_lines.gfa 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg giraffe with GFA input with W-lines"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for vg giraffe from GFA input"
vg convert -g graphs/gfa_with_w_lines.gfa -x > g.xg
vg sim -x g.xg -n 20 -a -l 4 | vg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -z auto.shortread.zipcodes -d auto.dist -G - > /dev/null
is $(echo $?) 0 "GFA autoindexing results can be used by vg giraffe"

rm auto.*
rm g.xg

vg autoindex -p auto -w giraffe -g graphs/named_with_walk.gfa 
is $(echo $?) 0 "autoindexing successfully completes on a GFA with named segments and W-lines"
printf '@read\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATTACACATTAGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n' > read.fq
vg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -d auto.dist -f read.fq --named-coordinates > read.gam
is "$(vg view -aj read.gam | jq -r '.path.mapping[].position.name')" "Ishmael" "GFA segment names are available in output GAM when a walk exists"

rm auto.*
rm read.fq read.gam

vg autoindex -p auto -w giraffe -g graphs/named.gfa 
is $(echo $?) 0 "autoindexing successfully completes on a GFA with named segments and no W-lines"
printf '@read\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATTACACATTAGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n' > read.fq
vg giraffe -Z auto.giraffe.gbz -m auto.shortread.withzip.min -d auto.dist -f read.fq --named-coordinates > read.gam
is "$(vg view -aj read.gam | jq -r '.path.mapping[].position.name')" "Ishmael" "GFA segment names are available in output GAM when no walks exist"

# reduce disk limit to 1MB to trigger it during k-mer generation
is "$(vg autoindex -p auto -w map --gcsa-size-limit 1000000 -g graphs/linked_cycles.gfa 2>&1 | grep Rewind | wc -l)" 1 "Running out of room during k-mer enumeration triggers a rewind"
is "$(echo $?)" 0 "Indexing is successful after rewinding from k-mer generation"

rm -f auto.gcsa auto.gcsa.lcp

# reduce disk limit to 2MB to trigger it during doubling steps
is "$(vg autoindex -p auto -w map --gcsa-size-limit 2000000 -g graphs/linked_cycles.gfa 2>&1 | grep Rewind | wc -l)" 1 "Running out of room during GCSA2 indexing triggers a rewind"
is "$(echo $?)" 0 "Indexing is successful after rewinding from GCSA2 indexing"

rm -f auto.gcsa auto.gcsa.lcp

# use the memory limit to trigger a rewide
is "$(vg autoindex -p auto -w map -M 512M -g graphs/linked_cycles.gfa 2>&1 | grep Rewind | wc -l)" 1 "Running out of memory during GCSA2 indexing triggers a rewind"
is "$(echo $?)" 0 "Indexing is successful after rewinding from GCSA2 indexing"

rm auto.*
rm read.fq read.gam

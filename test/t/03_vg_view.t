#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 13

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -d - | wc -l) 505 "view produces the expected number of lines of dot output"
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - | wc -l) 641 "view produces the expected number of lines of GFA output"
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view - | head | vg view -Fv - | vg view - | wc -l) 10 "view converts back and forth between GFA and vg format"

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -a - | wc -l) $(samtools view -u minigiab/NA12878.chr22.tiny.bam | samtools view - | wc -l) "view can convert BAM to GAM"

is "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | jq -c --sort-keys . | sort | md5sum)" "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | vg view -JGa - | vg view -aj - | jq -c --sort-keys . | sort | md5sum)" "view can round-trip JSON and GAM"

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg view -j x.vg | jq . | vg view -Jv - | diff x.vg -
is $? 0 "view can reconstruct a VG graph from JSON"

vg view -v x.vg | cmp -s - x.vg
is $? 0 "view can pass through VG"

rm -f x.vg

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -a - | jq .sample_name | grep -v ^\"1\"$ | wc -l ) 0 "view parses sample names"

is $(vg view -f ./small/x.fa_1.fastq  ./small/x.fa_2.fastq | vg view -a - | wc -l) 2000 "view can handle fastq input"

is $(vg view -Jv ./cyclic/two_node.json | vg view -j - | jq ".edge | length") 4 "view can translate graphs with 2-node cycles"

is $(vg view -g ./cyclic/all.vg | tr '\t' ' ' | grep "4 + 4 -" | wc -l) 1 "view outputs properly oriented GFA"

is $(vg view -d ./cyclic/all.vg | wc -l) 23 "view produces the expected number of lines of dot output from a cyclic graph"

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
is $(cat x.vg x.vg x.vg x.vg | vg view -c - | wc -l) 4 "streaming JSON output produces the expected number of chunks"
rm x.vg

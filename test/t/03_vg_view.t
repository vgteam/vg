#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 20

is $(vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | vg view -d - | wc -l) 505 "view produces the expected number of lines of dot output"
is $(vg construct -m 1000 -r small/x.fa -v small/x.vcf.gz | vg view -g - | wc -l) 503 "view produces the expected number of lines of GFA output"
# This one may throw warnings related to dangling edges because we just take an arbitrary subset of the GFA
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view - | head | vg view -Fv - | vg view - | wc -l) 10 "view converts back and forth between GFA and vg format"

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -a - | wc -l) $(samtools view -u minigiab/NA12878.chr22.tiny.bam | samtools view - | wc -l) "view can convert BAM to GAM"

is "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | jq -c --sort-keys . | sort | md5sum)" "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | vg view -JGa - | vg view -aj - | jq -c --sort-keys . | sort | md5sum)" "view can round-trip JSON and GAM"

# Look at only forward strand reads here to avoid having to reverse qualities.
is "$(samtools view -F 16 minigiab/NA12878.chr22.tiny.bam | head -n1 | cut -f11)" "$(samtools view -F 16 -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | vg view -JGa - | vg view -aj - | head -n1 | jq -r '.quality' | base64 -d | tr '\000-\077' '\041-\150')" "view can round-trip qualities and encodes tham as PHRED+0 base-64-encoded strings in JSON"

# We need to run through GFA because vg construct doesn't necessarily chunk the
# graph the way vg view wants to. We also need to treat as vg::VG to preserve ordering.
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -Vg - | vg view -Fv - >x.vg
vg view -Vj x.vg | jq . | vg view -Jv - | diff x.vg -
is $? 0 "view can reconstruct a VG graph from JSON"

vg view -vV x.vg | cmp -s - x.vg
is $? 0 "view can pass through VG when loading as vg::VG"

vg view x.vg | sort > x.gfa.sorted
vg view -v x.vg | vg view - | sort | cmp -s - x.gfa.sorted
is $? 0 "view can pass through semantically identical VG normally"

rm -f x.vg x.gfa.sorted

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -a - | jq .sample_name | grep -v '^"1"$' | wc -l ) 0 "view parses sample names"

is $(vg view -f ./small/x.fa_1.fastq  ./small/x.fa_2.fastq | vg view -a - | wc -l) 2000 "view can handle fastq input"

is $(vg view -Jv ./cyclic/two_node.json | vg view -j - | jq ".edge | length") 4 "view can translate graphs with 2-node cycles"

is $(vg view -g ./cyclic/all.vg | tr '\t' ' ' | grep "4 + 4 -" | wc -l) 1 "view outputs properly oriented GFA"

is $(vg view -d ./cyclic/all.vg | wc -l) 23 "view produces the expected number of lines of dot output from a cyclic graph"

# We need to make a single-chunk graph
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -v - >x.vg
is $(cat x.vg x.vg x.vg x.vg | vg view -c - | wc -l) 4 "streaming JSON output produces the expected number of chunks"

is "$(cat x.vg x.vg | vg view -vVD - 2>&1 > /dev/null | wc -l)" 0 "duplicate warnings can be suppressed when loading as vg::VG"

rm x.vg

vg view -Fv overlaps/two_snvs_assembly1.gfa >/dev/null 2>errors.txt
is "${?}" "1" "gfa graphs with overlaps are rejected"
is "$(cat errors.txt | wc -l)" "2" "GFA import produces a concise error message when overlaps are present"

vg view -Fv overlaps/incorrect_overlap.gfa >/dev/null 2>errors.txt
is "$?" "1" "GFA import rejects a GFA file with an overlap that goes beyond its sequences"
is "$(cat errors.txt | wc -l)" "2" "GFA import produces a concise error message in that case"

rm -f errors.txt



#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 4

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -d - | wc -l) 504 "view produces the expected number of lines of dot output"
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - | wc -l) 641 "view produces the expected number of lines of GFA output"
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view - | head | vg view -v - | vg view - | wc -l) 10 "view converts back and forth between GFA and vg format"
is $(vg construct -r small/x.fa -v small/x.vcf.gz >/dev/null -P x; vg view -p x | wc -l) 139 "view can dump paths in tabular format"

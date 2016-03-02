#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 7

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -t -r 'http://example.org' - | wc -l) 90 "vg view produces the expected number of lines of turtle"
is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -t -r 'http://example.org/' - | vg view -t -T -r  'http://example.org/' - | wc -l) 90 "vg view produces the expected number of lines of turtle"
is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -t -r 'http://example.org/' - | rapper -c --input turtle -I "http://example.org/vg" -; echo $?) 0 "rapper passed"
is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -tC -r 'http://example.org/' - | rapper -c --input turtle -I "http://example.org/vg" -; echo $?) 0 "rapper passed"
is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -tC -r 'http://example.org' - | wc -l) 5 "vg view produces the expected number of lines of turtle"
is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -tC -r 'http://example.org/' - | vg view -tC -T -r  'http://example.org/' - | wc -l) 5 "vg view produces the expected number of lines of turtle"

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg view -tC -r 'http://example.org/' - | roqet -e "SELECT (COUNT (DISTINCT ?s) AS ?c) WHERE  {?s ?p ?o}" -F turtle -D /dev/stdin -q -r tsv | grep -c 25) 1 "There are 25 distinct subjects in the tiny.ttl tested with SPARQL"

#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 9

# Make sure there's no existing index or its reads will get scooped up.
rm -f tiny.gam.index

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.vg.xg -g tiny.vg.gcsa -k 16 tiny.vg
vg sim -a -n 100 -x tiny.vg.xg -l 30 > reads.gam
vg map -G reads.gam -x tiny.vg.xg -g tiny.vg.gcsa > tiny.gam
vg index -d tiny.gam.index -N tiny.gam
# Since there are no edits, the graph will not change at all.
# Snarls are 1-6, 6-9, 9-12, and 12-15.
vg genotype tiny.vg tiny.gam.index > tiny.loci
is $(vg locify -g tiny.gam.index -x tiny.vg.xg -l tiny.loci -f -n -s loci.sorted | vg view -a - | wc -l) 100 "locify produces output for each input alignment"
is $(cat loci.sorted | wc -l) $(vg stats -R tiny.vg | grep ultrabubble | wc -l) "the sorted list of loci has one locus per snarl"
is $(head -1 loci.sorted) "1+0_6+0" "the first locus is as expected"
is $(head -2 loci.sorted | tail -1) "6+0_9+0" "a middle locus is as expected"
is $(tail -1 loci.sorted) "12+0_15+0" "the last locus is as expected"
rm -rf tiny.gam.index

vg construct -r tiny/tiny.fa -v tiny/multi.vcf.gz >tiny.vg
# This has snarls 1-7, 7-12, 12-15, and 15-19
vg index -x tiny.vg.xg -g tiny.vg.gcsa -k 16 tiny.vg
vg map -G tiny/tiny-s1337-n500-l30.gam -x tiny.vg.xg -g tiny.vg.gcsa > tiny.gam
vg index -d tiny.gam.index -N tiny.gam
vg genotype tiny.vg tiny.gam.index >tiny.loci
is $(vg locify -g tiny.gam.index -b 2 -x tiny.vg.xg -l tiny.loci -f -n -s loci.sorted | vg view -a -| jq -c '.locus[] | [.name, .allele[].name]' | grep "7+0_12+0" | sort | uniq | wc -l) 2 "limitation to 2-best works"
is $(vg locify -g tiny.gam.index -b 3 -x tiny.vg.xg -l tiny.loci -f -n -s loci.sorted | vg view -a -| jq -c '.locus[] | [.name, .allele[].name]' | grep "7+0_12+0" | sort | uniq | wc -l) 3 "limitation to 3-best works"
is $(vg locify -g tiny.gam.index -b 4 -x tiny.vg.xg -l tiny.loci -f -n -s loci.sorted | vg view -a -| jq -c '.locus[] | [.name, .allele[].name]' | grep "7+0_12+0" | sort | uniq | wc -l) 4 "limitation to 4-best works"

vg locify -g tiny.gam.index -b 2 -x tiny.vg.xg -l tiny.loci -f -n -o out.loci >/dev/null
is $(vg view -q out.loci | jq '.allele | length' | sort | uniq -c | wc -l) 1 "we always get one allele when all the reads match the graph"

rm -rf tiny.vg tiny.vg.xg tiny.vg.gcsa tiny.vg.gcsa.lcp tiny.vg tiny.gam tiny.gam.index tiny.loci loci.sorted out.loci

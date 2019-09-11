#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 15

vg construct -m 1000 -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam -e
is $(vg pack -x flat.xg -di 2snp.gam.cx -e | tail -n+2 | cut -f 5 | grep -v ^0$ | wc -l) 2 "allele observation packing detects 2 SNPs"

vg augment -a pileup flat.vg 2snp.gam -P 2snp.gam.vgpu >/dev/null
is $(vg view -l 2snp.gam.vgpu|  jq '.node_pileups[].base_pileup[] | (.num_bases // 0)' | awk '{ print NR-1, $0 }' | head | md5sum | cut -f 1 -d\ )\
   $(vg pack -x flat.xg -di 2snp.gam.cx -e | awk '{ print $3, $4 }' | tail -n+2 | head | md5sum | cut -f 1 -d\ ) "pileup packs agree with graph coverage"
   
vg pack -x flat.xg -i 2snp.gam.cx -i 2snp.gam.cx -i 2snp.gam.cx -o 2snp.gam.cx.3x

is $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx.3x | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")"/3 | bc) $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")" | bc) "graph coverages are merged from multiple .cx indexes"

is $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx.3x | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")"/3 | bc) $(echo "("$(vg pack -x 2snp.xg -di 2snp.gam.cx | tail -n+2 | awk '{ print $4 }' | paste -sd+ - )")" | bc) "edit records are merged from multiple .cx indexes"

x=$(vg pack -x flat.xg -di 2snp.gam.cx | wc -c )
vg pack -x flat.xg -o 2snp.gam.cx -b 10 -g 2snp.gam
y=$(vg pack -x flat.xg -di 2snp.gam.cx | wc -c )
is $x $y "binned edit accumulation does not affect the result"

x=$(vg pack -x flat.xg -di 2snp.gam.cx -n 1 | wc -c)
y=$(vg pack -x flat.xg -di 2snp.gam.cx | wc -c)
is $x $y "pack records are filtered by node id"

x=$(vg pack -x flat.xg -di 2snp.gam.cx -n 1 | cut -f 2 | grep -v "1" | wc -c)
y=$(vg pack -x flat.xg -di 2snp.gam.cx | cut -f 2 | head -n 1 | wc -c)
is $x $y "pack records are filtered by node id"

vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx.3x -i 2snp.gam.cx -i 2snp.gam.cx -i 2snp.gam.cx
x=$(vg pack -x flat.xg -di 2snp.gam.cx.3x | wc -c)
cat 2snp.gam 2snp.gam 2snp.gam | vg pack -x flat.xg -o 2snp.gam.cx -g -
y=$(vg pack -x flat.xg -di 2snp.gam.cx.3x | wc -c)

is $x $y "pack index merging produces the expected result"

vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx.3x -i 2snp.gam.cx -i 2snp.gam.cx -i 2snp.gam.cx
x=$(vg pack -x flat.xg -Di 2snp.gam.cx.3x | wc -c)
cat 2snp.gam 2snp.gam 2snp.gam | vg pack -x flat.xg -o 2snp.gam.cx -g -
y=$(vg pack -x flat.xg -Di 2snp.gam.cx.3x | wc -c)

is $x $y "pack index merging produces the expected result for edges"

rm -f flat.vg 2snp.vg 2snp.xg 2snp.sim flat.gcsa flat.gcsa.lcp flat.xg 2snp.xg 2snp.gam 2snp.gam.cx 2snp.gam.cx.3x 2snp.gam.vgpu

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > tiny.vg
vg index tiny.vg -x tiny.xg
vg sim -l 10 -x tiny.xg -n 50 -a  -s 23 > tiny.gam
vg augment -a pileup tiny.vg tiny.gam -P tiny.vgpu > /dev/null
x=$(vg view -l tiny.vgpu | jq  '.edge_pileups' | grep num_reads | awk '{print $2}' | sed -e 's/\,//' | awk '{sum+=$1} END {print sum}')
y=$(vg pack -x tiny.xg -g tiny.gam -D -o tiny.pack | grep -v from | awk '{sum+=$5} END {print sum}')
is $x $y "pack computes the same total edge coverage as pileup"
x=$(vg pack -x tiny.xg -i tiny.pack -D | grep -v from | awk '{sum+=$5} END {print sum}')
is $x $y "pack stores the correct edge pileup to disk"

rm -f tiny.vg tiny.xg tiny.gam tiny.vgpu tiny.pack

vg construct -m 20 -r tiny/tiny.fa >flat.vg
printf '@forward\nCAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG\n+\n<B<BBB!BBBB<BBBBBBBBBBBBBBBBBBB<BBBBBBBBBBBBB<B7BB\n' > reads.fq
printf '@reverses\nCAGAGAGTTGGAATATAATAGAACTCCAGAAAATTTCCAAGCCTTATTTG\n+\nBB7B<BBBBBBBBBBBBB<BBBBBBBBBBBBBBBBBBB<BBBB!BBB<B<\n' >> reads.fq
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -f reads.fq -k 8 > reads.gam
vg pack -x flat.xg -o reads.gam.cx -g reads.gam -q
is $(vg pack -x flat.xg -di reads.gam.cx | tail -n+2 | cut -f 4 | grep ^0$ | wc -l) 1 "qual-adjust packing detects 1 base with 0 quality"
is $(vg pack -x flat.xg -Di reads.gam.cx | tail | cut -f 5 | grep ^59$ | wc -l) 1 "qual-adjust packing gets correct edge support"

rm -f flat.vg flat.xg flat.gcsa reads.fq reads.gam reads.gam.cx

vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index -x x.xg x.vg
vg sim -s 1 -n 1000 -l 150 -x x.xg -a > sim.gam
vg pack -x x.xg -g sim.gam -o x.xg.cx
vg pack -x x.vg -g sim.gam -o x.vg.cx
vg pack -x x.xg -i x.xg.cx -d | awk '!($1="")' | sort > node-table.xg.tsv
vg pack -x x.vg -i x.vg.cx -d | awk '!($1="")' | sort > node-table.vg.tsv
diff node-table.xg.tsv node-table.vg.tsv
is "$?" 0 "node packs on vg same as xg"

vg pack -x x.xg -i x.xg.cx -D | sort > edge-table.xg.tsv
vg pack -x x.vg -i x.vg.cx -D | sort > edge-table.vg.tsv
diff edge-table.xg.tsv edge-table.vg.tsv
is "$?" 0 "edge packs on vg same as xg"

rm -f x.vg x.xg sim.gam x.xg.cx x.vg.cx node-table.vg.tsv node-table.xg.tsv edge-table.vg.tsv edge-table.xg.tsv

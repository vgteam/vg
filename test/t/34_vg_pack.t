#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 20

vg construct -m 1000 -r tiny/tiny.fa >flat.vg
vg view flat.vg| sed 's/CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG/CAAATAAGGCTTGGAAATTTTCTGGAGATCTATTATACTCCAACTCTCTG/' | vg view -Fv - >2snp.vg
vg index -x 2snp.xg 2snp.vg
vg sim -l 30 -x 2snp.xg -n 30 -a >2snp.sim
vg index -x flat.xg -g flat.gcsa -k 16 flat.vg
vg map -g flat.gcsa -x flat.xg -G 2snp.sim -k 8 >2snp.gam
vg pack -x flat.xg -o 2snp.gam.cx -g 2snp.gam -e
is $(vg pack -x flat.xg -di 2snp.gam.cx -e | tail -n+2 | cut -f 5 | grep -v ^0$ | wc -l) 2 "allele observation packing detects 2 SNPs"

# we replace the comparison to the pileup output, to a snapshot of what it used to be (as pileup is gone)
vg pack -x flat.xg -o 2snp.gam.snapshot.cx -g pileup/2snp.gam -e
is $(cat pileup/2snp.gam.vgpu.json | jq '.node_pileups[].base_pileup[] | (.num_bases // 0)' | awk '{ print NR-1, $0 }' | head | md5sum | cut -f 1 -d\ )\
   $(vg pack -x flat.xg -di 2snp.gam.snapshot.cx -e | awk '{ print $3, $4 }' | tail -n+2 | head | md5sum | cut -f 1 -d\ ) "pileup packs agree with graph coverage"
rm -f 2snp.gam.snapshot.cx
   
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
# now that pileup doesnt exist, we use the snapshots in pileup/ instead of recomputing
x=$(cat pileup/tiny.vgpu.json | jq  '.edge_pileups' | grep num_reads | awk '{print $2}' | sed -e 's/\,//' | awk '{sum+=$1} END {print sum}')
y=$(vg pack -x tiny.xg -g pileup/tiny.gam -D -o tiny.pack | grep -v from | awk '{sum+=$5} END {print sum}')
is $x $y "pack computes the same total edge coverage as pileup"
x=$(vg pack -x tiny.xg -i tiny.pack -D | grep -v from | awk '{sum+=$5} END {print sum}')
is $x $y "pack stores the correct edge pileup to disk"

rm -f tiny.vg tiny.xg tiny.gam tiny.vgpu tiny.pack

vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index -x x.xg x.vg
vg sim -s 1 -n 1000 -l 150 -x x.xg -a > sim.gam
vg pack -x x.xg -g sim.gam -o x.xg.cx
vg pack -x x.vg -g sim.gam -o x.vg.cx -t 1
vg pack -x x.xg -i x.xg.cx -d | awk '!($1="")' | sort > node-table.xg.tsv
vg pack -x x.vg -i x.vg.cx -d | awk '!($1="")' | sort > node-table.vg.tsv
diff node-table.xg.tsv node-table.vg.tsv
is "$?" 0 "node packs on vg same as xg"

vg pack -x x.xg -i x.xg.cx -D | sort > edge-table.xg.tsv
vg pack -x x.vg -i x.vg.cx -D | sort > edge-table.vg.tsv
diff edge-table.xg.tsv edge-table.vg.tsv
is "$?" 0 "edge packs on vg same as xg"

vg pack -x x.vg -g sim.gam -d -t 3 | awk '!($1="")' | sort > node-table.vg.t3.tsv
diff node-table.vg.tsv node-table.vg.t3.tsv
is "$?" 0 "node packs same on vg when using 3 threads as when using 1"

vg pack -x x.vg -g sim.gam -D -t 2 | sort > edge-table.vg.t3.tsv
diff edge-table.vg.tsv edge-table.vg.t3.tsv
is "$?" 0 "edge packs same on vg when using 2 threads as when using 1"

vg convert x.vg -G sim.gam | bgzip | vg pack -x x.vg -a - -o x.vg.gaf.cx
vg pack -x x.vg -i x.vg.gaf.cx -d | awk '!($1="")' | sort > node-table.vg.gaf.tsv
diff node-table.vg.gaf.tsv node-table.vg.tsv
is "$?" 0 "node packs on gaf same as gam"

vg pack -x x.vg -i x.vg.gaf.cx -D | sort > edge-table.vg.gaf.tsv
diff edge-table.vg.gaf.tsv edge-table.vg.tsv
is "$?" 0 "edge packs on gaf same as gam"

rm -f x.vg x.xg sim.gam x.xg.cx x.vg.cx node-table.vg.tsv node-table.xg.tsv edge-table.vg.tsv edge-table.xg.tsv edge-table.vg.t3.tsv node-table.vg.t3.tsv x.vg.gaf.cx node-table.vg.gaf.tsv edge-table.vg.gaf.tsv

vg construct -m 5 -r tiny/tiny.fa >flat.vg
vg index flat.vg -g flat.gcsa
# add a 20 mapq to nodes 1 and 2
vg map -x flat.vg -g flat.gcsa -s CAAATAAGG | vg view -a - | sed -e 's/60/20/g' | vg view -JaG - > flat.gam
# add a 10 mapq to nodes 2 and 3 and 4
vg map -x flat.vg -g flat.gcsa -s GGCTTGGAA | vg view -a - | sed -e 's/60/10/g' | vg view -JaG - >> flat.gam
# add a 60 mapq to node 9 and 10
vg map -x flat.vg -g flat.gcsa -s AACTCTCTG | vg view -a - | vg view -JaG - >> flat.gam
vg pack -x flat.vg -o flat.cx -g flat.gam
is $(vg pack -x flat.vg -i flat.cx -u | awk ' NR>1 {print $2 "\t" $3}' | sort -g | awk '{print $2}' | tr '\n' '-') 20-15-10-10-0-0-0-0-60-60- "average node qualities are correct"

rm -f flat.gam flat.cx 

vg map -x flat.vg -g flat.gcsa -s CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG > span2.gam
vg map -x flat.vg -g flat.gcsa -s CAGAGAGTTGGAATATAATAGAACTCCAGAAAATTTCCAAGCCTTATTTG >> span2.gam
vg pack -x flat.vg -g span2.gam -o span2.pack
vg pack -x flat.vg -g span2.gam -s 10 -o span2.s10.pack
rm -f span2.tsv
for i in `seq 10`; do echo 0 >> span2.tsv ; done
vg pack -x flat.vg -i span2.pack -d | sort -n -k 2 | tail -40 | head -30 | awk '{print $4}' >> span2.tsv
for i in `seq 10`; do echo 0 >> span2.tsv ; done
vg pack -x flat.vg -i span2.s10.pack -d | sort -n -k 2 | awk '{print $4}' | grep -v coverage >  span2.s10.tsv
diff span2.tsv span2.s10.tsv
is "$?" 0 "pack -s 10 sets fist and last 10bp of coverage to 0 like expected"

# another quick test to check if results consisten when node size deoesn't exactly divide length
vg pack -x flat.vg -g span2.gam -s 9 -o span2.s9.pack
rm -f span2.tsv
for i in `seq 9`; do echo 0 >> span2.tsv ; done
vg pack -x flat.vg -i span2.pack -d | sort -n -k 2 | tail -41 | head -32 | awk '{print $4}' >> span2.tsv
for i in `seq 9`; do echo 0 >> span2.tsv ; done
vg pack -x flat.vg -i span2.s9.pack -d | sort -n -k 2 | awk '{print $4}' | grep -v coverage >  span2.s9.tsv
diff span2.tsv span2.s9.tsv
is "$?" 0 "pack -s 9 sets fist and last 9bp of coverage to 0 like expected"

rm -f flat.vg flat.gcsa span2.gam span2.pack span2.s10.pack span2.tsv span2.s10.tsv span2.s9.tsv

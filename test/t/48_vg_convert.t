#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 13

vg construct -r complex/c.fa -v complex/c.vcf.gz > c.vg
cat <(vg view c.vg | grep ^S | sort) <(vg view c.vg | grep L | uniq | wc -l) <(vg paths -v c.vg -E) > c.info

vg convert c.vg -x > c.xg
vg convert c.xg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout xg conversion"

rm -f c.xg c1.vg c1.info

vg convert c.vg -a > c.hg
vg convert c.hg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout hash-graph conversion"

rm -f c.hg c1.vg c1.info

vg convert c.vg -p > c.pg
vg convert c.pg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout packed-graph conversion"

rm -f c.pg c1.vg c1.info

vg convert c.vg -o > c.odgi
vg convert c.odgi -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout ODGI conversion"

rm -f c.vg c.odgi c1.vg c.info c1.info

# some less rigorous tests I made without noticing that the earlier ones had already been written
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg view x.vg > x.gfa

is "$(vg convert -a x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "hash graph conversion looks good"
is "$(vg convert -p x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "packed graph conversion looks good"
is "$(vg convert -v x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "vg conversion looks good"
is "$(vg convert -o x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "odgi conversion looks good"
is "$(vg convert -x x.vg | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "xg conversion looks good"

is "$(vg convert -g -a x.gfa | vg view - | wc -l)" "$(wc -l < x.gfa)" "on disk gfa conversion looks good"
is "$(cat x.gfa | vg convert -g -a - | vg view - | wc -l)" "$(wc -l < x.gfa)" "streaming gfa conversion looks good"
is "$(vg convert -g -x x.gfa | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "gfa to xg conversion looks good"

rm x.vg x.gfa
rm -f c.vg c.pg c1.vg c.info c1.info

vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index x.vg -g x.gcsa
vg sim -x x.vg -n 10 -s 23 -a > sim.gam
vg map -x x.vg -g x.gcsa -G sim.gam > sim-rm.gam
vg convert x.vg -G sim-rm.gam > sim-rm.gaf
vg convert x.vg -F sim-rm.gaf | vg convert x.vg -G - > sim-rm2.gaf
diff sim-rm.gaf sim-rm2.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf makes same gaf twice"

rm -f x.vg x.gcsa sim.gam sim-rm.gam sim-rm.gaf sim-rm2.gaf

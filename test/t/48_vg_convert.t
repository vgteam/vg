#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 3

vg construct -r complex/c.fa -v complex/c.vcf.gz > c.vg

vg convert c.vg -x > c.xg
vg convert c.xg -v > c1.vg
cat <(vg view c.vg | grep ^S | sort) <(vg view c.vg | grep L | uniq | wc -l) <(vg paths -v c.vg -E) > c.info
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout xg conversion"

rm -f c.xg c1.vg c.info c1.info

vg convert c.vg -a > c.hg
vg convert c.hg -v > c1.vg
cat <(vg view c.vg | grep ^S | sort) <(vg view c.vg | grep L | uniq | wc -l) <(vg paths -v c.vg -E) > c.info
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout hash-graph conversion"

rm -f c.hg c1.vg c.info c1.info

vg convert c.vg -p > c.pg
vg convert c.pg -v > c1.vg
cat <(vg view c.vg | grep ^S | sort) <(vg view c.vg | grep L | uniq | wc -l) <(vg paths -v c.vg -E) > c.info
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout packed-graph conversion"

rm -f c.vg c.pg c1.vg c.info c1.info


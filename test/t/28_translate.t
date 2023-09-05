#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 2

vg construct -m 1000 -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.xg -g tiny.gcsa -k 16 tiny.vg
vg sim -n 5 -e 0.01 -i 0.005 -x tiny.xg -l 30 -a -s 1 | vg view -a - | sort | vg view -JGa - > tiny.sim
vg map -G tiny.sim -x tiny.xg -g tiny.gcsa -t 1 > tiny.gam
vg augment -Z tiny.trans -i tiny.vg tiny.gam >tiny.mod.vg
vg paths -v tiny.mod.vg -X | vg view -a - | grep -v x | sort | vg view -JGa - >tiny.paths.gam
vg translate -a tiny.paths.gam tiny.trans | vg view -a - | sort | vg view -JGa - >tiny.paths.trans.gam
vg augment -Z tiny.trans.1 -i tiny.vg tiny.paths.trans.gam >tiny.mod.vg.1

is $(vg mod -U 10 tiny.mod.vg | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg mod -U 10 tiny.mod.vg.1 | vg mod -c - | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) "alignments used to modify a graph may be projected back to the original graph and used to regenerate the same graph"

rm -Rf tiny.vg tiny.xg tiny.gcsa tiny.gcsa.lcp tiny.sim tiny.gam tiny.trans tiny.mod.vg tiny.trans.1 tiny.paths.gam tiny.paths.trans.gam tiny.mod.vg.1

vg construct -m 1000 -r tiny/tiny.fa >flat.vg
vg index -x flat.xg -g flat.gcsa -k 8 flat.vg
vg map -x flat.xg -g flat.gcsa -G tiny/flat-s69-n1-l50-e0.05.gam >flat.gam
vg augment -i -Z flat1.trans flat.vg flat.gam >flat1.vg
vg index -x flat1.xg -g flat1.gcsa -k 8 flat1.vg
vg map -x flat1.xg -g flat1.gcsa -G tiny/flat1-s77-n1-l50-e0.05.gam >flat1.gam
vg augment -i -Z flat2.trans flat1.vg flat1.gam >flat2.vg
vg translate -o flat2.trans flat1.trans >flatover.trans
vg paths -v flat2.vg -X | vg view -a - | grep -v x | vg view -JGa - >flat2.paths.gam
vg translate -a flat2.paths.gam flatover.trans >flatback.gam
vg augment -i flat.vg flatback.gam >flat2back.vg
is $(vg view flat2back.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 3ce4390abe667b24eb5148e216a2f5ec "translation overlay works and produces a sane result"

rm -f flat.vg flat.xg flat.gcsa.lcp flat.gcsa flat.gam flat1.vg flat1.trans flat1.xg flat1.gcsa.lcp flat1.gcsa flat1.gam flat2.vg flat2.trans flatover.trans flat2.paths.gam flatback.gam flat2back.vg

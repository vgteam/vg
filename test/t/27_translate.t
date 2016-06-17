#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg index -x tiny.vg.xg tiny.vg
vg sim -s 420 -n 5 -e 0.1 -i 0.2 -x tiny.vg.xg -l 30 -a | vg view -a - | grep -v 6f1cc97082a37070 | sort | vg view -JGa - > tiny.sim
vg map -G tiny.sim -k 8 -V tiny.vg -t 1 > tiny.gam
vg mod -Z tiny.trans -i tiny.gam tiny.vg >tiny.mod.vg
vg paths -x tiny.mod.vg | vg view -a - | grep -v x | sort | vg view -JGa - >tiny.paths.gam
vg translate -a tiny.paths.gam tiny.trans | vg view -a - | sort | vg view -JGa - >tiny.paths.trans.gam
vg mod -Z tiny.trans.1 -i tiny.paths.trans.gam tiny.vg >tiny.mod.vg.1

is $(vg view tiny.mod.vg | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg view tiny.mod.vg.1 | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) "alignments used to modify a graph may be projected back to the original graph and used to regenerate the same graph"

rm -Rf tiny.vg tiny.vg.xg tiny.sim tiny.gam tiny.trans tiny.mod.vg tiny.trans.1 tiny.paths.gam tiny.paths.trans.gam tiny.mod.vg.1




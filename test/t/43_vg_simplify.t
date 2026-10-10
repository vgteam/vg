#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 16

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg

# Make sure to discard all the warnings about having removed alt paths.
vg simplify --algorithm small x.vg > x.small.vg 2>/dev/null
is "${?}" "0" "vg simplify runs through when popping small bubbles"
is "$(vg paths -d -v x.small.vg | vg mod --unchop - | vg stats -N -)" "1" "simplification pops all the bubbles in a simple graph"

# Same test but new path simplifier logic
vg simplify --algorithm small x.vg -P x > x.small.vg 2>/dev/null
is "${?}" "0" "vg path simplify runs through when popping small bubbles"
is "$(vg paths -d -v x.small.vg | vg mod --unchop - | vg stats -N -)" "1" "path simplification pops all the bubbles in a simple graph"

vg simplify --algorithm rare --min-count 2 -v small/x.vcf.gz x.vg > x.rare.vg
is "${?}" "0" "vg simplify runs through when removing rare variants"

vg validate x.rare.vg
is "${?}" "0" "the graph is valid after removing rare variants"

is "$(vg paths -d -v x.rare.vg | vg mod --unchop - | vg stats -N -)" "118" "simplification keeps only some variants"

rm -f x.vg x.small.vg x.rare.vg
 

# -L clustering: a pure deletion used to be unclusterable, so this path could never merge one.
# vg simplify is the only -L consumer that mutates the graph, so the merged allele's nodes and
# edges are actually removed.  -m 0 is required or the length filter fires first.
vg view -Fv nesting/simplify_del_absorbs.gfa > simplify_del.vg
vg simplify --algorithm small -P x -m 0 -L 0.6 simplify_del.vg > simplify_del_L.vg 2>/dev/null
is "$(vg stats -E simplify_del_L.vg)" "3" "-L removes the edge of an allele merged into a pure deletion"
# the node count falls from 4 to 3 as a consequence: with that edge gone unchop fuses the chain.
# the edge assertion above is the one that pins the merge itself.
is "$(vg stats -N simplify_del_L.vg)" "3" "and the resulting chain is unchopped"
vg validate simplify_del_L.vg 2>/dev/null
is "$?" 0 "the graph is still valid after -L simplification"
vg simplify --algorithm small -P x -m 0 -L 1.0 simplify_del.vg > simplify_del_noop.vg 2>/dev/null
is "$(vg stats -N simplify_del_noop.vg)" "4" "-L 1.0 changes nothing"

# The scale a pure deletion is measured against must come from the -P reference, not from whichever
# spanning path sorts first: at a top-level snarl every spanning path ties on rank, so picking the
# first by name made the merge depend on how the haplotypes happen to be named.  del59_vs_del60.gfa
# is simplify_del_absorbs.gfa with its alt paths renamed to sort BEFORE the reference.
vg view -Fv nesting/del59_vs_del60.gfa > simplify_del_ref_last.vg
vg simplify --algorithm small -P x -m 0 -L 0.6 simplify_del_ref_last.vg > simplify_del_rl.vg 2>/dev/null
is "$(vg stats -E simplify_del_rl.vg)" "3" "-L merges the same site when the reference path does not sort first"
is "$(vg paths -L -v simplify_del_rl.vg | grep -c '^x$')" "1" "and the -P reference path survives the merge intact"

rm -f simplify_del.vg simplify_del_L.vg simplify_del_noop.vg simplify_del_ref_last.vg simplify_del_rl.vg

# The -P reference must survive simplification whatever it happens to be called.  The traversal that
# seeds the set of nodes and edges kept is ref_trav_candidates[0], and at a top-level snarl every
# spanning path ties on rank, so that is just the alphabetically-first path.  When a haplotype sorted
# ahead of the reference, the reference's own nodes were deleted and its path left in fragments.
# No -L involved -- this is the plain -m length filter.
R40=$(printf 'A%.0s' $(seq 40)); E40=$(printf 'T%.0s' $(seq 40))
{ printf "H\tVN:Z:1.0\n"
  printf "S\t1\t$R40\nS\t2\tGGGGGGGGG\nS\t3\tCCCCCCCCC\nS\t4\t$E40\n"
  printf "L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n"
  printf "P\tzzREF\t1+,2+,4+\t*,*,*\n"
  printf "P\taAlt\t1+,3+,4+\t*,*,*\n"; } > refsort.gfa
vg simplify --algorithm small -P zzREF -m 10 -k refsort.gfa > refsort.vg 2>/dev/null
is "$(vg paths -E -x refsort.vg | awk '$1=="zzREF"{print $2}')" "89" "the -P reference keeps its full length when a haplotype sorts ahead of it"
is "$(vg paths -L -x refsort.vg | grep -c '^zzREF$')" "1" "and its path is not left in fragments"
is "$(vg stats -N refsort.vg)" "2" "while the alt allele is still collapsed"
rm -f refsort.gfa refsort.vg

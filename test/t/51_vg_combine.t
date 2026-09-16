#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 115

C=combine
mkdir -p $C

# Construct chunks holding a reference plus two haplotypes.
RS='H\tVN:Z:1.1\tRS:Z:GRCh38\n'

# Left chunk including HG002 and HG005#1
printf "$RS"'S\t1\tAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tGGGGG\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+,4+\t*\nP\tHG002#1#chr1#0[0-10]\t1+,2+,4+\t*\nP\tHG002#2#chr1#0[0-10]\t1+,3+,4+\t*\nP\tHG005#1#chr1#0[0-10]\t1+,3+,4+\t*\n' > $C/ml.gfa
# Right chunk including only HG002
printf "$RS"'S\t5\tTTTT\nS\t6\tA\nS\t7\tC\nS\t8\tGGGGG\nL\t5\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t8\t+\t0M\nL\t7\t+\t8\t+\t0M\nP\tGRCh38#0#chr1[10-20]\t5+,6+,8+\t*\nP\tHG002#1#chr1#0[10-20]\t5+,7+,8+\t*\nP\tHG002#2#chr1#0[10-20]\t5+,6+,8+\t*\n' > $C/mr.gfa
vg convert -p $C/ml.gfa > $C/ml.pg 2>/dev/null
vg convert -p $C/mr.gfa > $C/mr.pg 2>/dev/null

is "$(vg paths -M -v $C/ml.pg | tail -n +2 | cut -f2 | grep -c HAPLOTYPE)" "3" "the chunk carries three haplotype paths"
is "$(vg paths -M -v $C/ml.pg | tail -n +2 | cut -f2 | grep -c REFERENCE)" "1" "the chunk carries one reference path"

vg combine $C/ml.pg $C/mr.pg > $C/o_def.vg 2>/dev/null
is "$(vg paths -Ev $C/o_def.vg | wc -l)" "4" "every path is considered"
is "$(vg paths -Ev $C/o_def.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-20],HG002#1#chr1#0[0-20],HG002#2#chr1#0[0-20],HG005#1#chr1#0[0-10]" "each sample and haplotype merges within its own group"
is "$(vg paths -Ev $C/o_def.vg | cut -f2 | sort -n | paste -sd,)" "10,20,20,20" "the merged paths span both chunks and the lone one does not"
is "$(vg stats -N $C/o_def.vg)" "8" "merging keeps every node of both chunks, 4 for each"
is "$(vg stats -E $C/o_def.vg)" "9" "the three paths cross the seam on one shared edge, so the edges count increase by 1(8 -> 9)"
is "$(vg stats -s $C/o_def.vg | wc -l)" "1" "the two chunks end up in one component"

# Only fragments of the same logical paths are merged.
is "$(vg paths -Fv $C/o_def.vg | grep -A1 '^>GRCh38#0#chr1\[0-20\]$' | tail -1)" "AAAACGGGGGTTTTAGGGGG" "the reference only merge its own fragments"
is "$(vg paths -Fv $C/o_def.vg | grep -A1 '^>HG002#1#chr1#0\[0-20\]$' | tail -1)" "AAAACGGGGGTTTTCGGGGG" "HG002#1#chr1 only merge its own fragments"
is "$(vg paths -Fv $C/o_def.vg | grep -A1 '^>HG002#2#chr1#0\[0-20\]$' | tail -1)" "AAAAGGGGGGTTTTAGGGGG" "HG002#2#chr1 only merge its own fragments"
is "$(vg paths -Fv $C/o_def.vg | grep -A1 '^>HG005#1#chr1#0\[0-10\]$' | tail -1)" "AAAAGGGGGG" "HG005#1#chr1 is the lone fragment untouched"

# renumber is the default, so naming it explicitly must change nothing.
vg combine -s renumber $C/ml.pg $C/mr.pg > $C/o_ren.vg 2>/dev/null
is "$(vg paths -Ev $C/o_ren.vg)" "$(vg paths -Ev $C/o_def.vg)" "--seam renumber is the default"
is "$(vg paths -Fv $C/o_ren.vg | sort)" "$(vg paths -Fv $C/o_def.vg | sort)" "and spells out the same sequences"

# A path with no subrange is not a fragment of anything, so its name must come through untouched.
printf "$RS"'S\t1\tAAAACCCC\nP\tGRCh38#0#chr1\t1+\t*\n' > $C/solo.gfa
vg convert -p $C/solo.gfa > $C/solo.pg 2>/dev/null

vg combine $C/solo.pg > $C/o_solo.vg 2>/dev/null
is "$(vg paths -Ev $C/o_solo.vg)" "$(printf 'GRCh38#0#chr1\t8')" "a lone path with no subrange keeps its name"

# Fragments from different haplotype must not be merged.
printf "$RS"'S\t1\tAAAA\nS\t2\tCGGGGG\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+\t*\nP\tHG002#1#chr1#0[0-10]\t1+,2+\t*\n' > $C/h1.gfa
printf "$RS"'S\t3\tTTTT\nS\t4\tAGGGGG\nL\t3\t+\t4\t+\t0M\nP\tGRCh38#0#chr1[10-20]\t3+,4+\t*\nP\tHG002#2#chr1#0[10-20]\t3+,4+\t*\n' > $C/h2.gfa
vg convert -p $C/h1.gfa > $C/h1.pg 2>/dev/null
vg convert -p $C/h2.gfa > $C/h2.pg 2>/dev/null

vg combine $C/h1.pg $C/h2.pg > $C/o_hap.vg 2>/dev/null
is "$(vg paths -Ev $C/o_hap.vg | wc -l)" "3" "only the reference path has fragments to merge"
is "$(vg paths -Ev $C/o_hap.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-20],HG002#1#chr1#0[0-10],HG002#2#chr1#0[10-20]" "fragments of different haplotypes stay apart"

# Shuffled fragments of two contigs must not splice into each other.
printf "$RS"'S\t1\tAAAA\nS\t2\tCCCCCC\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+\t*\n' > $C/mc_chr1l.gfa
printf "$RS"'S\t1\tGGGG\nS\t2\tTTTTTT\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[10-20]\t1+,2+\t*\n' > $C/mc_chr1r.gfa
printf "$RS"'S\t1\tTTTT\nS\t2\tGGGGGG\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr2[0-10]\t1+,2+\t*\n' > $C/mc_chr2l.gfa
printf "$RS"'S\t1\tCCCC\nS\t2\tAAAAAA\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr2[10-20]\t1+,2+\t*\n' > $C/mc_chr2r.gfa
for g in mc_chr1l mc_chr1r mc_chr2l mc_chr2r; do vg convert -p $C/$g.gfa > $C/$g.pg 2>/dev/null; done

vg combine $C/mc_chr2r.pg $C/mc_chr1l.pg $C/mc_chr2l.pg $C/mc_chr1r.pg > $C/o_multi.vg 2>/dev/null
is "$(vg paths -Ev $C/o_multi.vg | wc -l)" "2" "each contig yields its own merged path"
is "$(vg paths -Ev $C/o_multi.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-20],GRCh38#0#chr2[0-20]" "both contigs span both of their chunks"
is "$(vg paths -Ev $C/o_multi.vg | cut -f2 | sort -u)" "20" "neither contig grew by taking the other's chunks"
is "$(vg stats -s $C/o_multi.vg | wc -l)" "2" "the contigs stay in separate components"
is "$(vg stats -E $C/o_multi.vg)" "6" "only the two intra-contig seams are edged"
is "$(vg paths -Fv $C/o_multi.vg | grep -A1 '^>GRCh38#0#chr1\[0-20\]$' | tail -1)" "AAAACCCCCCGGGGTTTTTT" "chr1 spells out only its own chunks"
is "$(vg paths -Fv $C/o_multi.vg | grep -A1 '^>GRCh38#0#chr2\[0-20\]$' | tail -1)" "TTTTGGGGGGCCCCAAAAAA" "chr2 spells out only its own chunks"

# --seam shared: the two chunks sharing nodes 4 and 5.
printf "$RS"'S\t1\tAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tGGGGG\nS\t5\tTTTT\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\nP\tGRCh38#0#chr1[0-14]\t1+,2+,4+,5+\t*\nP\tHG002#1#chr1#0[0-14]\t1+,2+,4+,5+\t*\nP\tHG002#2#chr1#0[0-14]\t1+,3+,4+,5+\t*\n' > $C/sl.gfa
printf "$RS"'S\t4\tGGGGG\nS\t5\tTTTT\nS\t6\tA\nS\t7\tC\nS\t8\tGGGGG\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t8\t+\t0M\nL\t7\t+\t8\t+\t0M\nP\tGRCh38#0#chr1[5-20]\t4+,5+,6+,8+\t*\nP\tHG002#1#chr1#0[5-20]\t4+,5+,7+,8+\t*\nP\tHG002#2#chr1#0[5-20]\t4+,5+,6+,8+\t*\n' > $C/sr.gfa
vg convert -p $C/sl.gfa > $C/sl.pg 2>/dev/null
vg convert -p $C/sr.gfa > $C/sr.pg 2>/dev/null

vg combine -s shared $C/sl.pg $C/sr.pg > $C/o_sh.vg 2>/dev/null
is "$(vg paths -Ev $C/o_sh.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-20],HG002#1#chr1#0[0-20],HG002#2#chr1#0[0-20]" "--seam shared merges all three paths over their shared nodes"
is "$(vg paths -Ev $C/o_sh.vg | cut -f2 | sort -u)" "20" "no path counts the same nodes twice"
is "$(vg stats -N $C/o_sh.vg)" "8" "the shared nodes are deduplicated rather than copied"
is "$(vg stats -E $C/o_sh.vg)" "9" "the edges between them are also deduplicated"
is "$(vg paths -Fv $C/o_sh.vg | grep -A1 '^>HG002#1#chr1#0\[0-20\]$' | tail -1)" "AAAACGGGGGTTTTCGGGGG" "HG002#1#chr1 keeps its own sequences over the shared nodes"
is "$(vg paths -Fv $C/o_sh.vg | grep -A1 '^>HG002#2#chr1#0\[0-20\]$' | tail -1)" "AAAAGGGGGGTTTTAGGGGG" "HG002#2#chr1 keeps its own sequences over the shared nodes"

# --seam shared: different paths share different nodes.
printf "$RS"'S\t4\tGGGGG\nS\t5\tTTTT\nS\t6\tA\nS\t7\tC\nS\t8\tGGGGG\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t8\t+\t0M\nL\t7\t+\t8\t+\t0M\nP\tGRCh38#0#chr1[5-20]\t4+,5+,6+,8+\t*\nP\tHG002#1#chr1#0[5-20]\t4+,5+,7+,8+\t*\nP\tHG002#2#chr1#0[10-20]\t5+,6+,8+\t*\n' > $C/sr2.gfa
vg convert -p $C/sr2.gfa > $C/sr2.pg 2>/dev/null

vg combine -s shared $C/sl.pg $C/sr2.pg > $C/o_sh2.vg 2>/dev/null
is "$(vg paths -Ev $C/o_sh2.vg | cut -f2 | sort -u)" "20" "HG002#2#chr1 has a shorter overlap but it still keeps the same total length"
is "$(vg paths -Fv $C/o_sh2.vg | grep -A1 '^>HG002#2#chr1#0\[0-20\]$' | tail -1)" "AAAAGGGGGGTTTTAGGGGG" "HG002#2#chr1 also keeps the same sequence"

# --seam trim: on a 2bp overlap.
printf "$RS"'S\t1\tAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tGGGGG\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+,4+\t*\nP\tHG002#1#chr1#0[0-10]\t1+,2+,4+\t*\nP\tHG002#2#chr1#0[0-10]\t1+,3+,4+\t*\n' > $C/tl.gfa
printf "$RS"'S\t5\tGGTT\nS\t6\tA\nS\t7\tC\nS\t8\tGGGGG\nL\t5\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t8\t+\t0M\nL\t7\t+\t8\t+\t0M\nP\tGRCh38#0#chr1[8-18]\t5+,6+,8+\t*\nP\tHG002#1#chr1#0[8-18]\t5+,7+,8+\t*\nP\tHG002#2#chr1#0[8-18]\t5+,6+,8+\t*\n' > $C/tr.gfa
vg convert -p $C/tl.gfa > $C/tl.pg 2>/dev/null
vg convert -p $C/tr.gfa > $C/tr.pg 2>/dev/null

vg combine -s trim $C/tl.pg $C/tr.pg > $C/o_trim.vg 2>/dev/null
is "$(vg paths -Ev $C/o_trim.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-18],HG002#1#chr1#0[0-18],HG002#2#chr1#0[0-18]" "all three paths report the same trimmed span"
is "$(vg paths -Ev $C/o_trim.vg | cut -f2 | sort -u)" "18" "the overlap of three paths has been trimmed"
is "$(vg stats -N $C/o_trim.vg)" "8" "trimming splits the boundary node and drops the overlapping half, leaving the node counts unchanged"
is "$(vg stats -E $C/o_trim.vg)" "9" "seaming create one edge at the boundary"
is "$(vg paths -Fv $C/o_trim.vg | grep -A1 '^>GRCh38#0#chr1\[0-18\]$' | tail -1)" "AAAACGGGGGTTAGGGGG" "the overlap of the reference path is trimmed"
is "$(vg paths -Fv $C/o_trim.vg | grep -A1 '^>HG002#1#chr1#0\[0-18\]$' | tail -1)" "AAAACGGGGGTTCGGGGG" "the overlap of the HG002#1#chr1 is trimmed"
is "$(vg paths -Fv $C/o_trim.vg | grep -A1 '^>HG002#2#chr1#0\[0-18\]$' | tail -1)" "AAAAGGGGGGTTAGGGGG" "the overlap of the HG002#2#chr1 is trimmed"

# --seam trim --fuse: must rewrite the step of every path, replace the original nodes with the fused nodes.
vg combine -s trim -u $C/tl.pg $C/tr.pg > $C/o_fuse.vg 2>/dev/null
is "$(vg stats -N $C/o_fuse.vg)" "7" "the two boundary nodes become one"
is "$(vg stats -E $C/o_fuse.vg)" "8" "the edge between two boundary nodes is removed"
is "$(vg paths -Ev $C/o_fuse.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-18],HG002#1#chr1#0[0-18],HG002#2#chr1#0[0-18]" "welding does not change any path's coordinate"
is "$(vg paths -Ev $C/o_fuse.vg | cut -f2 | sort -u)" "18" "welding doesn't change any path's length"
is "$(vg paths -Fv $C/o_fuse.vg | sort)" "$(vg paths -Fv $C/o_trim.vg | sort)" "welding doesn't change any path's sequence"

# --seam trim on two contigs interleaved: each must trim to its own span.
printf "$RS"'S\t1\tAAAAA\nS\t2\tCCCCC\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+\t*\n' > $C/mt_chr1l.gfa
printf "$RS"'S\t1\tCCGGG\nS\t2\tTTTTT\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[8-18]\t1+,2+\t*\n' > $C/mt_chr1r.gfa
printf "$RS"'S\t1\tTTTTT\nS\t2\tGGGGG\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr2[0-10]\t1+,2+\t*\n' > $C/mt_chr2l.gfa
printf "$RS"'S\t1\tGGAAA\nS\t2\tCCCCC\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr2[8-18]\t1+,2+\t*\n' > $C/mt_chr2r.gfa
for g in mt_chr1l mt_chr1r mt_chr2l mt_chr2r; do vg convert -p $C/$g.gfa > $C/$g.pg 2>/dev/null; done

vg combine -s trim $C/mt_chr2r.pg $C/mt_chr1l.pg $C/mt_chr2l.pg $C/mt_chr1r.pg > $C/o_mt.pg 2>/dev/null
is "$(vg paths -Ev $C/o_mt.pg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-18],GRCh38#0#chr2[0-18]" "each contig trims to its own span"
is "$(vg paths -Ev $C/o_mt.pg | cut -f2 | sort -u)" "18" "each contig lost only its own overlap"
is "$(vg stats -s $C/o_mt.pg | wc -l)" "2" "the contigs stay in separate components"
is "$(vg stats -E $C/o_mt.pg)" "6" "only the two intra-contig seams are edged"
is "$(vg paths -Fv $C/o_mt.pg | grep -A1 '^>GRCh38#0#chr1\[0-18\]$' | tail -1)" "AAAAACCCCCGGGTTTTT" "chr1 keeps its own sequence through the trim"
is "$(vg paths -Fv $C/o_mt.pg | grep -A1 '^>GRCh38#0#chr2\[0-18\]$' | tail -1)" "TTTTTGGGGGAAACCCCC" "chr2 keeps its own"

# --fuse welds one seam per contig, not one overall.
vg combine -s trim -u $C/mt_chr2r.pg $C/mt_chr1l.pg $C/mt_chr2l.pg $C/mt_chr1r.pg > $C/o_mtu.pg 2>/dev/null
is "$(vg stats -N $C/o_mtu.pg)" "6" "each contig's boundary nodes weld into one"
is "$(vg stats -E $C/o_mtu.pg)" "4" "and each contig's seam edge goes with them"
is "$(vg paths -Fv $C/o_mtu.pg | sort)" "$(vg paths -Fv $C/o_mt.pg | sort)" "welding leaves both contigs spelling the same sequence"

# Another single chunk.
printf "$RS"'S\t1\tAAAA\nS\t2\tCCC\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[0-7]\t1+,2+\t*\n' > $C/one.gfa
vg convert -p $C/one.gfa > $C/one.pg 2>/dev/null

# Error: nothing to combine.
vg combine -s renumber 2>$C/none.err >/dev/null
is "$?" "1" "combining nothing is rejected"
is "$(grep -c 'At least one input graph' $C/none.err)" "1" "the error says at least one input graph is required"

# Error: --seam trim needs more than one chunk.
vg combine -s trim $C/one.pg 2>$C/lone.err >/dev/null
is "$?" "1" "--seam trim on a single input is rejected"
is "$(grep -c 'at least two input graphs' $C/lone.err)" "1" "the error says trim needs at least two inputs"

# Error: --fuse only makes sense when --seam trim leaves a boundary to weld.
vg combine -u $C/tl.pg $C/tr.pg 2>$C/fusearg.err >/dev/null
is "$?" "1" "--fuse without --seam trim is rejected"
is "$(grep -c 'requires --seam trim' $C/fusearg.err)" "1" "the error says --fuse needs --seam trim"

# Error: --seam only takes the policies it knows.
vg combine -s bogus $C/tl.pg $C/tr.pg 2>$C/seam.err >/dev/null
is "$?" "1" "an unrecognized --seam policy is rejected"
is "$(grep -c 'renumber, shared, trim' $C/seam.err)" "1" "the error lists the policies that do exist"

# Error: vg combine joins pieces from separate inputs, not pieces inside one input.
printf "$RS"'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nP\tGRCh38#0#chr1[0-10]\t1+\t*\nP\tGRCh38#0#chr1[10-20]\t2+\t*\n' > $C/split.gfa
# Same input, plus the edge that already joins the two pieces.
printf "$RS"'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+\t*\nP\tGRCh38#0#chr1[10-20]\t2+\t*\n' > $C/split_joined.gfa
vg convert -p $C/split.gfa > $C/split.pg 2>/dev/null
vg convert -p $C/split_joined.gfa > $C/split_joined.pg 2>/dev/null

vg combine $C/split.pg 2>$C/split.err >/dev/null
is "$?" "1" "an input holding a path in pieces is rejected"
is "$(grep -c 'GRCh38#0#chr1\[0-10\].*GRCh38#0#chr1\[10-20\]' $C/split.err)" "1" "the error reports the pieces"
is "$(grep -c 'split.pg' $C/split.err)" "1" "the error reports the input holding them"

vg combine $C/split_joined.pg 2>/dev/null >/dev/null
is "$?" "1" "and is rejected even when an edge already joins the pieces"

# Error: a read-only graph format cannot be combined; VPKG reports this itself.
vg gbwt -g $C/m.gbz --gbz-format -G graphs/gfa_with_reference.gfa 2>/dev/null
vg combine $C/m.gbz $C/one.pg 2>/dev/null >/dev/null
is "$?" "1" "a GBZ input is rejected rather than crashing"

# Error: --seam trim joins chunks end to end, so a gap between them is rejected.
printf "$RS"'S\t1\tAAAAAAAAAA\nP\tGRCh38#0#chr1[0-10]\t1+\t*\n' > $C/gapl.gfa
printf "$RS"'S\t2\tCCCCCCCCCC\nP\tGRCh38#0#chr1[20-30]\t2+\t*\n' > $C/gapr.gfa
vg convert -p $C/gapl.gfa > $C/gapl.pg 2>/dev/null
vg convert -p $C/gapr.gfa > $C/gapr.pg 2>/dev/null
vg combine -s trim $C/gapl.pg $C/gapr.pg 2>$C/tgap.err >/dev/null
is "$?" "1" "--seam trim rejects a gap between chunks"
is "$(grep -c 'REFERENCE offset gap of 10 bp' $C/tgap.err)" "1" "--seam trim also reports the size of the gap"

# Error: the default policy refuses the same gap, caught by a different check.
vg combine $C/gapl.pg $C/gapr.pg 2>$C/gap.err >/dev/null
is "$?" "1" "the default policy rejects a gap between fragments too"
is "$(grep -c 'Gap of 10 bp' $C/gap.err)" "1" "the default policy reports the size of the gap"

# Error: --seam trim needs at least one reference.
printf 'S\t1\tAAAA\nP\tplain_contig\t1+\t*\n' > $C/noref.gfa
vg convert -p $C/noref.gfa > $C/noref.pg 2>/dev/null
vg combine -s trim $C/noref.pg $C/one.pg 2>$C/noref.err >/dev/null
is "$?" "1" "--seam trim rejects a chunk with no REFERENCE path"
is "$(grep -c 'No REFERENCE-sense path found' $C/noref.err)" "1" "the error says a reference path is needed"

# Error: --seam trim rejects more than one reference.
printf "$RS"'S\t1\tAAAA\nP\tGRCh38#0#chr1[0-4]\t1+\t*\nP\tGRCh38#0#chr2[0-4]\t1+\t*\n' > $C/tworef.gfa
vg convert -p $C/tworef.gfa > $C/tworef.pg 2>/dev/null
vg combine -s trim $C/tworef.pg $C/one.pg 2>$C/tworef.err >/dev/null
is "$?" "1" "--seam trim rejects a chunk with two REFERENCE paths"
is "$(grep -c '2 REFERENCE-sense paths found' $C/tworef.err)" "1" "the error counts the paths it found"

# Error: --seam trim requires trimming on the forward strand of the boundary node.
printf "$RS"'S\t1\tAAAA\nS\t2\tCCC\nL\t1\t-\t2\t+\t0M\nP\tGRCh38#0#chr1[0-7]\t1-,2+\t*\n' > $C/rev.gfa
vg convert -p $C/rev.gfa > $C/rev.pg 2>/dev/null
vg combine -s trim $C/rev.pg $C/one.pg 2>$C/rev.err >/dev/null
is "$?" "1" "--seam trim rejects a reverse-strand boundary"
is "$(grep -c 'forward-strand boundaries' $C/rev.err)" "1" "the error says only forward boundaries work"

# The left chunk shared by three --seam trim error cases below.
printf "$RS"'S\t1\tAAAAAAAAAA\nP\tGRCh38#0#chr1[0-10]\t1+\t*\n' > $C/bigl.gfa
vg convert -p $C/bigl.gfa > $C/bigl.pg 2>/dev/null

# Error: --seam trim requires the left and right chunks to share the same sequence.
printf "$RS"'S\t2\tAACCC\nS\t3\tGGGGG\nL\t2\t+\t3\t+\t0M\nP\tGRCh38#0#chr1[6-16]\t2+,3+\t*\n' > $C/mmr.gfa
vg convert -p $C/mmr.gfa > $C/mmr.pg 2>/dev/null
vg combine -s trim $C/bigl.pg $C/mmr.pg 2>$C/mm.err >/dev/null
is "$?" "1" "chunks that disagree on the overlapping sequence are rejected"
is "$(grep -c 'but disagree there' $C/mm.err)" "1" "the error says the overlap sequence does not match"
is "$(grep -c 'at offset 8' $C/mm.err)" "1" "the error reports the first position they differ at"

# Error: --seam trim requires that all paths start from the reference path.
printf "$RS"'S\t5\tGGTT\nS\t6\tA\nS\t7\tC\nS\t8\tGGGGG\nL\t5\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t8\t+\t0M\nL\t7\t+\t8\t+\t0M\nP\tGRCh38#0#chr1[8-18]\t5+,6+,8+\t*\nP\tHG002#1#chr1#0[8-18]\t5+,7+,8+\t*\nP\tHG002#2#chr1#0[12-18]\t6+,8+\t*\n' > $C/trb.gfa
printf "$RS"'S\t1\tAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tGGGGG\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+,4+\t*\nP\tHG002#1#chr1#0[0-10]\t1+,2+,4+\t*\nP\tHG002#2#chr1#0[0-5]\t1+,3+\t*\n' > $C/tlb.gfa
vg convert -p $C/trb.gfa > $C/trb.pg 2>/dev/null
vg convert -p $C/tlb.gfa > $C/tlb.pg 2>/dev/null

vg combine -s trim $C/tl.pg $C/trb.pg >/dev/null 2>$C/start.err
is "$?" "1" "a path that does not start at the REFERENCE start node is rejected"
is "$(grep -c 'does not start at the REFERENCE start node (id 5)' $C/start.err)" "1" "the error reports the path and the REFERENCE start node it should have started on"

# Error: --seam trim requires that all paths end at the reference path.
vg combine -s trim $C/tlb.pg $C/tr.pg >/dev/null 2>$C/end.err
is "$?" "1" "a path that does not end at the REFERENCE end node is rejected"
is "$(grep -c 'does not end at the REFERENCE end node (id 4)' $C/end.err)" "1" "the error reports the path and the REFERENCE end node it should have ended on"

# Error: --seam trim does not allow the overlap to exceed the length of the first node of the right chunk.
printf "$RS"'S\t2\tAA\nS\t3\tAAAAGG\nL\t2\t+\t3\t+\t0M\nP\tGRCh38#0#chr1[4-12]\t2+,3+\t*\n' > $C/bigr.gfa
vg convert -p $C/bigr.gfa > $C/bigr.pg 2>/dev/null
vg combine -s trim $C/bigl.pg $C/bigr.pg 2>$C/big.err >/dev/null
is "$?" "1" "an overlap longer than the first node is rejected"
is "$(grep -c 'exceeds the length of its first node' $C/big.err)" "1" "the error says the overlap does not fit"

# Error: --seam trim does not allow the overlap to equal the length of the first node of the right chunk.
printf "$RS"'S\t2\tAAA\nS\t3\tGGGGGG\nL\t2\t+\t3\t+\t0M\nP\tGRCh38#0#chr1[7-16]\t2+,3+\t*\n' > $C/eqr.gfa
vg convert -p $C/eqr.gfa > $C/eqr.pg 2>/dev/null
vg combine -s trim $C/bigl.pg $C/eqr.pg 2>$C/eq.err >/dev/null
is "$?" "1" "an overlap that would empty the first node is rejected"
is "$(grep -c 'would leave a zero-length' $C/eq.err)" "1" "the error says the node would be left empty"

# Error: --seam trim needs a clean chunk boundary.
printf "$RS"'S\t1\tAAAAA\nS\t2\tCCCCC\nL\t1\t+\t2\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+\t*\n' > $C/bl.gfa
printf "$RS"'S\t3\tCCGGG\nS\t4\tTTTTT\nL\t3\t+\t4\t+\t0M\nP\tGRCh38#0#chr1[8-18]\t3+,4+\t*\n' > $C/br.gfa
# Same right chunk, plus an edge back into the node the trim want to remove.
printf "$RS"'S\t3\tCCGGG\nS\t4\tTTTTT\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t3\t+\t0M\nP\tGRCh38#0#chr1[8-18]\t3+,4+\t*\n' > $C/br_back.gfa
# Same left chunk, plus a node hanging off the right of its boundary node.
printf "$RS"'S\t1\tAAAAA\nS\t2\tCCCCC\nS\t5\tAC\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t5\t+\t0M\nP\tGRCh38#0#chr1[0-10]\t1+,2+\t*\n' > $C/bl_dangle.gfa
for g in bl br br_back bl_dangle; do vg convert -p $C/$g.gfa > $C/$g.pg 2>/dev/null; done

vg combine -s trim $C/bl.pg $C/br.pg > $C/o_b.pg 2>/dev/null
is "$(vg paths -Ev $C/o_b.pg)" "$(printf 'GRCh38#0#chr1[0-18]\t18')" "the two chunks are joined into one path"
is "$(vg paths -Fv $C/o_b.pg | tail -1)" "AAAAACCCCCGGGTTTTT" "the overlapped CC is trimmed away and the two chunks are joined into one sequence"

vg combine -s trim $C/bl.pg $C/br_back.pg >/dev/null 2>$C/back.err
is "$?" "1" "an edge into the node that would be trimmed is rejected"
is "$(grep -c 'would be lost' $C/back.err)" "1" "the error says the edge would be lost"

# An edge off the left boundary is fine when seaming, but not for fusing.
vg combine -s trim $C/bl_dangle.pg $C/br.pg > $C/o_bd.pg 2>/dev/null
is "$(vg paths -Ev $C/o_bd.pg)" "$(printf 'GRCh38#0#chr1[0-18]\t18')" "an edge off the left boundary is fine for seaming"

vg combine -s trim -u $C/bl_dangle.pg $C/br.pg >/dev/null 2>$C/fuse.err
is "$?" "1" "--fuse rejects a boundary node with an edge on the welded side"
is "$(grep -c 'weld would have to drop' $C/fuse.err)" "1" "the error says what the weld would drop"

vg combine -s trim -u $C/bl.pg $C/br.pg > $C/o_bu.pg 2>/dev/null
is "$(vg stats -N $C/o_bu.pg)" "3" "a clean boundary is required for fusing into one node"

# Error: --seam shared requires that all shared nodes have the same sequence.
printf "$RS"'S\t1\tAAAA\nP\tGRCh38#0#chr1[0-4]\t1+\t*\n' > $C/c1.gfa
printf "$RS"'S\t1\tCCCC\nP\tGRCh38#0#chr1[4-8]\t1+\t*\n' > $C/c2.gfa
vg convert -p $C/c1.gfa > $C/c1.pg 2>/dev/null
vg convert -p $C/c2.gfa > $C/c2.pg 2>/dev/null
vg combine -s shared $C/c1.pg $C/c2.pg 2>$C/conflict.err >/dev/null
is "$?" "1" "--seam shared rejects same ID carrying different sequences"
is "$(grep -c 'has different sequences' $C/conflict.err)" "1" "the error reports the node that disagrees"

# Error: pieces of one path are told apart by name, so the same name cannot arrive twice.
printf "$RS"'S\t2\tGGGG\nP\tGRCh38#0#chr1[0-4]\t2+\t*\n' > $C/dup.gfa
vg convert -p $C/dup.gfa > $C/dup.pg 2>/dev/null
vg combine $C/c1.pg $C/dup.pg 2>$C/dup.err >/dev/null
is "$?" "1" "a path with the same name in both chunks is rejected"
is "$(grep -c 'in an earlier input' $C/dup.err)" "1" "the error reports the path that was duplicated"

# Phase-block fragments: blocks 0 and 10 read as offsets, block 1 does not.
printf 'H\tVN:Z:1.0\nS\t1\tACGTA\nS\t2\tCCGGT\nP\tHG002#1#x#0\t1+,2+\t*\nL\t1\t+\t2\t+\t0M\n' > $C/pb1.gfa
printf 'H\tVN:Z:1.0\nS\t3\tTTACG\nS\t4\tGGCAT\nP\tHG002#1#x#10\t3+,4+\t*\nL\t3\t+\t4\t+\t0M\n' > $C/pb2.gfa
printf 'H\tVN:Z:1.0\nS\t3\tTTACG\nS\t4\tGGCAT\nP\tHG002#1#x#1\t3+,4+\t*\nL\t3\t+\t4\t+\t0M\n' > $C/pbnum.gfa
vg convert -p $C/pb1.gfa > $C/pb1.pg 2>/dev/null
vg convert -p $C/pb2.gfa > $C/pb2.pg 2>/dev/null
vg convert -p $C/pbnum.gfa > $C/pbnum.pg 2>/dev/null

# By default a phase block is just a phase block, not an offset, so these fragments do not merge.
vg combine $C/pb1.pg $C/pb2.pg > $C/o_pb.vg 2>/dev/null
is "$(vg paths -Ev $C/o_pb.vg | wc -l)" "2" "phase blocks is not merged by default"
is "$(vg paths -Ev $C/o_pb.vg | cut -f1 | sort | paste -sd,)" "HG002#1#x#0,HG002#1#x#10" "the untouched fragments keep their names"

# -P reads the phase block as a start offset instead.
vg combine -P $C/pb1.pg $C/pb2.pg > $C/o_pbP.vg 2>/dev/null
is "$(vg paths -Ev $C/o_pbP.vg | cut -f1)" "HG002#1#x#0[0-20]" "-P merges by phase block into one spanning path"
is "$(vg paths -Ev $C/o_pbP.vg | cut -f2)" "20" "the merged length matches the span the blocks imply"

# The phase block as a start offset is not continuous.
vg combine -P $C/pb1.pg $C/pbnum.pg > $C/o_pbnum.vg 2>$C/pb.err
is "$?" "1" "-P reads the phase block as a start offset, but they are not continuous"
is "$(grep -c 'rather than offsets' $C/pb.err)" "1" "the error asks whether these are numbers, not offsets"

# A phase block that is not 0 has to survive the merge intact.
printf "$RS"'S\t1\tACGTA\nS\t2\tCCGGT\nL\t1\t+\t2\t+\t0M\nP\tHG002#2#chr22#4[0-10]\t1+,2+\t*\n' > $C/hb1.gfa
printf "$RS"'S\t3\tTTACG\nS\t4\tGGCAT\nL\t3\t+\t4\t+\t0M\nP\tHG002#2#chr22#4[10-20]\t3+,4+\t*\n' > $C/hb2.gfa
vg convert -p $C/hb1.gfa > $C/hb1.pg 2>/dev/null
vg convert -p $C/hb2.gfa > $C/hb2.pg 2>/dev/null

is "$(vg paths -M -v $C/hb1.pg | tail -1 | cut -f2,6)" "$(printf 'HAPLOTYPE\t4')" "the fixture really is a HAPLOTYPE path with phase block 4"
vg combine $C/hb1.pg $C/hb2.pg > $C/o_hb.pg 2>/dev/null
is "$(vg paths -Ev $C/o_hb.pg)" "$(printf 'HG002#2#chr22#4[0-20]\t20')" "merging keeps the phase block and spans both fragments"

# Error: a merged name that is already taken is rejected.
printf "$RS"'S\t1\tAAAACCCC\nP\tGRCh38#0#chr1[0-8]\t1+\t*\n' > $C/whole.gfa
vg convert -v $C/whole.gfa > $C/whole_lin.vg 2>/dev/null
vg circularize -p 'GRCh38#0#chr1[0-8]' $C/whole_lin.vg > $C/whole.vg 2>/dev/null
printf "$RS"'S\t1\tAAAA\nP\tGRCh38#0#chr1[0-4]\t1+\t*\n' > $C/fa.gfa
printf "$RS"'S\t1\tCCCC\nP\tGRCh38#0#chr1[4-8]\t1+\t*\n' > $C/fb.gfa
vg convert -p $C/fa.gfa > $C/fa.pg 2>/dev/null
vg convert -p $C/fb.gfa > $C/fb.pg 2>/dev/null
vg combine $C/whole.vg $C/fa.pg $C/fb.pg 2>$C/collide.err >/dev/null
is "$?" "1" "a merged name that is already taken is rejected"
is "$(grep -c 'already an unrelated path' $C/collide.err)" "1" "the error names the path standing in the way"

# vg combine has to work on graphs vg itself produces, not just hand-built GFA.
sed 's/^>x/>SAMPLE#0#x/' small/x.fa > $C/pansn.fa
vg construct -r $C/pansn.fa > $C/base.vg
vg index -x $C/base.xg $C/base.vg
is "$(vg paths -M -v $C/base.vg | tail -1 | cut -f2)" "REFERENCE" "PanSN naming yields a REFERENCE-sense path"

# Chunking by node range gives two chunks that meet end to end and share no node IDs,
# which is what the default policy expects.
vg chunk -x $C/base.xg -r 1:16  -c 0 > $C/adj1.vg
vg chunk -x $C/base.xg -r 17:32 -c 0 > $C/adj2.vg
vg combine $C/adj1.vg $C/adj2.vg > $C/o_real.vg 2>/dev/null
is "$(vg paths -Ev $C/o_real.vg)" "$(printf 'SAMPLE#0#x[0-1001]\t1001')" "two node-range chunks merge into one contig-spanning path"
is "$(vg stats -N $C/o_real.vg)" "32" "the merged graph keeps every input node"
is "$(vg stats -E $C/o_real.vg)" "31" "merging adds exactly one inter-fragment edge"

# Whatever VPKG can load, combine can merge.
vg convert -p $C/adj1.vg > $C/adj1.pg
vg convert -a $C/adj2.vg > $C/adj2.hg
vg combine $C/adj1.pg $C/adj2.hg > $C/o_fmt.vg 2>/dev/null
is "$(vg paths -Ev $C/o_fmt.vg)" "$(vg paths -Ev $C/o_real.vg)" "PackedGraph and HashGraph inputs merge the same way"
is "$(vg stats -N $C/o_fmt.vg)" "32" "and keep the same nodes"

rm -rf $C

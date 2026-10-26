#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 45

# Scratch directory for this test's files; removed at the end.
SCRATCH=combine_scratch
mkdir -p $SCRATCH

# GFA fixtures are written with spaces between fields; tr turns them into the tabs GFA needs.

# Construct chunks holding a reference plus two haplotypes.
# Left chunk including HG002 and HG005#1
tr ' ' '\t' > $SCRATCH/ml.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 1 AAAA
S 2 C
S 3 G
S 4 GGGGG
L 1 + 2 + 0M
L 1 + 3 + 0M
L 2 + 4 + 0M
L 3 + 4 + 0M
P GRCh38#0#chr1[0-10] 1+,2+,4+ *
P HG002#1#chr1#0[0-10] 1+,2+,4+ *
P HG002#2#chr1#0[0-10] 1+,3+,4+ *
P HG005#1#chr1#0[0-10] 1+,3+,4+ *
EOF
# Right chunk including only HG002
tr ' ' '\t' > $SCRATCH/mr.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 5 TTTT
S 6 A
S 7 C
S 8 GGGGG
L 5 + 6 + 0M
L 5 + 7 + 0M
L 6 + 8 + 0M
L 7 + 8 + 0M
P GRCh38#0#chr1[10-20] 5+,6+,8+ *
P HG002#1#chr1#0[10-20] 5+,7+,8+ *
P HG002#2#chr1#0[10-20] 5+,6+,8+ *
EOF
vg convert -p $SCRATCH/ml.gfa > $SCRATCH/ml.pg 2>/dev/null
vg convert -p $SCRATCH/mr.gfa > $SCRATCH/mr.pg 2>/dev/null

is "$(vg paths -M -v $SCRATCH/ml.pg | tail -n +2 | cut -f2 | grep -c HAPLOTYPE)" "3" "the chunk carries three haplotype paths"
is "$(vg paths -M -v $SCRATCH/ml.pg | tail -n +2 | cut -f2 | grep -c REFERENCE)" "1" "the chunk carries one reference path"

vg combine $SCRATCH/ml.pg $SCRATCH/mr.pg > $SCRATCH/o_def.vg 2>/dev/null
is "$(vg paths -Ev $SCRATCH/o_def.vg | wc -l)" "4" "every path is considered"
is "$(vg paths -Ev $SCRATCH/o_def.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-20],HG002#1#chr1#0[0-20],HG002#2#chr1#0[0-20],HG005#1#chr1#0[0-10]" "each sample and haplotype merges within its own group"
is "$(vg paths -Ev $SCRATCH/o_def.vg | cut -f2 | sort -n | paste -sd,)" "10,20,20,20" "the merged paths span both chunks and the lone one does not"
is "$(vg stats -N $SCRATCH/o_def.vg)" "8" "merging keeps every node of both chunks, 4 for each"
is "$(vg stats -E $SCRATCH/o_def.vg)" "9" "the three paths cross the seam on one shared edge, so the edges count increase by 1(8 -> 9)"
is "$(vg stats -s $SCRATCH/o_def.vg | wc -l)" "1" "the two chunks end up in one component"

# renumber is the default, so naming it explicitly must change nothing.
vg combine -s renumber $SCRATCH/ml.pg $SCRATCH/mr.pg > $SCRATCH/o_ren.vg 2>/dev/null
is "$(vg paths -Ev $SCRATCH/o_ren.vg)" "$(vg paths -Ev $SCRATCH/o_def.vg)" "--seam renumber is the default"
is "$(vg paths -Fv $SCRATCH/o_ren.vg | sort)" "$(vg paths -Fv $SCRATCH/o_def.vg | sort)" "and spells out the same sequences"

# --seam shared: the two chunks sharing nodes 4 and 5.
tr ' ' '\t' > $SCRATCH/sl.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 1 AAAA
S 2 C
S 3 G
S 4 GGGGG
S 5 TTTT
L 1 + 2 + 0M
L 1 + 3 + 0M
L 2 + 4 + 0M
L 3 + 4 + 0M
L 4 + 5 + 0M
P GRCh38#0#chr1[0-14] 1+,2+,4+,5+ *
P HG002#1#chr1#0[0-14] 1+,2+,4+,5+ *
P HG002#2#chr1#0[0-14] 1+,3+,4+,5+ *
EOF
tr ' ' '\t' > $SCRATCH/sr.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 4 GGGGG
S 5 TTTT
S 6 A
S 7 C
S 8 GGGGG
L 4 + 5 + 0M
L 5 + 6 + 0M
L 5 + 7 + 0M
L 6 + 8 + 0M
L 7 + 8 + 0M
P GRCh38#0#chr1[5-20] 4+,5+,6+,8+ *
P HG002#1#chr1#0[5-20] 4+,5+,7+,8+ *
P HG002#2#chr1#0[5-20] 4+,5+,6+,8+ *
EOF
vg convert -p $SCRATCH/sl.gfa > $SCRATCH/sl.pg 2>/dev/null
vg convert -p $SCRATCH/sr.gfa > $SCRATCH/sr.pg 2>/dev/null

vg combine -s shared $SCRATCH/sl.pg $SCRATCH/sr.pg > $SCRATCH/o_sh.vg 2>/dev/null
is "$(vg paths -Ev $SCRATCH/o_sh.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-20],HG002#1#chr1#0[0-20],HG002#2#chr1#0[0-20]" "--seam shared merges all three paths over their shared nodes"
is "$(vg paths -Ev $SCRATCH/o_sh.vg | cut -f2 | sort -u)" "20" "no path counts the same nodes twice"
is "$(vg stats -N $SCRATCH/o_sh.vg)" "8" "the shared nodes are deduplicated rather than copied"
is "$(vg stats -E $SCRATCH/o_sh.vg)" "9" "the edges between them are also deduplicated"
is "$(vg paths -Fv $SCRATCH/o_sh.vg | grep -A1 '^>HG002#1#chr1#0\[0-20\]$' | tail -1)" "AAAACGGGGGTTTTCGGGGG" "HG002#1#chr1 keeps its own sequences over the shared nodes"
is "$(vg paths -Fv $SCRATCH/o_sh.vg | grep -A1 '^>HG002#2#chr1#0\[0-20\]$' | tail -1)" "AAAAGGGGGGTTTTAGGGGG" "HG002#2#chr1 keeps its own sequences over the shared nodes"

# --seam trim: on a 2bp overlap.
tr ' ' '\t' > $SCRATCH/tl.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 1 AAAA
S 2 C
S 3 G
S 4 GGGGG
L 1 + 2 + 0M
L 1 + 3 + 0M
L 2 + 4 + 0M
L 3 + 4 + 0M
P GRCh38#0#chr1[0-10] 1+,2+,4+ *
P HG002#1#chr1#0[0-10] 1+,2+,4+ *
P HG002#2#chr1#0[0-10] 1+,3+,4+ *
EOF
tr ' ' '\t' > $SCRATCH/tr.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 5 GGTT
S 6 A
S 7 C
S 8 GGGGG
L 5 + 6 + 0M
L 5 + 7 + 0M
L 6 + 8 + 0M
L 7 + 8 + 0M
P GRCh38#0#chr1[8-18] 5+,6+,8+ *
P HG002#1#chr1#0[8-18] 5+,7+,8+ *
P HG002#2#chr1#0[8-18] 5+,6+,8+ *
EOF
vg convert -p $SCRATCH/tl.gfa > $SCRATCH/tl.pg 2>/dev/null
vg convert -p $SCRATCH/tr.gfa > $SCRATCH/tr.pg 2>/dev/null

vg combine -s trim $SCRATCH/tl.pg $SCRATCH/tr.pg > $SCRATCH/o_trim.vg 2>/dev/null
is "$(vg paths -Ev $SCRATCH/o_trim.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-18],HG002#1#chr1#0[0-18],HG002#2#chr1#0[0-18]" "all three paths report the same trimmed span"
is "$(vg paths -Ev $SCRATCH/o_trim.vg | cut -f2 | sort -u)" "18" "the overlap of three paths has been trimmed"
is "$(vg stats -N $SCRATCH/o_trim.vg)" "8" "trimming splits the boundary node and drops the overlapping half, leaving the node counts unchanged"
is "$(vg stats -E $SCRATCH/o_trim.vg)" "9" "seaming create one edge at the boundary"
is "$(vg paths -Fv $SCRATCH/o_trim.vg | grep -A1 '^>GRCh38#0#chr1\[0-18\]$' | tail -1)" "AAAACGGGGGTTAGGGGG" "the overlap of the reference path is trimmed"
is "$(vg paths -Fv $SCRATCH/o_trim.vg | grep -A1 '^>HG002#1#chr1#0\[0-18\]$' | tail -1)" "AAAACGGGGGTTCGGGGG" "the overlap of the HG002#1#chr1 is trimmed"
is "$(vg paths -Fv $SCRATCH/o_trim.vg | grep -A1 '^>HG002#2#chr1#0\[0-18\]$' | tail -1)" "AAAAGGGGGGTTAGGGGG" "the overlap of the HG002#2#chr1 is trimmed"

# --seam trim --fuse: must rewrite the step of every path, replace the original nodes with the fused nodes.
vg combine -s trim -u $SCRATCH/tl.pg $SCRATCH/tr.pg > $SCRATCH/o_fuse.vg 2>/dev/null
is "$(vg stats -N $SCRATCH/o_fuse.vg)" "7" "the two boundary nodes become one"
is "$(vg stats -E $SCRATCH/o_fuse.vg)" "8" "the edge between two boundary nodes is removed"
is "$(vg paths -Ev $SCRATCH/o_fuse.vg | cut -f1 | sort | paste -sd,)" "GRCh38#0#chr1[0-18],HG002#1#chr1#0[0-18],HG002#2#chr1#0[0-18]" "fusing does not change any path's coordinate"
is "$(vg paths -Ev $SCRATCH/o_fuse.vg | cut -f2 | sort -u)" "18" "fusing doesn't change any path's length"
is "$(vg paths -Fv $SCRATCH/o_fuse.vg | sort)" "$(vg paths -Fv $SCRATCH/o_trim.vg | sort)" "fusing doesn't change any path's sequence"

# Another single chunk.
tr ' ' '\t' > $SCRATCH/one.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 1 AAAA
S 2 CCC
L 1 + 2 + 0M
P GRCh38#0#chr1[0-7] 1+,2+ *
EOF
vg convert -p $SCRATCH/one.gfa > $SCRATCH/one.pg 2>/dev/null

# Error: nothing to combine.
vg combine -s renumber 2>$SCRATCH/none.err >/dev/null
is "$?" "1" "combining nothing is rejected"
is "$(grep -c 'At least one input graph' $SCRATCH/none.err)" "1" "the error says at least one input graph is required"

# Error: --seam trim needs more than one chunk.
vg combine -s trim $SCRATCH/one.pg 2>$SCRATCH/lone.err >/dev/null
is "$?" "1" "--seam trim on a single input is rejected"
is "$(grep -c 'at least two input graphs' $SCRATCH/lone.err)" "1" "the error says trim needs at least two inputs"

# Error: --fuse requires --seam trim.
vg combine -u $SCRATCH/tl.pg $SCRATCH/tr.pg 2>$SCRATCH/fusearg.err >/dev/null
is "$?" "1" "--fuse without --seam trim is rejected"
is "$(grep -c 'requires --seam trim' $SCRATCH/fusearg.err)" "1" "the error says --fuse needs --seam trim"

# Error: --seam only takes the policies it knows.
vg combine -s bogus $SCRATCH/tl.pg $SCRATCH/tr.pg 2>$SCRATCH/seam.err >/dev/null
is "$?" "1" "an unrecognized --seam policy is rejected"
is "$(grep -c 'renumber, shared, trim' $SCRATCH/seam.err)" "1" "the error lists the policies that do exist"

# Error: a read-only graph format cannot be combined; VPKG reports this itself.
vg gbwt -g $SCRATCH/m.gbz --gbz-format -G graphs/gfa_with_reference.gfa 2>/dev/null
vg combine $SCRATCH/m.gbz $SCRATCH/one.pg 2>/dev/null >/dev/null
is "$?" "1" "a GBZ input is rejected rather than crashing"

# Error: --seam trim joins chunks end to end, so a gap between them is rejected.
# This also checks that an error from the combining code reaches the command line.
tr ' ' '\t' > $SCRATCH/gapl.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 1 AAAAAAAAAA
P GRCh38#0#chr1[0-10] 1+ *
EOF
tr ' ' '\t' > $SCRATCH/gapr.gfa <<'EOF'
H VN:Z:1.1 RS:Z:GRCh38
S 2 CCCCCCCCCC
P GRCh38#0#chr1[20-30] 2+ *
EOF
vg convert -p $SCRATCH/gapl.gfa > $SCRATCH/gapl.pg 2>/dev/null
vg convert -p $SCRATCH/gapr.gfa > $SCRATCH/gapr.pg 2>/dev/null
vg combine -s trim $SCRATCH/gapl.pg $SCRATCH/gapr.pg 2>$SCRATCH/tgap.err >/dev/null
is "$?" "1" "--seam trim rejects a gap between chunks"
is "$(grep -c 'offset gap of 10 bp' $SCRATCH/tgap.err)" "1" "--seam trim also reports the size of the gap"

# Two fragments of one haplotype whose phase blocks, 0 and 10, are their start offsets.
tr ' ' '\t' > $SCRATCH/pb1.gfa <<'EOF'
H VN:Z:1.0
S 1 ACGTA
S 2 CCGGT
P HG002#1#x#0 1+,2+ *
L 1 + 2 + 0M
EOF
tr ' ' '\t' > $SCRATCH/pb2.gfa <<'EOF'
H VN:Z:1.0
S 3 TTACG
S 4 GGCAT
P HG002#1#x#10 3+,4+ *
L 3 + 4 + 0M
EOF
vg convert -p $SCRATCH/pb1.gfa > $SCRATCH/pb1.pg 2>/dev/null
vg convert -p $SCRATCH/pb2.gfa > $SCRATCH/pb2.pg 2>/dev/null

# -P reads the phase block as a start offset.
vg combine -P $SCRATCH/pb1.pg $SCRATCH/pb2.pg > $SCRATCH/o_pbP.vg 2>/dev/null
is "$(vg paths -Ev $SCRATCH/o_pbP.vg | cut -f1)" "HG002#1#x#0[0-20]" "-P merges by phase block into one spanning path"
is "$(vg paths -Ev $SCRATCH/o_pbP.vg | cut -f2)" "20" "the merged length matches the span the blocks imply"

# Combine chunks made by vg construct and vg chunk, not just hand-written GFA.
sed 's/^>x/>SAMPLE#0#x/' small/x.fa > $SCRATCH/pansn.fa
vg construct -r $SCRATCH/pansn.fa > $SCRATCH/base.vg
vg index -x $SCRATCH/base.xg $SCRATCH/base.vg
is "$(vg paths -M -v $SCRATCH/base.vg | tail -1 | cut -f2)" "REFERENCE" "PanSN naming yields a REFERENCE-sense path"

# Chunking by node range gives two chunks that meet end to end and share no node IDs,
# which is what the default policy expects.
vg chunk -x $SCRATCH/base.xg -r 1:16  -c 0 > $SCRATCH/adj1.vg
vg chunk -x $SCRATCH/base.xg -r 17:32 -c 0 > $SCRATCH/adj2.vg
vg combine $SCRATCH/adj1.vg $SCRATCH/adj2.vg > $SCRATCH/o_real.vg 2>/dev/null
is "$(vg paths -Ev $SCRATCH/o_real.vg)" "$(printf 'SAMPLE#0#x[0-1001]\t1001')" "two node-range chunks merge into one contig-spanning path"
is "$(vg stats -N $SCRATCH/o_real.vg)" "32" "the merged graph keeps every input node"
is "$(vg stats -E $SCRATCH/o_real.vg)" "31" "merging adds exactly one inter-fragment edge"

rm -rf $SCRATCH

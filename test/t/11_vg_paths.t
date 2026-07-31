#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 123

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg
vg construct -r small/x.fa -v small/x.vcf.gz > x2.vg
vg index -x x.xg x.vg
vg gbwt -v small/x.vcf.gz -o x.gbwt -x x.vg

# List path/thread names from various input formats
is "$(vg paths --list -v x2.vg)" "x" "path listing works from vg"
is "$(vg paths --list -x x.xg)" "x" "path listing works from XG"
is "$(vg paths --list -x x.xg -G)" "x" "generic path listing works from XG"
is $(vg paths --list -g x.gbwt | wc -l) 2 "thread listing works from GBWT"
is $(vg paths --list -g x.gbwt -H | wc -l) 2 "haplotype thread listing works from GBWT"

# Select threads by name
is $(vg paths --list -Q "1#0#x#" -g x.gbwt | wc -l) 1 "thread selection by name prefix works correctly"
is $(vg paths --list -S 1 -g x.gbwt | wc -l) 2 "thread selection by sample name works correctly"
vg paths --list -S 2 -g x.gbwt > out.txt 2> err.txt
is $(cat out.txt | wc -l) 0 "no threads are reported for invalid samples"
is $(grep "no matching" err.txt | wc -l) 1 "warning provided when 0 threads are matched"

# Exclude sample from selection
is $(vg paths --list -g x.gbwt --exclude-sample 1 2>/dev/null | wc -l) 0 "excluding sample removes all matching threads from GBWT listing"
is $(vg paths --list -g x.gbwt --exclude-sample nonexistent | wc -l) 2 "excluding nonexistent sample has no effect on GBWT listing"
is $(vg paths --list -Q "1#0#x#" -g x.gbwt --exclude-sample 1 2>/dev/null | wc -l) 0 "exclude-sample combined with -Q prefix selection works"
is $(vg paths --list -x x.xg --exclude-sample nonexistent | wc -l) 1 "excluding nonexistent sample has no effect on graph path listing"

# Extract threads as alignments
is $(vg paths -x x.xg -g x.gbwt -X | vg view -a -  | wc -l) 2 "vg paths may be used to extract threads"

# Extract threads as GAF alignments
is $(vg paths -x x.xg -g x.gbwt -A | grep -v "^@" | wc -l) 2 "vg paths may be used to extract threads as GAF"

# Extract paths as fasta
vg paths -x x.xg -Q x -F > x_from_xg.fa
diff x_from_xg.fa small/x.fa
is $? 0 "Fasta extracted from xg is the same as the input fasta"
vg paths -v x.vg -Q x -F > x_from_vg.fa
diff x_from_vg.fa small/x.fa
is $? 0 "Fasta extracted from vg is the same as the input fasta"
is $(vg paths -x x.xg -g x.gbwt -F | wc -l) 28 "Fasta extracted from threads has correct number of lines"
vg paths --paths-by fakename -v x.vg -F > out.txt 2> err.txt
is $(cat out.txt | wc -l) 0 "no paths are reported for invalid path name"
is $(grep "no matching" err.txt | wc -l) 1 "warning provided when 0 paths are matched"

touch empty.fa
vg construct -r empty.fa > empty.vg
vg gbwt --index-paths -x empty.vg -o empty.gbwt

vg paths --list -g empty.gbwt 2> err.txt
is $? 1 "vg paths exits with error when no paths are found"
is $(grep "does not contain" err.txt | wc -l) 1 "useful error provided when no paths are found in gbwt"

vg paths --list -x empty.vg 2> err.txt
is $? 1 "vg paths exits with error when no paths are found"
is $(grep "does not contain" err.txt | wc -l) 1 "useful error provided when no paths are found in vg"

is $(vg paths -v msgas/s.vg -r -Q s1 | vg view - | grep ^P | cut -f 3 | sort | uniq | wc -l) 1 "a single path may be retained"

is $(vg paths -v msgas/s.vg -r -Q s1 | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) $(vg view msgas/s.vg | grep -v ^P | md5sum | cut -f 1 -d\ ) "path filtering does not modify the graph"

is $(vg construct -a -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg paths -d -a -v - | vg paths -L -v - | wc -l) 1 "alt allele paths can be dropped"

rm -f x.xg x.gbwt x.vg x2.vg x_from_xg.fa x_from_vg.fa

is $(vg paths -cv msgas/q.vg | awk '{print NF; exit}') 4 "vg path coverage has correct number of columns"
is $(vg paths -cv msgas/q.vg | wc -l) 4 "vg path coverage has correct number of rows"

# note: coverage doesn't include cycles at moment, so s2 path will not have full length
vg paths -Q s2 -v msgas/q.vg -d | vg paths -cv - | grep -v ^Path | awk '{print $1 "\t" $2}' > q.cov.len
vg paths -Q s2 -v msgas/q.vg -d | vg paths -Ev - > q.len
is $(cat q.len | wc -l) 2 "vg paths found correct number of lengths"
diff q.cov.len q.len
is $? 0 "vg path coverage reports correct lengths in first column"

rm -f q.cov.len q.len

# redundant paths are x2,x3,x4,x5 but not x1
vg paths -x graphs/path_norm_test.gfa -n -Q x2 > norm_x2.gfa
vg validate norm_x2.gfa
is $? 0 "path normalizer produces valid graph"

vg paths -x graphs/path_norm_test.gfa -F | sort > original.fa
vg paths -x norm_x2.gfa -F | sort > norm_x2.fa
diff original.fa norm_x2.fa
is $? 0 "path normalizer doesnt alter path sequences"

grep x2 norm_x2.gfa | awk '{print $3}' > x2.path
grep x2 norm_x2.gfa | awk '{print $3}' >> x2.path
grep x3 norm_x2.gfa | awk '{print $3}'> x2.norm.path
grep x5 norm_x2.gfa | awk '{print $3}' >> x2.norm.path
diff x2.path x2.norm.path
is $? 0 "path normalizer correctly snapped all equivalent paths to x2"

rm -f norm_x2.gfa original.fa norm_x2.fa x2.path x2.norm.path

vg paths -x graphs/path_norm_test.gfa -n -Q x4 > norm_x4.gfa
vg validate norm_x4.gfa
is $? 0 "path normalizer produces valid graph"

vg paths -x graphs/path_norm_test.gfa -F | sort > original.fa
vg paths -x norm_x4.gfa -F | sort > norm_x4.fa
diff original.fa norm_x4.fa
is $? 0 "path normalizer doesn't alter path sequences"

# note: x3 is x4 in reverse, so we key on that
grep x3 norm_x4.gfa | awk '{print $3}' > x4.path
grep x3 norm_x4.gfa | awk '{print $3}' >> x4.path
grep x3 norm_x4.gfa | awk '{print $3}'> x4.norm.path
grep x5 norm_x4.gfa | awk '{print $3}' >> x4.norm.path
diff x4.path x4.norm.path
is $? 0 "path normalizer correctly snapped all equivalent paths to x4"

vg convert --gfa-in graphs/named_with_walk.gfa > x.pg
vg gbwt --index-paths -x x.pg -o x.gbwt
vg gbwt --gbz-format -x x.pg x.gbwt -g x.gbz
vg paths --list -x x.gbz >out.txt

is "${?}" "0" "vg paths can list paths from a GBZ with only haplotypes"
is "$(cat out.txt | wc -l)" "1" "vg paths sees the haplotype path in a GBZ with only haplotypes"
diff <(vg paths -x tiny/tiny.gfa -F | sort) <(vg paths -x tiny/tiny.gfaz -F | sort)
is "${?}" "0" "vg paths emits matching path FASTA for equivalent GFA and GFAZ inputs"

rm -f norm_x4.gfa original.fa norm_x4.fa x4.path x4.norm.path out.txt err.txt
rm -f empty.vg empty.gbwt empty.fa

# Gref (graph reference path) computation tests
vg paths -x nesting/nested_snp_in_ins.gfa -Q x --compute-gref --min-gref-len 1 > gref_test.vg
vg validate gref_test.vg
is $? 0 "gref computation produces valid graph"

is $(vg paths -x gref_test.vg -L | grep "_alt$" | wc -l) 2 "gref computation creates expected number of gref paths"

is $(vg paths -x gref_test.vg -L | grep "^x$" | wc -l) 1 "original reference path is preserved after gref computation"

is $(vg paths -x gref_test.vg -L | grep -c "^gref_x$") 1 "reference path is copied into the gref namespace"

# Test gref naming convention matches pattern gref_x_{N}_alt
is $(vg paths -x gref_test.vg -L | grep -E "^gref_x_[0-9]+_alt$" | wc -l) 2 "gref paths follow naming convention gref_{path}_{N}_alt"

# Test with triple_nested.gfa which has more complex structure
vg paths -x nesting/triple_nested.gfa -Q x --compute-gref --min-gref-len 1 > triple_gref.vg
vg validate triple_gref.vg
is $? 0 "gref computation works on complex nested structure"

is $(vg paths -x triple_gref.vg -L | grep "_alt$" | wc -l) 2 "correct number of gref paths for triple nested graph"

# Test minimum length filter
vg paths -x nesting/triple_nested.gfa -Q x --compute-gref --min-gref-len 100 > triple_gref_long.vg
is $(vg paths -x triple_gref_long.vg -L | grep "_alt$" | wc -l) 0 "min-gref-len filters out short fragments"

# Test second pass coverage of dangling nodes (nodes outside snarls but on haplotype paths)
vg paths -x nesting/dangling_node.gfa -Q x --compute-gref --min-gref-len 1 > dangling_gref.vg
vg validate dangling_gref.vg
is $? 0 "gref computation handles dangling nodes outside snarls"

is $(vg paths -x dangling_gref.vg -L | grep "_alt$" | wc -l) 2 "gref second pass covers dangling nodes"

# Verify the dangling node (node 5, 8bp) is covered by checking gref path lengths include 8bp
# Note: use -E with grep instead of -Q since gref paths are filtered from prefix matching
is $(vg paths -x dangling_gref.vg -E | grep "_alt" | awk '{sum+=$2} END {print sum}') 12 "gref paths cover both snarl node and dangling node"

# Test --gref-segs option for writing segment table
vg paths -x nesting/nested_snp_in_ins.gfa -Q x --compute-gref --min-gref-len 1 --gref-segs gref_test.segs > gref_segs_test.vg
is $? 0 "gref-segs option produces no error"

is $(wc -l < gref_test.segs) 2 "gref-segs produces correct number of lines"

is $(cut -f4 gref_test.segs | grep -cE "^gref_x_[0-9]+_alt$") 2 "gref-segs contains gref path names"

is $(cut -f1 gref_test.segs | grep -c "#") 2 "gref-segs contains source path names with metadata"

is $(cut -f5 gref_test.segs | grep -c "^x$") 2 "gref-segs contains reference path name"

# Test that gref-segs requires compute-gref
vg paths -x nesting/nested_snp_in_ins.gfa -Q x -L --gref-segs gref_test.segs 2>&1 | grep -q "requires --compute-gref"
is $? 0 "gref-segs requires compute-gref option"

# The gref namespace is a convention, not an option: --compute-gref always writes it
vg paths -x nesting/nested_snp_in_ins.gfa -Q x --compute-gref --min-gref-len 1 --gref-segs gref_sample_test.segs > gref_sample_test.vg
is $(cut -f4 gref_sample_test.segs | grep -c "^gref_") 2 "gref-segs names are all in the gref namespace"

is $(vg paths -x gref_sample_test.vg -L | grep -c "^gref_") 3 "compute-gref writes the base copy and its fragments together"

# Test cross-path interval merging (left merge: new interval absorbs previous from different path)
vg paths -x nesting/cross_path_merge.gfa -Q x --compute-gref --min-gref-len 1 > cross_merge_test.vg
vg validate cross_merge_test.vg
is $? 0 "cross-path merge: gref computation produces valid graph"

# Cross-path merge should combine snarl interval [2,3,4] + dangling [9] into one path on hap3
# Without merging: 3 gref paths. With merging: 2 gref paths.
is $(vg paths -x cross_merge_test.vg -L | grep "_alt$" | wc -l) 2 "cross-path left merge reduces gref path count"

# Test cross-path interval merging (right merge: new interval absorbs following from different path)
vg paths -x nesting/cross_path_merge_right.gfa -Q x --compute-gref --min-gref-len 1 > cross_merge_right_test.vg
vg validate cross_merge_right_test.vg
is $? 0 "cross-path right merge: gref computation produces valid graph"

# Cross-path merge should combine dangling [9] + snarl interval [2,3,4] into one path on hap3
is $(vg paths -x cross_merge_right_test.vg -L | grep "_alt$" | wc -l) 2 "cross-path right merge reduces gref path count"

# A haplotype that walks the whole reference path with extra sequence before it must not
# be able to swallow the rank-0 reference interval during cross-path merging.  When it
# did, the reference slot was taken over by the haplotype: fragments got named after it
# (a#1#y0#0_1_alt) and the extra sequence lost its cover entirely.
vg paths -x nesting/hap_extends_ref_start.gfa -Q x --compute-gref --min-gref-len 1 --gref-segs ref_start.segs > ref_start_test.vg
vg validate ref_start_test.vg
is $? 0 "haplotype extending past reference start: gref computation produces valid graph"

is $(vg paths -x ref_start_test.vg -L | grep "_alt$" | wc -l) 2 "haplotype extending past reference start covers both off-reference nodes"

is $(vg paths -x ref_start_test.vg -L | grep -cE "^gref_x_[0-9]+_alt$") 2 "gref paths stay named after the reference, not the haplotype that spans it"

is "$(vg paths -x ref_start_test.vg -E | grep "_alt" | awk '{sum+=$2} END {print sum+0}')" "16" "gref cover includes the sequence before the reference start"

is $(cut -f5 ref_start.segs | grep -c "^x$") 2 "gref-segs reference column stays on the reference path"

# Same, with the extra sequence after the reference path's last node (the other merge branch)
vg paths -x nesting/hap_extends_ref_end.gfa -Q x --compute-gref --min-gref-len 1 > ref_end_test.vg
vg validate ref_end_test.vg
is $? 0 "haplotype extending past reference end: gref computation produces valid graph"

is $(vg paths -x ref_end_test.vg -L | grep "_alt$" | wc -l) 2 "haplotype extending past reference end covers both off-reference nodes"

is $(vg paths -x ref_end_test.vg -L | grep -cE "^gref_x_[0-9]+_alt$") 2 "gref paths after reference end stay named after the reference"

# 16bp = the two off-reference nodes only.  A gref path that had absorbed the reference
# interval would span it too and come to 24bp.
is "$(vg paths -x ref_end_test.vg -E | grep "_alt" | awk '{sum+=$2} END {print sum+0}')" "16" "gref paths do not absorb the reference path itself"

# Intervals that abut across an orientation flip must not be merged: the result would be a
# mixed-orientation interval, which apply() and write_gref_segments() both skip, silently
# dropping sequence that the cover still counts as covered.
vg paths -x nesting/orientation_flip.gfa -Q x --compute-gref --min-gref-len 1 --gref-segs flip.segs > flip_test.vg
vg validate flip_test.vg
is $? 0 "orientation flip: gref computation produces valid graph"

is $(vg paths -x flip_test.vg -L | grep "_alt$" | wc -l) 2 "intervals on either side of an orientation flip are kept separate"

is "$(vg paths -x flip_test.vg -E | grep "_alt" | awk '{sum+=$2} END {print sum+0}')" "32" "no sequence is dropped at an orientation flip"

is $(wc -l < flip.segs) 2 "gref-segs describes every emitted fragment at an orientation flip"

# Fragments in a component with no reference path get named after their source path.  The
# name still has to be a valid path name: the "_{N}_alt" suffix must land on the locus, not
# after a phase block, or the result parses as GENERIC with a '#' inside the locus and drops
# out of the gref sample entirely.
vg paths -x nesting/unanchored_component.gfa -Q GRCh38 --compute-gref --min-gref-len 1 --gref-segs unanchored.segs > unanchored_test.vg
vg validate unanchored_test.vg
is $? 0 "reference-disconnected component: gref computation produces valid graph"

is $(vg paths -x unanchored_test.vg -M | grep "_alt" | cut -f2 | sort -u) "REFERENCE" "gref paths in a reference-disconnected component are still reference sense"

is $(vg paths -x unanchored_test.vg -M | grep "_alt" | cut -f5 | grep -c "#") 0 "gref path names never leave a separator inside the locus"

is $(vg paths -x unanchored_test.vg -S gref_GRCh38 -L | wc -l) 2 "the anchored fragment and its base copy share the reference's gref sample"

is $(vg paths -x unanchored_test.vg -L | grep -cE "^gref_HG[12]#1#ctgZ_[0-9]+_alt$") 2 "fragments with no reference to reach are namespaced under the path they came from"

diff <(cut -f4 unanchored.segs | sort) <(vg paths -x unanchored_test.vg -L | grep "_alt" | sort)
is $? 0 "gref-segs names match the gref paths that were created"

# A PanSN reference read from a GFA without an RS header comes in as haplotype sense, so its
# name carries a phase block (GRCh38#0#chr1#0) that must not survive into the gref name.
vg paths -x nesting/haplotype_sense_ref.gfa -Q GRCh38 --compute-gref --min-gref-len 1 > hap_sense_test.vg
is $(vg paths -x hap_sense_test.vg -L | grep -c "^gref_GRCh38#0#chr1_[0-9]*_alt$") 2 "gref names off a haplotype-sense reference drop the phase block"

is $(vg paths -x hap_sense_test.vg -M | grep "_alt" | cut -f2 | sort -u) "REFERENCE" "gref paths off a haplotype-sense reference are reference sense"

# Subpaths of one contig must stay distinct in the gref namespace.  The fragment base
# name drops the subrange (the "_{N}_alt" suffix has to land on the locus), but the copy
# keeps it: collapsing the copies would publish only one subpath and silently drop the
# other one's sequence.
vg paths -x nesting/subranged_ref.gfa -Q GRCh38 --compute-gref --min-gref-len 1 > subranged_test.vg
vg validate subranged_test.vg
is $? 0 "subranged reference: gref computation produces valid graph"

is $(vg paths -x subranged_test.vg -L | grep -cE "^gref_GRCh38#0#chr1\[[0-9]+-[0-9]+\]$") 2 "each reference subpath gets its own gref copy"

is $(vg paths -x subranged_test.vg -L | grep -cE "^gref_GRCh38#0#chr1_[0-9]+_alt$") 2 "fragments off different subpaths get distinct names from the shared counter"

rm -f subranged_test.vg

# An inverted allele that another haplotype walks partly forward must still be covered in
# one piece.  Preferring a forward fragment over a longer reverse one shatters it into two
# sub-minimum pieces, and --min-gref-len then deletes both, so the sequence disappears from
# the cover entirely -- and adding a haplotype to a graph would remove sequence from it.
vg paths -x nesting/inverted_allele.gfa -Q GRCh38 --compute-gref --min-gref-len 50 > inverted_test.vg
is $(vg paths -x inverted_test.vg -L | grep -c "_alt$") 1 "an inverted allele survives the length filter as one fragment"

is "$(vg paths -x inverted_test.vg -E | grep "_alt" | cut -f2)" "60" "the whole inverted allele is covered, not just the part one haplotype walks forward"

# Same graph with no length filter: still one piece, and every gref path is nodes-forward
vg paths -x nesting/inverted_allele.gfa -Q GRCh38 --compute-gref --min-gref-len 1 > inverted_l1.vg
is $(vg convert -f inverted_l1.vg 2>/dev/null | grep -E "^[PW]" | grep "gref_" | grep -c "<") 0 "gref paths never walk a node backwards"

rm -f inverted_test.vg inverted_l1.vg

# --- Characterization tests -------------------------------------------------
# These pin gref cover behaviour that previously existed only in code, ahead of
# the interval-merging refactor (plan-gref-cover-refactor.md).  They are written
# against the current implementation and must keep passing through it.

# A cyclic reference path is a hard error, not a corrupt cover.  The cover
# requires disjoint acyclic reference paths, and a path visiting a node twice
# would give one node two rank-0 owners.
vg paths -x nesting/cyclic_ref_multiple_variants.gfa -Q x --compute-gref --min-gref-len 1 > cyclic_ref.vg 2> cyclic_ref.err
is $? 1 "a cyclic reference path is rejected rather than producing a bad cover"

is $(grep -c "disjoint acyclic reference paths" cyclic_ref.err) 1 "the cyclic reference error explains the requirement"

rm -f cyclic_ref.vg cyclic_ref.err

# Recomputing the cover on a graph that already has one must reproduce it
# exactly: clear() drops every gref path first, and gref paths are never used as
# fragment sources.  Idempotence is what makes the cover safe to recompute.
vg paths -x nesting/triple_nested.gfa -Q x --compute-gref --min-gref-len 1 > idem_1.vg
vg paths -x idem_1.vg -Q x --compute-gref --min-gref-len 1 > idem_2.vg
is $? 0 "recomputing a gref cover over an existing one succeeds"

is $(vg paths -x idem_2.vg -L | grep -c "^gref_") $(vg paths -x idem_1.vg -L | grep -c "^gref_") "recomputing the cover does not accumulate gref paths"

diff <(vg paths -x idem_1.vg -L | grep "^gref_" | sort) <(vg paths -x idem_2.vg -L | grep "^gref_" | sort)
is $? 0 "recomputing the cover reproduces the same gref path names"

diff <(vg paths -x idem_1.vg -E | grep "_alt" | cut -f2 | sort) <(vg paths -x idem_2.vg -E | grep "_alt" | cut -f2 | sort)
is $? 0 "recomputing the cover reproduces the same fragment lengths"

rm -f idem_1.vg idem_2.vg

# A node on no path at all cannot be covered -- the cover only ever emits
# substrings of existing paths.  It must warn and carry on, not crash or drop
# the rest of the cover.
vg paths -x nesting/unpathed_node.gfa -Q x --compute-gref --min-gref-len 1 > unpathed.vg 2> unpathed.err
is $? 0 "a node on no path does not stop the cover"

vg validate unpathed.vg
is $? 0 "a graph with an unpathed node still produces a valid cover"

is $(grep -c "not covered by gref paths" unpathed.err) 1 "an unpathed node is reported as uncovered"

is $(vg paths -x unpathed.vg -L | grep -c "_alt$") 1 "the pathed off-reference node is still covered"

rm -f unpathed.vg unpathed.err

# Cross-path merging must be refused when the two paths diverge inside the
# stretch being merged.  cross_path_merge.gfa covers the accept case; here
# a#3#y2 reaches node 9 through 5,3 while a#1#y0 reaches 3 through 2, so no fragment may
# span both routes: that would emit a path which is not a substring of any haplotype.
# Selecting longest-first takes a#3#y2's whole 5,3,9 run (16 bp) and then a#1#y0's leftover
# node 2, so the four off-reference nodes are covered in two real walks.  Name-ordered
# greedy selection used to take a#1#y0's 2,3 first and strand 5 and 9 as singletons, giving
# three.  Either way nothing is spliced; the count is what changed.
vg paths -x nesting/cross_path_merge_reject.gfa -Q x --compute-gref --min-gref-len 1 > reject_test.vg
vg validate reject_test.vg
is $? 0 "cross-path merge reject: gref computation produces valid graph"

is $(vg paths -x reject_test.vg -L | grep -c "_alt$") 2 "diverging paths are covered by whole runs, not spliced together"
# The real invariant: every fragment is a contiguous walk of one source path.  5,3,9 is
# a#3#y2 exactly; a fragment containing both 2 and 5 would be the splice this guards against.
is $(vg convert -f reject_test.vg 2>/dev/null | grep -E "^[PW]" | grep "_alt" | grep -c "2+.*5+\|5+.*2+") 0 "no fragment splices the two divergent routes together"

is "$(vg paths -x reject_test.vg -E | grep "_alt" | awk '{sum+=$2} END {print sum+0}')" "20" "refusing the merge still covers every off-reference node"

rm -f reject_test.vg

# A fragment that coalesces candidates from more than one snarl still gets one
# segment line, with reference coordinates spanning the enclosing snarl.
vg paths -x nesting/consecutive_nested.gfa -Q x --compute-gref --min-gref-len 1 --gref-segs consec.segs > consec.vg
is $(vg paths -x consec.vg -L | grep -c "_alt$") 1 "consecutive nested snarls are covered by a single fragment"

is $(wc -l < consec.segs) 1 "a fragment spanning two nested snarls gets one segment line"

is "$(cut -f5,6,7 consec.segs)" "$(printf 'x\t0\t2')" "the segment reference interval spans the enclosing snarl"

rm -f consec.vg consec.segs

# --compute-gref writes paths into the graph, so a read-only input has to say so
# and say what to do about it, rather than reporting a generic failure.
vg convert -p nesting/nested_snp_in_ins.gfa > gref_mut.pg
vg gbwt -G nesting/nested_snp_in_ins.gfa --gbz-format -g gref_ro.gbz 2>/dev/null
vg paths -x gref_ro.gbz -Q x --compute-gref --min-gref-len 1 > /dev/null 2> ro.err
is $? 1 "--compute-gref refuses a read-only graph format"

is $(grep -c "read-only format" ro.err) 1 "the read-only error names the problem"

is $(grep -c "vg convert -p" ro.err) 1 "the read-only error says how to fix it"

vg paths -x gref_mut.pg -o -Q x --compute-gref --min-gref-len 1 > /dev/null 2> ov.err
is $? 1 "--compute-gref refuses -o/--overlay rather than blaming the input"

is $(grep -c "overlay" ov.err) 1 "the overlay error names the overlay"

rm -f gref_mut.pg gref_ro.gbz ro.err ov.err

# --- Thread determinism -----------------------------------------------------
# The cover is computed by an OMP task loop over top-level snarls
# (SnarlManager::for_each_top_level_snarl_parallel), so which thread sees which
# snarl is picked by the runtime task scheduler and varies between thread counts
# *and* between runs at a fixed thread count.  The per-thread interval lists
# therefore reach the fold in a different order every time, and the sort at the
# top of the fold in GrefCover::compute() is the only thing that makes the
# result canonical.  apply() then numbers the fragments in fold order, so a
# broken sort renames every fragment -- it does not fail quietly.
#
# thread_determinism.gfa exists because none of the fixtures above can catch
# that.  They have one or two top-level snarls, so there is nothing to spread
# over threads; a graph built from small/x.vcf.gz has 70, but each is an
# isolated bubble contributing one interval, so its cover comes out the same
# whatever order the fold sees.  This graph chains 32 star clusters (three
# competing traversals each, every other one inverted) along a reference that
# skips all of them, giving 32 independently scheduled snarls whose intervals
# do compete.  Neutering the sort comparator to a constant leaves the graph from
# x.vcf.gz byte-identical at every thread count; this graph then comes out
# differently at nearly every thread count, and differently between runs at the
# same one.
vg paths -x nesting/thread_determinism.gfa -Q x --compute-gref --min-gref-len 1 -t 1 --gref-segs det_t1.segs > det_t1.vg 2>/dev/null
vg validate det_t1.vg
is $? 0 "the thread determinism graph produces a valid cover"

is $(vg paths -x det_t1.vg -L | grep -c "_alt$") 112 "the thread determinism graph covers enough snarls to schedule across threads"

# The same binary on the same input at a different thread count has to write the
# same bytes, so compare the serialized graphs directly.  That covers fragment
# contents, fragment names, and the order they were added to the graph in one
# shot, and it stays honest under refactoring: a change to the cover shifts every
# thread count together, only a thread-order dependence splits them apart.
for t in 2 4 8; do
    vg paths -x nesting/thread_determinism.gfa -Q x --compute-gref --min-gref-len 1 -t $t --gref-segs det_t$t.segs > det_t$t.vg 2>/dev/null

    cmp -s det_t1.vg det_t$t.vg
    is $? 0 "the gref cover at -t $t is identical to the one at -t 1"

    cmp -s det_t1.segs det_t$t.segs
    is $? 0 "the gref segments at -t $t are identical to those at -t 1"
done

# If the byte comparison above ever fails, this one says which fragment moved
# instead of just that some byte did.
diff <(vg convert -f det_t1.vg 2>/dev/null | awk '$1=="P"||$1=="W"' | grep gref_ | sort) \
     <(vg convert -f det_t8.vg 2>/dev/null | awk '$1=="P"||$1=="W"' | grep gref_ | sort)
is $? 0 "every gref path at -t 8 walks the same nodes as its -t 1 counterpart"

# Same thread count, run again.  OMP hands the snarl tasks out afresh on every
# run, so this catches order dependence that a fixed thread count would hide.
vg paths -x nesting/thread_determinism.gfa -Q x --compute-gref --min-gref-len 1 -t 8 --gref-segs det_rep.segs > det_rep.vg 2>/dev/null

cmp -s det_t8.vg det_rep.vg
is $? 0 "repeating the same -t 8 run reproduces the cover exactly"

cmp -s det_t8.segs det_rep.segs
is $? 0 "repeating the same -t 8 run reproduces the same segments"

rm -f det_t1.vg det_t2.vg det_t4.vg det_t8.vg det_rep.vg
rm -f det_t1.segs det_t2.segs det_t4.segs det_t8.segs det_rep.segs

rm -f gref_test.vg triple_gref.vg triple_gref_long.vg dangling_gref.vg x.pg x.gbwt x.gbz
rm -f gref_test.segs gref_segs_test.vg gref_sample_test.segs gref_sample_test.vg
rm -f cross_merge_test.vg cross_merge_right_test.vg
rm -f ref_start_test.vg ref_start.segs ref_end_test.vg flip_test.vg flip.segs
rm -f unanchored_test.vg unanchored.segs hap_sense_test.vg




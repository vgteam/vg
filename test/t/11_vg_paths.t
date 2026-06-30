#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 60

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

rm -f norm_x4.gfa original.fa norm_x4.fa x4.path x4.norm.path out.txt err.txt
rm -f empty.vg empty.gbwt empty.fa

# Gref (graph reference path) computation tests
vg paths -x nesting/nested_snp_in_ins.gfa -Q x --compute-gref --min-gref-len 1 > gref_test.vg
vg validate gref_test.vg
is $? 0 "gref computation produces valid graph"

is $(vg paths -x gref_test.vg -L | grep "_alt$" | wc -l) 2 "gref computation creates expected number of gref paths"

is $(vg paths -x gref_test.vg -L | grep "^x$" | wc -l) 1 "original reference path is preserved after gref computation"

# Test gref naming convention matches pattern x_{N}_alt
is $(vg paths -x gref_test.vg -L | grep -E "^x_[0-9]+_alt$" | wc -l) 2 "gref paths follow naming convention path_{N}_alt"

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

is $(cut -f4 gref_test.segs | grep -c "x_.*_alt") 2 "gref-segs contains gref path names"

is $(cut -f1 gref_test.segs | grep -c "#") 2 "gref-segs contains source path names with metadata"

is $(cut -f5 gref_test.segs | grep -c "^x$") 2 "gref-segs contains reference path name"

# Test that gref-segs requires compute-gref
vg paths -x nesting/nested_snp_in_ins.gfa -Q x -L --gref-segs gref_test.segs 2>&1 | grep -q "requires --compute-gref"
is $? 0 "gref-segs requires compute-gref option"

# Test gref-segs with gref-sample option
vg paths -x nesting/nested_snp_in_ins.gfa -Q x --compute-gref --min-gref-len 1 --gref-sample TESTSAMPLE --gref-segs gref_sample_test.segs > gref_sample_test.vg
is $(cut -f4 gref_sample_test.segs | grep -c "TESTSAMPLE") 2 "gref-segs uses gref-sample for path names"

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

rm -f gref_test.vg triple_gref.vg triple_gref_long.vg dangling_gref.vg x.pg x.gbwt x.gbz
rm -f gref_test.segs gref_segs_test.vg gref_sample_test.segs gref_sample_test.vg
rm -f cross_merge_test.vg cross_merge_right_test.vg





#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 26

vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg
vg construct -r small/x.fa -v small/x.vcf.gz > x2.vg
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg

# List path/thread names from various input formats
is "$(vg paths --list -v x2.vg)" "x" "path listing works from vg"
is "$(vg paths --list -x x.xg)" "x" "path listing works from XG"
is "$(vg paths --list -x x.xg -G)" "x" "generic path listing works from XG"
is $(vg paths --list -g x.gbwt | wc -l) 2 "thread listing works from GBWT"
is $(vg paths --list -g x.gbwt -H | wc -l) 2 "haplotype thread listing works from GBWT"

# Select threads by name
is $(vg paths --list -Q "1#0#x#" -g x.gbwt | wc -l) 1 "thread selection by name prefix works correctly"
is $(vg paths --list -S 1 -g x.gbwt | wc -l) 2 "thread selection by sample name works correctly"
is $(vg paths --list -S 2 -g x.gbwt | wc -l) 0 "no threads are reported for invalid samples"

# Extract threads as alignments
is $(vg paths -x x.xg -g x.gbwt -X | vg view -a -  | wc -l) 2 "vg paths may be used to extract threads"

# Extract threads as GAF alignments
is $(vg paths -x x.xg -g x.gbwt -A | wc -l) 2 "vg paths may be used to extract threads as GAF"

# Extract paths as fasta
vg paths -x x.xg -Q x -F > x_from_xg.fa
diff x_from_xg.fa small/x.fa
is $? 0 "Fasta extracted from xg is the same as the input fasta"
vg paths -v x.vg -Q x -F > x_from_vg.fa
diff x_from_vg.fa small/x.fa
is $? 0 "Fasta extracted from vg is the same as the input fasta"
is $(vg paths -x x.xg -g x.gbwt -F | wc -l) 28 "Fasta extracted from threads has correct number of lines"

is $(vg msga -w 20 -f msgas/s.fa | vg paths -v - -r -Q s1 | vg view - | grep ^P | cut -f 3 | sort | uniq | wc -l) 1 "a single path may be retained"

is $(vg msga -w 20 -f msgas/s.fa  | vg paths -v - -r -Q s1 | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) $(vg msga -w 20 -f msgas/s.fa  | vg view - | grep -v ^P | md5sum | cut -f 1 -d\ ) "path filtering does not modify the graph"

is $(vg construct -a -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg paths -d -a -v - | vg paths -L -v - | wc -l) 1 "alt allele paths can be dropped"

rm -f x.xg x.gbwt x.vg x2.vg x_from_xg.fa x_from_vg.fa q.vg

vg msga -w 20 -f msgas/q.fa  > q.vg
is $(vg paths -cv q.vg | awk '{print NF; exit}') 4 "vg path coverage has correct number of columns"
is $(vg paths -cv q.vg | wc -l) 4 "vg path coverage has correct number of rows"

# note: coverage doesn't include cycles at moment, so s2 path will not have full length
vg paths -Q s2 -v q.vg -d | vg paths -cv - | grep -v ^Path | awk '{print $1 "\t" $2}' > q.cov.len
vg paths -Q s2 -v q.vg -d | vg paths -Ev - > q.len
is $(cat q.len | wc -l) 2 "vg paths found correct number of lengths"
diff q.cov.len q.len
is $? 0 "vg path coverage reports correct lengths in first column"

rm -f q.vg q.cov.len q.len

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
is $? 0 "path normalizere correctly snapped all equivalent paths to x2"

rm -f norm_x2.gfa original.fa norm_x2.fa x2.path x2.norm.path

vg paths -x graphs/path_norm_test.gfa -n -Q x4 > norm_x4.gfa
vg validate norm_x4.gfa
is $? 0 "path normalizer produces valid graph"

vg paths -x graphs/path_norm_test.gfa -F | sort > original.fa
vg paths -x norm_x4.gfa -F | sort > norm_x4.fa
diff original.fa norm_x4.fa
is $? 0 "path normalizer doesnt alter path sequences"

# note: x3 is x4 in reverse, so we key on that
grep x3 norm_x2.gfa | awk '{print $3}' > x4.path
grep x3 norm_x2.gfa | awk '{print $3}' >> x4.path
grep x3 norm_x2.gfa | awk '{print $3}'> x4.norm.path
grep x5 norm_x2.gfa | awk '{print $3}' >> x4.norm.path
diff x4.path x4.norm.path
is $? 0 "path normalizere correctly snapped all equivalent paths to x4"

rm -f norm_x4.gfa original.fa norm_x4.fa x4.path x4.norm.path


   
   




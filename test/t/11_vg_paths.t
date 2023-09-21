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

vg paths -v rgfa/rgfa_tiny.gfa -R 1 -Q x | vg convert - -fW > rgfa_tiny.rgfa
printf "P	_rGFA_#y[33]	10+	*
P	_rGFA_#y[38]	13+	*
P	_rGFA_#y[8]	2+,4+	*\n" > rgfa_tiny_expected_fragments.rgfa
grep ^P rgfa_tiny.rgfa | grep rGFA | sort > rgfa_tiny_fragments.rgfa
diff rgfa_tiny_fragments.rgfa rgfa_tiny_expected_fragments.rgfa
is $? 0 "Found the expected rgfa SNP cover of tiny graph"

rm -f rgfa_tiny.rgfa rgfa_tiny_expected_fragments.rgfa rgfa_tiny_fragments.rgfa

vg paths -v rgfa/rgfa_ins.gfa -R 5 -Q x | vg convert - -fW >  rgfa_ins.rgfa
printf "P	_rGFA_#z[8]	2+,3+,4+	*\n" > rgfa_ins_expected_fragments.rgfa
grep ^P rgfa_ins.rgfa | grep _rGFA_ | sort > rgfa_ins_fragments.rgfa
diff rgfa_ins_fragments.rgfa rgfa_ins_expected_fragments.rgfa
is $? 0 "Found the expected rgfa cover for simple nested insertion"

rm -f rgfa_ins.rgfa rgfa_ins_expected_fragments.rgfa rgfa_ins_fragments.rgfa

vg paths -v rgfa/rgfa_ins2.gfa -R 3 -Q x | vg convert - -fW > rgfa_ins2.rgfa
printf "P	_rGFA_#y[8]	2+,6+,4+	*
P	_rGFA_#z[11]	3+	*\n" > rgfa_ins2_expected_fragments.rgfa
grep ^P rgfa_ins2.rgfa | grep _rGFA_ | sort > rgfa_ins2_fragments.rgfa
diff rgfa_ins2_fragments.rgfa rgfa_ins2_expected_fragments.rgfa
is $? 0 "Found the expected rgfa cover for simple nested insertion that requires two fragments"

rm -f rgfa_ins2.rgfa rgfa_ins2_expected_fragments.rgfa rgfa_ins2_fragments.rgfa

vg paths -v rgfa/rgfa_ins2.gfa -R 3 -Q x | grep ^S | sort -nk2 > rgfa_ins2.rgfa
printf "S	1	CAAATAAG	SN:Z:x	SO:i:0	SR:i:0
S	2	TTT	SN:Z:y	SO:i:8	SR:i:1
S	3	TTT	SN:Z:z	SO:i:11	SR:i:2
S	4	TTT	SN:Z:y	SO:i:21	SR:i:1
S	5	TTT	SN:Z:x	SO:i:8	SR:i:0
S	6	TTTTTTTTTT	SN:Z:y	SO:i:11	SR:i:1\n" > rgfa_ins2_expected_segs.rgfa
diff rgfa_ins2_expected_segs.rgfa rgfa_ins2.rgfa
is $? 0 "Found the expected rgfa tags for simple nested insertion that requires two fragments"

rm -f rgfa_ins2_expected_segs.rgfa rgfa_ins2.rgfa

vg paths -v rgfa/rgfa_ins2.gfa -R 5 -Q x | vg convert - -fW > rgfa_ins2R5.rgfa
printf "P	_rGFA_#y[8]	2+,6+,4+	*\n" > rgfa_ins2R5_expected_fragments.rgfa
grep ^P rgfa_ins2R5.rgfa | grep _rGFA_ | sort > rgfa_ins2R5_fragments.rgfa
diff rgfa_ins2R5_fragments.rgfa rgfa_ins2R5_expected_fragments.rgfa
is $? 0 "rgfa Minimum fragment length filters out small fragment"

rm -f rgfa_ins2R5.rgfa rgfa_ins2R5_expected_fragments.rgfa rgfa_ins2R5_fragments.rgfa

vg paths -v rgfa/rgfa_ins3.gfa -R 3 -Q x | vg convert - -fW  > rgfa_ins3.rgfa
printf "P	x	1+,5+	*
P	_rGFA_#y[3]	4+,6+,2+	*
P	y	5-,4+,6+,2+,1-	*
P	_rGFA_#z[11]	3+	*
P	z	1+,2-,3+,4-,5+	*\n" | sort > rgfa_ins3_expected_fragments.rgfa
grep ^P rgfa_ins3.rgfa | sort > rgfa_ins3_fragments.rgfa
diff rgfa_ins3_fragments.rgfa rgfa_ins3_expected_fragments.rgfa
is $? 0 "Found the expected rgfa cover for simple nested insertion that has a reversed path that needs forwardization"

rm -f rgfa_ins3.rgfa rgfa_ins3_expected_fragments.rgfa rgfa_ins3_fragments.rgfa


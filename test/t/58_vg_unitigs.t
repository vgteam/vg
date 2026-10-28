#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order

plan tests 22

# canonical_kmers K FASTA...: every canonical k-mer of the sequences (the smaller of the
# k-mer and its reverse complement), one per line.  The test graphs have random sequence,
# so k-mers that are equal as sequences are also the same graph positions.
canonical_kmers() {
    local k=$1
    shift
    cat "$@" | awk '/^>/ { if (seq != "") print seq; seq = ""; next } { seq = seq $0 } END { if (seq != "") print seq }' \
        | awk -v k=$k '{ for (i = 1; i <= length($0) - k + 1; i++) print substr($0, i, k) }' > kmers.fwd
    rev kmers.fwd | tr ACGT TGCA > kmers.rev
    paste kmers.fwd kmers.rev | awk '{ print ($1 < $2) ? $1 : $2 }'
    rm -f kmers.fwd kmers.rev
}

# Plain unitigs on a graph with a triallelic SNP and no gref
vg gbwt -G nesting/core_snp_in_flanks.gfa --gbz-format -g snp.gbz
vg paths -x snp.gbz -F > snp.paths.fa
vg unitigs -k 11 -g snp.gaf snp.gbz > snp.tigs.fa
is $? 0 "vg unitigs runs on a GBZ"

is "$(grep -c '^>' snp.tigs.fa)" "5" "a triallelic SNP gives two flanking unitigs and one per allele"

is "$(grep -v '^>' snp.tigs.fa | awk 'length($0) == 21 { n++ } END { print n + 0 }')" "3" "each allele's unitig is 2k-1 bp"

diff <(canonical_kmers 11 snp.paths.fa | sort -u) <(canonical_kmers 11 snp.tigs.fa | sort -u) > /dev/null
is $? 0 "the unitigs contain exactly the k-mers of the paths"

is "$(canonical_kmers 11 snp.tigs.fa | sort | uniq -d | awk 'END { print NR }')" "0" "no k-mer is in two unitigs"

is "$(awk 'END { print NR }' snp.gaf)" "5" "the GAF has one line per unitig"

is "$(cut -f 13 snp.gaf | grep -c '^sp:i:')" "5" "every GAF line has the haplotype support"

vg unitigs -k 11 -t 1 snp.gbz > snp.t1.fa
vg unitigs -k 11 -t 4 snp.gbz > snp.t4.fa
cmp -s snp.t1.fa snp.t4.fa
is $? 0 "the output does not depend on the thread count"

vg unitigs -k 1 snp.gbz > /dev/null 2>&1
is $? 1 "k below 2 is rejected"

vg unitigs -r -k 11 snp.gbz 2>&1 > /dev/null | grep -q "no gref paths"
is $? 0 "--gref on a graph without gref paths is an error"

# Junction tigs on a gref with a 50 bp insertion, which gets its own _alt path
vg gbwt -G nesting/core_ins50.gfa --gbz-format -g ins.gbz
vg paths -x ins.gbz -F > ins.paths.fa
vg paths -x nesting/core_ins50.gfa -Q x --compute-gref > ins.gref.vg
vg convert -f ins.gref.vg > ins.gref.gfa
vg gbwt -G ins.gref.gfa --gbz-format -g ins.gref.gbz
vg paths -x ins.gref.gbz -F -Q gref_ > ins.gref.fa
vg unitigs -r -k 11 ins.gref.gbz > ins.tigs.fa 2> ins.tigs.err
is $? 0 "vg unitigs --gref runs on a GBZ with gref paths"

diff <(canonical_kmers 11 ins.paths.fa | sort -u) <(canonical_kmers 11 ins.gref.fa ins.tigs.fa | sort -u) > /dev/null
is $? 0 "the gref and the junction tigs together contain exactly the k-mers of the paths"

is "$(canonical_kmers 11 ins.gref.fa ins.tigs.fa | sort | uniq -d | awk 'END { print NR }')" "0" "no k-mer is in both the gref and a junction tig"

is "$(grep -c '^>.* gref_x:61-70:+:10=,gref_x_1_alt:1-10:+:10=$' ins.tigs.fa)" "1" "a junction tig is named by its placement on the reference and the alt path"

is "$(grep -c 'not contained in any other path' ins.tigs.err)" "0" "no warning while every gref path is contained in its source path"

# Small variants are named as edits against the reference
vg paths -x nesting/core_snp_in_flanks.gfa -Q x --compute-gref > snp.gref.vg
vg convert -f snp.gref.vg > snp.gref.gfa
vg gbwt -G snp.gref.gfa --gbz-format -g snp.gref.gbz
vg unitigs -r -k 11 snp.gref.gbz > snp.jt.fa
is "$(grep -c '^>.* gref_x:91-111:+:10=1X10=$' snp.jt.fa)" "2" "each SNP allele off the reference is a substitution"

vg paths -x nesting/core_del61.gfa -Q x --compute-gref > del.gref.vg
vg convert -f del.gref.vg > del.gref.gfa
vg gbwt -G del.gref.gfa --gbz-format -g del.gref.gbz
vg unitigs -r -k 11 del.gref.gbz > del.jt.fa
is "$(grep -c '^>.* gref_x:[0-9]*-[0-9]*:+:10=6[01]D10=$' del.jt.fa)" "2" "each deletion allele is a deletion"

# Plain unitigs count every path as a haplotype, gref paths included
vg paths -x ins.gref.gbz -F > ins.gref.paths.fa
vg unitigs -k 11 ins.gref.gbz > ins.plain.fa
diff <(canonical_kmers 11 ins.paths.fa | sort -u) <(canonical_kmers 11 ins.plain.fa | sort -u) > /dev/null
is $? 0 "plain unitigs on a gref GBZ have the same k-mers as the paths"

# A gref whose source paths are gone, as after subsampling: its k-mers are in no haplotype
grep -v -E "^(P|W)"$'\t'"hap" ins.gref.gfa > nohaps.gfa
vg gbwt -G nohaps.gfa --gbz-format -g nohaps.gbz
vg unitigs -r -k 11 nohaps.gbz > /dev/null 2> nohaps.err
is $? 0 "vg unitigs still runs when the gref outlives its source paths"

is "$(grep -c 'not contained in any other path' nohaps.err)" "1" "a gref path that no other path contains is reported"

is "$(grep -c 'gref_x_1_alt' nohaps.err)" "1" "the report names an unsupported gref path"

# Long k, where no unitig can form, gives an empty cover
vg unitigs -k 1000 snp.gbz > long.fa
[ -s long.fa ]
is $? 1 "k longer than every path gives no unitigs"

rm -f snp.gbz snp.paths.fa snp.tigs.fa snp.gaf snp.t1.fa snp.t4.fa
rm -f ins.gbz ins.paths.fa ins.gref.vg ins.gref.gfa ins.gref.gbz ins.gref.fa ins.tigs.fa ins.tigs.err
rm -f snp.gref.vg snp.gref.gfa snp.gref.gbz snp.jt.fa del.gref.vg del.gref.gfa del.gref.gbz del.jt.fa
rm -f ins.gref.paths.fa ins.plain.fa nohaps.gfa nohaps.gbz nohaps.err long.fa

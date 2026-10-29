#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order

plan tests 16

# A triallelic SNP with the haplotypes stored as W lines at offset 1000 of contig c,
# as Minigraph-Cactus stores haplotype fragments
awk -v OFS='\t' '$1 == "P" && $2 ~ /^hap/ {
    n = split($3, steps, ",")
    walk = ""
    for (i = 1; i <= n; i++) {
        walk = walk (substr(steps[i], length(steps[i])) == "+" ? ">" : "<") substr(steps[i], 1, length(steps[i]) - 1)
    }
    print "W", "s", ($2 == "hapA" ? 1 : 2), "c", 1000, 1201, walk
    next
} { print }' nesting/core_snp_in_flanks.gfa > lift.gfa
vg paths -x lift.gfa -Q x --compute-gref > lift.gref.vg
vg convert -f lift.gref.vg > lift.gref.gfa
vg gbwt -G lift.gref.gfa --gbz-format -g lift.gbz
vg unitigs -r -k 11 -g lift.tigs.gaf lift.gbz > lift.tigs.fa

# The same 5 bases on the gref and on both junction tigs, a stranded feature across the
# SNP, a feature with no BED name, and a duplicate line
printf 'gref_x\t92\t97\tg5\t0\t+\nu1\t2\t7\tg5\t0\t+\nu2\t2\t7\tg5\t0\t+\n' > lift.bed
printf 'u1\t8\t13\tsnpA\t5\t-\nu2\t8\t13\tsnpB\nu2\t8\t13\tsnpB\ngref_x\t0\t5\n' >> lift.bed
vg lift -u lift.tigs.gaf -a lift.asm.bed lift.gbz lift.bed > lift.gaf
is $? 0 "vg lift runs on a cover of gref paths and junction tigs"

is "$(awk 'END { print NR }' lift.gaf)" "4" "copies of a feature on several cover sequences are kept once"

is "$(grep -c $'\tsp:i:3\tgr:Z:gref_x:93-97:+:5=\tan:Z:g5$' lift.gaf)" "1" "a feature on the gref is placed at its gref coordinates"

is "$(grep -c $'<6<4<2\t.*\tgr:Z:gref_x:99-103:-:2=1X2=\tan:Z:snpA$' lift.gaf)" "1" "a feature on the - strand of a tig is lifted to the reverse walk, with the SNP as a substitution"

is "$(awk '$4 == "g5"' lift.asm.bed | cut -f 1-3,6 | sort | tr '\t\n' ' ;')" "s#1#c 1092 1097 +;s#2#c 1092 1097 +;x 92 97 +;" "a haplotype place is on its contig, past the fragment's offset"

is "$(awk '$4 == "snpA"' lift.asm.bed | cut -f 1-6 | tr '\t' ' ')" "s#1#c 1098 1103 snpA 5 -" "a stranded feature keeps its fields and gets the strand of the path"

vg paths -x lift.gbz -F -Q 's#1#c' | awk 'NR > 1' | tr -d '\n' | cut -c 99-103 > lift.place.txt
awk '$1 == ">u1" { getline; print substr($0, 9, 5) }' lift.tigs.fa > lift.feature.txt
is "$(cat lift.place.txt) $(awk '{ print length($0) }' lift.place.txt)" "$(cat lift.feature.txt) 5" "the bases at a place on a haplotype are those of the feature on the tig"

is "$(awk '$4 == "." && $5 == 0 && $6 == "."' lift.asm.bed | awk 'END { print NR }')" "3" "a feature without a name, score or strand gets ., 0 and ."

join -t $'\t' <(cut -f 1,13 lift.gaf | sed 's/sp:i://' | sort) <(awk -v OFS='\t' '{ print $NF }' lift.asm.bed | sort | uniq -c | awk -v OFS='\t' '{ print $2, $1 }') > lift.support.txt
is "$(awk -F '\t' '$2 == $3 { n++ } END { print n + 0 }' lift.support.txt) $(cut -f 1 lift.gaf | sort | uniq -d | awk 'END { print NR }')" "4 0" "each feature has a unique ID and its support counts its places on paths"

vg lift -u lift.tigs.gaf -a lift.s1.bed -Q 's#1' lift.gbz lift.bed > /dev/null
is "$(cut -f 1 lift.s1.bed | sort -u)" "s#1#c" "places can be limited to paths by name prefix"

vg lift -u lift.tigs.gaf -t 1 -a lift.t1.bed lift.gbz lift.bed > lift.t1.gaf
vg lift -u lift.tigs.gaf -t 4 -a lift.t4.bed lift.gbz lift.bed > lift.t4.gaf
cmp -s lift.t1.gaf lift.t4.gaf && cmp -s lift.t1.bed lift.t4.bed
is $? 0 "the output does not depend on the thread count"

printf 'u99\t0\t5\n' > lift.unknown.bed
vg lift -u lift.tigs.gaf lift.gbz lift.unknown.bed > /dev/null 2>&1
is $? 1 "a feature on an unknown sequence is an error"

printf 'u1\t10\t30\n' > lift.long.bed
vg lift -u lift.tigs.gaf lift.gbz lift.long.bed 2>&1 > /dev/null | grep -q "past the end"
is $? 0 "a feature past the end of its sequence is an error"

vg lift lift.gbz lift.bed 2>&1 > /dev/null | grep -q "use -u"
is $? 0 "a feature on a tig needs the tig GAF"

# A plain cover, without a gref
vg gbwt -G lift.gfa --gbz-format -g lift.plain.gbz
vg unitigs -k 11 -g lift.plain.gaf lift.plain.gbz > lift.plain.fa
printf 'u1\t0\t5\tp\nu3\t9\t14\tq\n' > lift.plain.bed
vg lift -u lift.plain.gaf -a lift.plain.asm.bed lift.plain.gbz lift.plain.bed > lift.plain.out.gaf
is "$(grep -c 'gr:Z:' lift.plain.out.gaf)" "0" "features on a plain cover have no gref placement"

is "$(awk 'END { print NR }' lift.plain.asm.bed)" "4" "features on a plain cover are lifted to every path that visits them"

rm -f lift.gfa lift.gref.vg lift.gref.gfa lift.gbz lift.tigs.gaf lift.tigs.fa lift.bed lift.gaf lift.asm.bed
rm -f lift.place.txt lift.feature.txt lift.support.txt lift.s1.bed lift.t1.bed lift.t4.bed lift.t1.gaf lift.t4.gaf
rm -f lift.unknown.bed lift.long.bed lift.plain.gbz lift.plain.gaf lift.plain.fa lift.plain.bed lift.plain.out.gaf lift.plain.asm.bed

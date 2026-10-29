#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order

plan tests 30

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

is "$(grep -c $'\tsp:i:3\tgr:Z:gref_x:93-97:+:5=\tpc:Z:gref_x:93-97:+\tan:Z:g5$' lift.gaf)" "1" "a feature on the gref is placed at its gref coordinates, which are also its path coordinate"

is "$(grep -c $'<6<4<2\t.*\tgr:Z:gref_x:99-103:-:2=1X2=\tpc:Z:s#1#c:1099-1103:-\tan:Z:snpA$' lift.gaf)" "1" "a feature on the - strand of a tig is lifted to the reverse walk, with the SNP as a substitution"

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

# Path coordinates: a feature on the gref and on tigs, features off the gref on tigs, and
# features on a path, along the gref and off it
cp lift.bed lift.names.input.bed
printf 'u2\t9\t11\tcg\nx\t95\t105\talonggref\ns#1#c\t1098\t1103\toffgref\n' >> lift.names.input.bed
vg lift -u lift.tigs.gaf -c lift.names.bed lift.gbz lift.names.input.bed > lift.names.gaf 2> lift.names.err
is "$(grep -c $'\tpc:Z:s#2#c:1100-1101:+\tan:Z:cg$' lift.names.gaf)" "1" "a feature off the gref is named on the path that carries it, in contig coordinates"

is "$(grep -c $'\tpc:Z:gref_x:96-105:+\tan:Z:alonggref$' lift.names.gaf)" "1" "a feature read from a path along the gref is named on the gref"

is "$(grep -c $'\tpc:Z:s#1#c:1099-1103:+\tan:Z:offgref$' lift.names.gaf) $(grep -c . lift.names.err)" "1 0" "a feature read from a path off the gref is named on a path that carries it"

is "$(cut -f 1-4 lift.names.bed | tr '\t\n' ' ;')" "gref_x 0 5 .;gref_x 92 97 g5;gref_x 95 105 alonggref;s#1#c 1098 1103 offgref;s#1#c 1098 1103 snpA;s#2#c 1098 1103 snpB;s#2#c 1099 1101 cg;" "the canonical BED has each feature once, at its path coordinate, in priority order"

grep -v '^@' lift.tigs.gaf > lift.tigs.noheader.gaf
vg lift -u lift.tigs.noheader.gaf lift.gbz lift.names.input.bed 2> /dev/null | cmp -s - lift.names.gaf
is $? 0 "a tig GAF without a header gives the same output"

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

# Path to path: features on paths, in full-path coordinates
printf 'x\t80\t120\tspan\t0\t+\ns#1#c\t1080\t1120\tfromA\t0\t-\n' > lift.paths.bed
vg lift -a lift.paths.exact.bed lift.gbz lift.paths.bed > /dev/null
is "$(awk '$4 == "span"' lift.paths.exact.bed | cut -f 1-3)" "$(printf 'x\t80\t120')" "a feature across a SNP has only its own exact place"

vg lift -m 0.9 -a lift.paths.lift.bed lift.gbz lift.paths.bed > /dev/null
is "$(awk '$4 == "span" { print $1, $2, $3, $6, $NF }' lift.paths.lift.bed | tr '\n' ';')" "x 80 120 + 1;s#1#c 1080 1120 + 0.975;s#2#c 1080 1120 + 0.975;" "a liftover reaches the haplotypes with a SNP, at their contig coordinates"

is "$(awk '$4 == "fromA" && $1 == "x" { print $2, $3, $6 }' lift.paths.lift.bed)" "80 120 -" "a feature on a haplotype contig lifts to the reference on its strand"

vg lift -m 0.9 -b -a lift.paths.blocks.bed lift.gbz lift.paths.bed > /dev/null
is "$(awk '$4 == "span" && $1 == "s#1#c" { print $7, $8, $10, $11, $12 }' lift.paths.blocks.bed)" "1080 1120 2 20,19, 0,21," "blocks leave out the SNP base"

printf 's#1#c\t10\t20\tbefore\n' > lift.outside.bed
vg lift lift.gbz lift.outside.bed 2>&1 > /dev/null | grep -q "1 features on paths are not within one fragment"
is $? 0 "a feature outside the fragments of its path is skipped with a warning"

# A 50 bp insertion on the haplotypes, after base 70 of x
vg gbwt -G nesting/core_ins50.gfa --gbz-format -g lift.ins.gbz
printf 'x\t50\t110\twide\n' > lift.ins.bed
vg lift -m 0.9 -b -a lift.ins.blocks.bed lift.ins.gbz lift.ins.bed > /dev/null
is "$(awk '$1 == "hapA" { print $2, $3, $10, $11, $12, $NF }' lift.ins.blocks.bed)" "50 160 2 20,40, 0,70, 1" "a liftover spans an insertion shorter than the feature"

printf 'x\t60\t100\tnarrow\n' > lift.ins.bed
vg lift -m 0.7 -a lift.ins.lift.bed lift.ins.gbz lift.ins.bed > /dev/null
is "$(awk '$1 == "hapA" { print $2, $3, $NF }' lift.ins.lift.bed)" "120 150 0.75" "an insertion longer than the feature splits it"

printf 'x\t0\t10\tshared\n' > lift.shared.bed
is "$(vg lift lift.ins.gbz lift.shared.bed | grep -o 'pc:Z:[^[:space:]]*') $(vg lift -P x lift.ins.gbz lift.shared.bed | grep -o 'pc:Z:[^[:space:]]*')" "pc:Z:hapA:1-10:+ pc:Z:x:1-10:+" "a feature carried by several paths is named on the first in priority order"

vg unitigs -k 11 -g lift.ins.tigs.gaf lift.ins.gbz > /dev/null
vg lift -u lift.ins.tigs.gaf lift.gbz lift.bed 2>&1 > /dev/null | grep -q "not compatible"
is $? 0 "tigs of another graph are rejected"

rm -f lift.paths.bed lift.paths.exact.bed lift.paths.lift.bed lift.paths.blocks.bed lift.outside.bed
rm -f lift.ins.gbz lift.ins.bed lift.ins.blocks.bed lift.ins.lift.bed lift.ins.tigs.gaf
rm -f lift.names.input.bed lift.names.bed lift.names.gaf lift.names.err lift.tigs.noheader.gaf lift.shared.bed

rm -f lift.gfa lift.gref.vg lift.gref.gfa lift.gbz lift.tigs.gaf lift.tigs.fa lift.bed lift.gaf lift.asm.bed
rm -f lift.place.txt lift.feature.txt lift.support.txt lift.s1.bed lift.t1.bed lift.t4.bed lift.t1.gaf lift.t4.gaf
rm -f lift.unknown.bed lift.long.bed lift.plain.gbz lift.plain.gaf lift.plain.fa lift.plain.bed lift.plain.out.gaf lift.plain.asm.bed

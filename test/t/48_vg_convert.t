#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 79

vg construct -r complex/c.fa -v complex/c.vcf.gz > c.vg
cat <(vg view c.vg | grep ^S | sort) <(vg view c.vg | grep L | uniq | wc -l) <(vg paths -v c.vg -E) > c.info

vg convert c.vg -x > c.xg
vg convert c.xg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout xg conversion"

rm -f c.xg c1.vg c1.info

vg convert c.vg -a > c.hg
vg convert c.hg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout hash-graph conversion"

rm -f c.hg c1.vg c1.info

vg convert c.vg -p > c.pg
vg convert c.pg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout packed-graph conversion"

rm -f c.pg c1.vg c1.info

vg convert c.vg -f > c.gfa
vg convert -g c.gfa -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout gfa conversion"

rm -f c.gfa c1.vg c1.info

vg convert c.vg -o > c.odgi
vg convert c.odgi -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout ODGI conversion"

rm -f c.vg c.odgi c1.vg c.info c1.info

# some less rigorous tests I made without noticing that the earlier ones had already been written
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg view x.vg > x.gfa

is "$(vg convert -a x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "hash graph conversion looks good"
is "$(vg convert -p x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "packed graph conversion looks good"
is "$(vg convert -v x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "vg conversion looks good"
is "$(vg convert -o x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "odgi conversion looks good"
is "$(vg convert -f x.vg | vg convert -g - | vg view - | wc -l)" "$(wc -l < x.gfa)" "gfa conversion looks good"
is "$(vg convert -x x.vg | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "xg conversion looks good"

is "$(vg convert -g -a x.gfa | vg view - | wc -l)" "$(wc -l < x.gfa)" "on disk gfa conversion looks good"
is "$(cat x.gfa | vg convert -g -a - | vg view - | wc -l)" "$(wc -l < x.gfa)" "streaming gfa conversion looks good"
is "$(vg convert -g -x x.gfa | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "gfa to xg conversion looks good"
is "$(vg convert -g -f x.gfa | vg convert -g - | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "gfa to gfa conversion looks good"

rm x.vg x.gfa
rm -f c.vg c.pg c1.vg c.info c1.info

vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index x.vg -g x.gcsa
vg sim -x x.vg -n 10 -s 23 -a > sim.gam
vg map -x x.vg -g x.gcsa -G sim.gam > sim-rm.gam
vg convert x.vg -G sim-rm.gam -t 1 > sim-rm.gaf
vg convert x.vg -F sim-rm.gaf -t 1 | vg convert x.vg -G - -t 1 > sim-rm2.gaf
diff sim-rm.gaf sim-rm2.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf makes same gaf twice"

vg convert x.vg -G sim-rm.gam | vg convert x.vg -F - | vg convert x.vg -G - | sort > sim-rm2-mt-sort.gaf
sort sim-rm2.gaf > sim-rm2-sort.gaf
diff sim-rm2-sort.gaf sim-rm2-mt-sort.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf gives same result multithreaded as with -t 1"

vg convert x.vg -G sim-rm.gam | bgzip | vg convert x.vg -F - | vg convert x.vg -G - | sort > sim-rm2-mtbg-sort.gaf
diff sim-rm2-sort.gaf sim-rm2-mtbg-sort.gaf
is "$?" 0 "vg convert gam -> gaf.gz -> gam -> gaf gives same result multithreaded as with -t 1"

# some snps and indels
vg map -s "TAATGGATATGTTAAGCTTTTTTTTTCTTTGATTTATTTGAAAAGACGTTTGACAATCTATCGGGTAATGTGGGGAAA" -x x.vg -g x.gcsa > mut.gam
# reverse complement of above
vg map -s "TTTCCCCACATTACCCGATAGATTGTCAAACGTCTTTTCAAATAAATCAAAGAAAAAAAAAGCTTAACATATCCATTA" -x x.vg -g x.gcsa >> mut.gam

vg convert x.vg -G mut.gam -t 1 > mut.gaf
vg convert x.vg -F mut.gaf -t 1 > mut-back.gam
vg view -a mut.gam | jq .path > mut.path
vg view -a mut-back.gam | jq .path > mut-back.path
# Json comparison that is not order dependent: https://stackoverflow.com/a/31933234
is $(jq --argfile a mut.path --argfile b mut-back.path -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "vg convert gam -> gaf -> gam produces same gam Paths with snps and indels"

vg view -a mut.gam | jq .sequence > mut.seq
vg view -a mut-back.gam | jq .sequence > mut-back.seq
diff mut.seq mut-back.seq
is "$?" 0 "vg convert gam -> gaf -> gam preserves sequence"

vg convert x.vg -G mut-back.gam -t 1 > mut-back.gaf
diff mut.gaf mut-back.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf makes same gaf twice in presence of indels and snps"
  
#hand-code cg example.  this is (for reference) mut.gaf:
printf "*	78	0	78	+	>20>21>23>24>26>27>29>30>32>33>35	102	22	101	71	79	60	AS:i:47	cs:Z::13*GA*GA:8+TTT:16*GA*TA:18-AC-T-AG:16\n" > mut.cs.gaf
#manually convert to cg:
printf "*	78	0	78	+	>20>21>23>24>26>27>29>30>32>33>35	102	22	101	71	79	60	AS:i:47	cg:Z:13M1X1X8M3I16M1X1X18M5D16M\n" > mut.cg.gaf
#this is what we expect back, mut.gaf where insertions and snps are converted to Ns:
printf "*	78	0	78	+	>20>21>23>24>26>27>29>30>32>33>35	102	22	101	71	79	60	AS:i:47	cs:Z::13*GN*GN:8+NNN:16*GN*TN:18-AC-T-AG:16\n" > mut.cs.exp.gaf
vg convert x.vg -F mut.cg.gaf -t 1 | vg convert x.vg -G - -t 1 > mut.cs.back.gaf
diff mut.cs.back.gaf mut.cs.exp.gaf
is "$?" 0 "vg convert cg-gaf -> gam -> cs-gaf gives expected output (snps converted to matches, insertion converted to Ns)"
rm -f mut.cs.gaf mut.cg.gaf mut.cs.exp.gaf

rm -f x.vg x.gcsa sim.gam sim-rm.gam sim-rm.gaf sim-rm2.gaf sim-rm2-mt-sort.gaf sim-rm2-mtbg-sort.gaf sim-rm2-sort.gaf mut.gam mut-back.gam mut.gaf mut-back.gaf mut.path mut-back.path mut.seq mut-back.seq

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz > z.vg 2> /dev/null
vg sim -n 10000 -s 23 -a -x z.vg > sim.gam
vg construct -r 1mb1kgp/z.fa > zflat.vg
vg index zflat.vg -g zflat.gcsa
vg map -x zflat.vg -g zflat.gcsa -G sim.gam > sim-map.gam
vg convert zflat.vg -G sim-map.gam | bgzip > sim-map.gaf.gz
vg convert zflat.vg -F sim-map.gaf.gz > sim-map-back.gam
vg view -a sim-map.gam | jq .sequence | sort > sim-map.sequence
vg view -a sim-map-back.gam | jq .sequence | sort > sim-map-back.sequence
diff sim-map.sequence sim-map-back.sequence
is "$?" 0 "vg convert gam -> gaf -> gam preserves sequences of 1mb1kgp simulated reads"

vg convert zflat.vg -G sim-map-back.gam | sort > sim-map-back.gaf
bgzip -dc sim-map.gaf.gz | sort > sim-map.gaf
diff sim-map-back.gaf sim-map.gaf
is "$?" 0 "vg convert gam -> gaf -> gam ->gaf makes same gaf each time on 1mb1kgp simulated reads"

printf '{"name": "split", "path": {"mapping": [{"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "1", "offset": "10"}}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "3", "offset": "5"}}]}}' | vg view -JaG - > split.gam
vg convert zflat.vg -G split.gam > split.gaf
is "$(awk '{print $13}' split.gaf)" "cs:Z::13-CCAGTGCTC-GCATC:2" "split alignment converted using deletions to represent internal offsets"
vg convert zflat.vg -F split.gaf | vg convert zflat.vg -G - > split-back.gaf
diff split.gaf split-back.gaf
is "$?" 0 "vg convert gam -> gaf ->gam -> gaf makes same gaf each time for split alignment"

rm -f z.vg zflat.vg sim.gam sim-map.gam sim-map-back.gam sim-map.gaf.gz sim-map.sequence sim-map-back.sequence sim-map-back.gaf sim-map.gaf split.gam split.gaf split-back.gaf

printf "H\tVN:Z:1.0
S\t73333\tGGTGGGCGAGGACCTCCACACGTGTCACCA
S\t73368\tGCCCCT
S\t72943\tGGCGACTCTTCAGCAAGCCCCTCCACACGTGT
S\t72940\tC
S\t72941\tGGCCAGGT
S\t73255\tACTCTTCAGCAGGCCCCTCTGGT
S\t72942\tGGGCGAGGACCTCCACACGTGTCACCAGGCCA
S\t73318\tTCAGCA
S\t73367\tA
S\t73271\tA
S\t73289\tC
S\t73317\tC\n" | vg convert -g - -p > soft.pg
printf '{"annotation": {"fragment_length": 242, "fragment_length_distribution": "-I 561.110526 -D 141.152986", "mapq_applied_cap": 23.832780374978935, "mapq_extended_cap": 15, "mapq_uncapped": 10.325140756048304, "secondary_scores": [187.15265582416521, 182.4063994408188]}, "identity": 0.89682539682539686, "mapping_quality": 10, "name": "ERR903030.2067", "path": {"mapping": [{"edit": [{"from_length": 6, "to_length": 6}], "position": {"is_reverse": true, "node_id": "72943", "offset": "26"}}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"is_reverse": true, "node_id": "72942"}, "rank": "1"}, {"edit": [{"from_length": 23, "to_length": 23}], "position": {"is_reverse": true, "node_id": "73255"}, "rank": "2"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "73271"}, "rank": "3"}, {"edit": [{"from_length": 8, "to_length": 8}], "position": {"is_reverse": true, "node_id": "72941"}, "rank": "4"}, {"edit": [{"from_length": 7, "to_length": 7}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 19, "to_length": 19}], "position": {"is_reverse": true, "node_id": "73333"}, "rank": "5"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "72940"}, "rank": "6"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "73289"}, "rank": "7"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"is_reverse": true, "node_id": "73368"}, "rank": "8"}, {"edit": [{"from_length": 1, "sequence": "G", "to_length": 1}], "position": {"is_reverse": true, "node_id": "73367"}, "rank": "9"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"is_reverse": true, "node_id": "73318"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "73317"}, "rank": "11"}, {"edit": [{"sequence": "GGGTGGCCTG", "to_length": 10}], "position": {"is_reverse": true, "node_id": "73317", "offset": "1"}, "rank": "12"}]}, "quality": "ISEhISEmJiUmJCYmJSYlJiYmJiYmJCYjJiYdIRsmIyUmJiYmJiElJg4PDx0PDyEkIhsYJCYiECQQJBAcGR0kHx0QGSQPJCImHSElHR0PJA4OGyQbIxwPIg0iHw8PGRwiIR0jAgICAgICAgICAgICAgICAgICAgICAgICAgIC", "sample_name": "HG00514_961a37c_gssw", "score": 106, "sequence": "GTCGCCTGGCCTGGTGACACGTGTGGAGGTCCTCGCCCACCAGAGGGGCCTGCTGAAGAGTTACCTGGCCTGGTGACCCGGGTGGAGGTCCTCGCCCACCGGAGGGGCGTGCTGAGGGGTGGCCTG"}' | vg view -JaG - > soft.gam
vg convert soft.pg -G soft.gam > soft.gaf
vg view -a soft.gam | jq .sequence > gam.sequence
vg convert soft.pg -F soft.gaf | vg view -a - | jq .sequence > gam2.sequence
diff gam.sequence gam2.sequence
is "$?" 0 "convert gam->gaf->gam on softclipped read preserves sequence"

vg convert soft.pg -F soft.gaf | vg convert soft.pg -G - > soft2.gaf
diff soft.gaf soft2.gaf
is "$?" 0 "convert gam->gaf->gam->gaf makes same gaf each time of soft clipped alignment" 

rm -f soft.pg soft.gam soft.gaf gam.sequence gam2.sequence soft2.gaf

printf "H\tVN:Z:1.0
S\t91194329\tAGGAAGGAGAGGGAG\n" | vg convert -g - -p > floating-ins.pg
printf '{"annotation": {"fragment_length": 1098, "fragment_length_distribution": "-I 542.973684 -D 141.206118", "mapq_applied_cap": 46.364361584299516, "mapq_extended_cap": "Infinity", "mapq_uncapped": 1.5051499783199018, "rescued": true, "secondary_scores": [99.415737573445895]}, "mapping_quality": 1, "name": "ERR903030.51990324", "path": {"mapping": [{"edit": [{"sequence": "GGGCACGGTGGCTCACAGCTGTCACCACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", "to_length": 126}], "position": {"is_reverse": true, "node_id": "91194329", "offset": "15"}, "rank": "1"}]}, "quality": "ISEhICEhJCIkJiYdJiUQJSYmHyMmIxAkJiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIC", "sample_name": "HG00514_961a37c", "sequence": "GGGCACGGTGGCTCACAGCTGTCACCACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}' | vg view -JaG - > floating-ins.gam
vg convert floating-ins.pg -G floating-ins.gam > floating-ins.gaf
vg view -a floating-ins.gam | jq .sequence > gam.sequence
vg convert floating-ins.pg -F floating-ins.gaf | vg view -a - | jq .sequence > gam2.sequence
diff gam.sequence gam2.sequence
is "$?" 0 "convert gam->gaf->gam on read with floating insertion preserves sequence"

vg convert floating-ins.pg -F floating-ins.gaf | vg convert floating-ins.pg -G - > floating-ins2.gaf
diff floating-ins.gaf floating-ins2.gaf
is "$?" 0 "convert gam->gaf->gam->gaf makes same gaf each time of floating insertion alignment" 

rm -f floating-ins.pg floating-ins.gam floating-ins.gaf gam.sequence gam2.sequence floating-ins2.gaf

vg convert -g tiny/tiny.gfa | vg view - | sort > tiny.gfa.gfa
vg convert -g tiny/tiny.rgfa | vg view - | sort > tiny.rgfa.gfa
diff tiny.gfa.gfa tiny.rgfa.gfa
is "$?" 0 "converting gfa and equivalent rgfa produces same output"

rm -f tiny.gfa.gfa tiny.rgfa.gfa

is "$(vg convert -g tiny/tiny.rgfa -r 1 | vg view - | grep y | awk '{print $1","$2","$3}')" "P,y[35],2+" "rank-1 rgfa subpath found"

vg convert -g tiny/tiny.rgfa -T gfa-id-mapping.tsv > /dev/null
grep ^S tiny/tiny.rgfa  | awk '{print $2}' | sort > rgfa_nodes
grep ^S tiny/tiny.gfa  | awk '{print $2}' | sort > gfa_nodes
awk '{print $2}' gfa-id-mapping.tsv | sort > rgfa_translated_nodes
awk '{print $3}' gfa-id-mapping.tsv | sort > gfa_translated_nodes
diff rgfa_nodes rgfa_translated_nodes
is "$?" 0 "2nd column of gfa id translation file contains all rgfa nodes"
diff gfa_nodes gfa_translated_nodes
is "$?" 0 "3rd column of gfa id translation file contains all gfa nodes"

rm -f  gfa-id-mapping.tsv rgfa_nodes gfa_nodes rgfa_translated_nodes gfa_translated_nodes

vg convert -g tiny/tiny.gfa -v | vg convert - -f -P x > tiny.gfa.rgfa
is "$(grep ^P tiny.gfa.rgfa | wc -l)" 0 "rgfa output wrote no P-lines"
vg convert -g tiny/tiny.gfa -v | vg convert - -f | sort > tiny.gfa.gfa
vg convert -g tiny.gfa.rgfa -f | sort > tiny.gfa.rgfa.gfa
diff tiny.gfa.gfa tiny.gfa.rgfa.gfa
is "$?" 0 "rgfa export roundtrips back to normal P-lines"

rm -f tiny.gfa.rgfa tiny.gfa.gfa tiny.gfa.rgfa.gfa


# GFA to GBWTGraph to HashGraph to GFA
vg gbwt -o components.gbwt -g components.gg -G graphs/components_walks.gfa
sort graphs/components_walks.gfa > correct.gfa
vg convert -b components.gbwt -a components.gg > components.hg 2> /dev/null
is $? 0 "GBWTGraph to HashGraph conversion"
grep "^S" graphs/components_walks.gfa | sort > sorted.gfa
vg view components.hg | grep "^S" | sort > converted.gfa
cmp sorted.gfa converted.gfa
is $? 0 "GFA -> GBWTGraph -> HashGraph -> GFA conversion maintains segments"

# GFA to GBZ to HashGraph to GFA
vg gbwt -g components.gbz --gbz-format -G graphs/components_walks.gfa
vg convert -a components.gbz > components.hg 2> /dev/null
is $? 0 "GBZ to HashGraph conversion"
vg view components.hg | grep "^S" | sort > converted.gfa
cmp sorted.gfa converted.gfa
is $? 0 "GFA -> GBZ -> HashGraph -> GFA conversion maintains segments"

# GBWTGraph to GFA with walks (needs 1 thread)
vg convert -b components.gbwt -f -t 1 components.gg > extracted.gfa
is $? 0 "GBWTGraph to GFA conversion with walks, GBWTGraph algorithm"
cmp extracted.gfa graphs/components_walks.gfa
is $? 0 "GBWTGraph to GFA conversion with GBWTGraph algorithm creates the correct normalized GFA file"
vg convert --vg-algorithm -b components.gbwt -f -t 1 components.gg | sort - > extracted.gfa
is $? 0 "GBWTGraph to GFA conversion with walks, vg algorithm"
cmp extracted.gfa correct.gfa
is $? 0 "GBWTGraph to GFA conversion with vg algorithm creates the correct possibly-unnormalized GFA file"

# GBZ to GFA with walks (needs 1 thread)
vg convert -f -t 1 components.gbz > extracted.gfa
is $? 0 "GBZ to GFA conversion with walks, GBWTGraph algorithm"
cmp extracted.gfa graphs/components_walks.gfa
is $? 0 "GBZ to GFA conversion with GBWTGraph algorithm creates the correct normalized GFA file"

rm -f components.gbwt components.gg components.gbz
rm -f components.hg
rm -f sorted.gfa converted.gfa correct.gfa
rm -f extracted.gfa


# GFA to GBWTGraph and GBZ with paths and walks
vg gbwt -o components.gbwt -g components.gg -G graphs/components_paths_walks.gfa
vg gbwt -g components.gbz --gbz-format -G graphs/components_paths_walks.gfa
vg convert -g -a graphs/components_paths_walks.gfa > direct.hg
vg paths --generic-paths -v direct.hg -A > correct_paths.gaf
sort graphs/components_paths_walks.gfa > correct.gfa

# GBWTGraph to HashGraph with paths and walks
vg convert -b components.gbwt -a components.gg > components.hg
is $? 0 "GBWTGraph to HashGraph conversion with generic paths"
vg paths --generic-paths -A -v components.hg > hg_paths.gaf
cmp hg_paths.gaf correct_paths.gaf
is $? 0 "GBWTGraph to HashGraph conversion creates the correct generic paths"

# GBZ to HashGraph with paths and walks
vg convert -a components.gbz > components.hg
is $? 0 "GBZ to HashGraph conversion with generic paths"
vg paths --generic-paths -A -v components.hg > gbz_hg_paths.gaf
cmp gbz_hg_paths.gaf correct_paths.gaf
is $? 0 "GBZ to HashGraph conversion creates the correct generic paths"

# GBWTGraph to XG with paths and walks
vg convert -b components.gbwt -x components.gg > components.xg
is $? 0 "GBWTGraph to XG conversion with generic paths"
vg paths --generic-paths -A -v components.xg > xg_paths.gaf
cmp xg_paths.gaf correct_paths.gaf
is $? 0 "GBWTGraph to XG conversion creates the correct generic paths"

# GBZ to XG with paths and walks
vg convert -x components.gbz > components.xg
is $? 0 "GBZ to XG conversion with generic paths"
vg paths --generic-paths -A -v components.xg > gbz_xg_paths.gaf
cmp gbz_xg_paths.gaf correct_paths.gaf
is $? 0 "GBZ to XG conversion creates the correct generic paths"

# GBWTGraph to GFA with paths and walks (needs 1 thread)
vg convert -b components.gbwt -f -t 1 components.gg > extracted.gfa
is $? 0 "GBWTGraph to GFA conversion with paths and walks, GBWTGraph algorithm"
cmp extracted.gfa graphs/components_paths_walks.gfa
is $? 0 "GBWTGraph to GFA conversion with GBWTGraph algorithm creates the expected normalized GFA file"
vg convert --vg-algorithm -b components.gbwt -f -t 1 components.gg | sort - > extracted.gfa
is $? 0 "GBWTGraph to GFA conversion with paths and walks, vg algorithm"
cmp extracted.gfa correct.gfa
is $? 0 "GBWTGraph to GFA conversion with vg algorithm creates the correct possibly-unnormalized GFA file"

# GBWTGraph to HashGraph to GFA with paths and walks
vg convert -b components.gbwt -t 1 components.gg -a > extracted.hg
is $? 0 "GBWTGraph to HashGraph conversion with paths and walks"
vg convert -f -t1 extracted.hg | sort - > extracted.gfa
is $? 0 "HashGraph to GFA conversion with paths and walks"
cmp extracted.gfa correct.gfa
is $? 0 "GBWTGraph to HashGraph to GFA conversion creates the correct possibly-unnormalized GFA file"
vg convert --no-wline -f -t1 extracted.hg > no-walks.gfa
is $? 0 "HashGraph to GFA conversion writing walks as paths"
is "$(grep "^W" no-walks.gfa | wc -l)" "0" "HashGraph to GFA conversion writing walks as paths produces no walks"
is "$(grep "^P" no-walks.gfa | wc -l)" ""$(grep "^[PW]" correct.gfa | wc -l)"" "HashGraph to GFA conversion writing walks as paths produces all expected paths"

# GBZ to GFA with paths and walks (needs 1 thread)
vg convert --gbwtgraph-algorithm  -f -t 1 components.gbz > gbz.gfa
is $? 0 "GBZ to GFA conversion with paths and walks, GBWTGraph algorithm"
cmp gbz.gfa graphs/components_paths_walks.gfa
is $? 0 "GBZ to GFA conversion with GBWTGraph algorithm creates the correct normalized GFA file"

# Multithreaded GBZ to GFA with paths and walks
vg convert -f components.gbz | sort > sorted.gfa
cmp sorted.gfa correct.gfa
is $? 0 "GBZ to GFA conversion works with multiple threads"

rm -f components.gbwt components.gg components.gbz
rm -f direct.hg correct_paths.gaf
rm -f components.hg hg_paths.gaf gbz_hg_paths.gaf
rm -f components.xg xg_paths.gaf gbz_xg_paths.gaf
rm -f extracted.gfa gbz.gfa extracted.hg
rm -f sorted.gfa correct.gfa


# GFA Streaming
vg convert -g tiny/tiny.gfa -p | vg convert -f - | sort > tiny.roundtrip.gfa
vg convert tiny/tiny.gfa -p | vg convert -f - | sort > tiny.roundtrip2.gfa
diff tiny.roundtrip.gfa tiny.roundtrip2.gfa
is $? 0 "No difference roundtripping a GFA if it's loaded as a GFA or HandleGraph"

grep -v "S	6" tiny/tiny.gfa > tiny.unsort.gfa
grep "S	6" tiny/tiny.gfa >> tiny.unsort.gfa
cat tiny.unsort.gfa | vg convert -p - 2> tiny.roundtrip3.stderr | vg convert -f - | sort > tiny.roundtrip3.gfa
cat tiny.roundtrip3.stderr
diff tiny.roundtrip.gfa tiny.roundtrip3.gfa
is $? 0 "Streaming an unsorted GFA gives same output as sorted"
is $(grep "warning:\[gfa\]" tiny.roundtrip3.stderr | wc -l) 1 "Warning given when falling back to temp GFA buffer file"

cat tiny/tiny.gfa | vg convert -p - 2> tiny.roundtrip4.stderr | vg convert -f - | sort > tiny.roundtrip4.gfa
cat tiny.roundtrip4.stderr
diff tiny.roundtrip.gfa tiny.roundtrip4.gfa
is $? 0 "Streaming an sorted GFA gives same output as reading from file"
is $(cat tiny.roundtrip4.stderr | wc -l) 0 "No warnings given when streamed GFA is sorted"

vg convert -g tiny/tiny.gfa | vg mod - -X 3 | vg convert -f - | vg ids -s - | sort > tiny.chop3.gfa
vg mod -X 3 tiny/tiny.gfa | vg ids -s - | sort > tiny.chop3.1.gfa
diff tiny.chop3.gfa tiny.chop3.1.gfa
is $? 0 "Modding GFA directly produces same output as going through convert"
cat tiny/tiny.gfa | vg mod -X 3 - | vg ids -s - | sort > tiny.chop3.2.gfa
diff tiny.chop3.gfa tiny.chop3.2.gfa
is $? 0 "Modding sorted GFA stream produces same output as going through convert"
cat tiny.unsort.gfa | vg mod -X 3 - 2> tiny.chop3.3.stderr | vg ids -s - | sort > tiny.chop3.3.gfa
cat tiny.chop3.3.stderr
diff tiny.chop3.gfa tiny.chop3.3.gfa
is $? 0 "Modding unsorted GFA stream produces same output as going through convert"
is $(grep "warning:\[gfa\]" tiny.chop3.3.stderr | wc -l) 1 "Warning given when falling back to temp GFA buffer file in mod"
vg mod -X 3 tiny.unsort.gfa 2> tiny.chop3.4.stderr | vg ids -s - | sort > tiny.chop3.4.gfa
cat tiny.chop3.4.stderr
diff tiny.chop3.gfa tiny.chop3.4.gfa
is $? 0 "Modding unsorted GFA file produces same output as going through convert"
is $(cat tiny.chop3.4.stderr | wc -l) 0 "No warnings given when input GFA file is unsorted"

rm -f tiny.roundtrip.gfa tiny.roundtrip2.gfa tiny.roundtrip3.gfa tiny.roundtrip4.gfa
rm -f tiny.roundtrip3.stderr tiny.roundtrip4.stderr
rm -f tiny.unsort.gfa
rm -f tiny.chop3.gfa tiny.chop3.1.gfa  tiny.chop3.2.gfa  tiny.chop3.3.gfa tiny.chop3.4.gfa
rm -f tiny.chop3.3.stderr tiny.chop3.4.stderr

vg view tiny/tiny.gfa | sort > tiny.rgfa.1
cat tiny/tiny.gfa | vg view - | sort > tiny.rgfa.2
diff tiny.rgfa.1 tiny.rgfa.2
is $? 0 "rGFA handled consistently when streaming as when loaded from file"




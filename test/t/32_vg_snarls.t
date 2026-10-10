#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 42

vg view -J -v snarls/snarls.json > snarls.vg
vg snarls -t 1 snarls.vg -r st.pb > snarls.pb
is $(vg view -R snarls.pb | wc -l) 3 "vg snarls made right number of protobuf Snarls"
is $(vg view -E st.pb | wc -l) 6 "vg snarls made right number of protobuf SnarlTraversals"
is $(vg view -R snarls.pb | jq -r '[(.start.node_id | tonumber), (.end.node_id | tonumber)] | min' | tr '\n' ',') "1,3,7," "vg snarls made snarls in the right order"

vg snarls -l -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "vg snarls does not allow --leaf-only to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
vg snarls -l -t 1 snarls.vg -r /dev/null > leaf_only.pb
diff leaf_only.pb snarls.pb
is $? 0 "vg snarls --leaf-only matches full snarls"

vg snarls -o -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "vg snarls does not allow --top-level to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
vg snarls -o -t 1 snarls.vg -r /dev/null > top_level.pb
diff top_level.pb snarls.pb
is $? 0 "vg snarls --top-level matches full snarls"

vg snarls -a -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "vg snarls does not allow --any-snarl-type to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
vg snarls -a -t 1 snarls.vg -r /dev/null > any.pb
diff any.pb snarls.pb
is $? 0 "vg snarls --any-snarl-type matches full snarls"

vg snarls -m 5 -t 1 snarls.vg > /dev/null 2> error.txt
is $? 1 "vg snarls does not allow --max-nodes to be used without --traversals"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"
vg snarls -m 5 -t 1 snarls.vg -r /dev/null > small.pb
diff small.pb snarls.pb
is $? 0 "vg snarls with --max-nodes changed matches full snarls"

rm -f snarls.pb st.pb leaf_only.pb top_level.pb any.pb small.pb error.txt

vg index snarls.vg -x snarls.xg
is $(vg snarls snarls.xg -r st.pb | vg view -R - | wc -l) 3 "vg snarls on xg made right number of protobuf Snarls"
is $(vg view -E st.pb | wc -l) 6 "vg snarls on xg made right number of protobuf SnarlTraversals"

rm -f snarls.vg snarls.xg st.pb

# vcf alt traversals in tiny graph
vg construct -Saf -v tiny/tiny.vcf.gz -r tiny/tiny.fa > tiny.vg
vg snarls tiny.vg -r tiny.exhaustive.trav > /dev/null
vg snarls tiny.vg -r tiny.vcf.trav -v tiny/tiny.vcf.gz > /dev/null
vg view -E tiny.exhaustive.trav | sort > tiny.exhaustive.trav.sort
vg view -E tiny.vcf.trav | sort > tiny.vcf.trav.sort
diff tiny.exhaustive.trav.sort tiny.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for tiny graph"

vg snarls tiny.vg -v tiny/tiny.vcf.gz > /dev/null 2> error.txt
is $? 1 "vg snarls requires --traversals to be given if --vcf is"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"

vg snarls tiny.vg -r tiny.vcf.trav -v tiny/tiny.vcf.gz -e > /dev/null 2> error.txt
is $? 1 "vg snarls does not allow both --vcf and --path-traversals"
grep "cannot be used with" error.txt
is $? 0 "the problem is explained"

vg snarls tiny.vg -e > /dev/null 2> error.txt
is $? 1 "vg snarls requires --traversals to be given if --path-traversals is"
grep "requires --traversals" error.txt
is $? 0 "the problem is explained"

vg snarls tiny.vg -e -m 5 -r /dev/null > /dev/null 2> warning.txt
grep "ignored" warning.txt
is $? 0 "warning about ignoring is given"

rm -f tiny.vg tiny.exhaustive.trav.sort tiny.exhaustive.trav tiny.vcf.trav.sort tiny.vcf.trav error.txt warning.txt

# vcf alt traversals on an inversion
vg construct -Saf -v sv/x.inv.vcf -r sv/x.fa > x.inv.vg
vg snarls x.inv.vg -a -r x.inv.exhaustive.trav > /dev/null
vg snarls x.inv.vg -a -r x.inv.vcf.trav -v sv/x.inv.vcf -f sv/x.fa > /dev/null
vg view -E x.inv.exhaustive.trav | sort > x.inv.exhaustive.trav.sort
vg view -E x.inv.vcf.trav | sort > x.inv.vcf.trav.sort
diff x.inv.exhaustive.trav.sort x.inv.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for inversion"

vg snarls x.inv.vg -a -r x.inv.vcf.trav -f sv/x.fa > /dev/null 2> error.txt
is $? 1 "vg snarls requires --vcf to be given if --fasta is"
grep "requires --vcf" error.txt
is $? 0 "the problem is explained"

vg snarls x.inv.vg -a -r x.inv.vcf.trav -i sv/x.fa > /dev/null 2> error.txt
is $? 1 "vg snarls requires --vcf to be given if --ins-fasta is"
grep "requires --vcf" error.txt
is $? 0 "the problem is explained"

rm -f x.inv.vg x.inv.exhaustive.trav.sort x.inv.exhaustive.trav x.inv.vcf.trav.sort x.inv.vcf.trav error.txt

# vcf alt traversals on ins_and_del
vg construct -Saf -v ins_and_del/ins_and_del.vcf.gz -r ins_and_del/ins_and_del.fa > ins_and_del.vg
vg snarls ins_and_del.vg -r ins_and_del.exhaustive.trav > /dev/null
vg snarls ins_and_del.vg -r ins_and_del.vcf.trav -v ins_and_del/ins_and_del.vcf.gz > /dev/null
vg view -E ins_and_del.exhaustive.trav | sort > ins_and_del.exhaustive.trav.sort
vg view -E ins_and_del.vcf.trav | sort > ins_and_del.vcf.trav.sort
diff ins_and_del.exhaustive.trav.sort ins_and_del.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for ins_and_del graph"

rm -f ins_and_del.vg ins_and_del.exhaustive.trav.sort ins_and_del.exhaustive.trav ins_and_del.vcf.trav.sort ins_and_del.vcf.trav

# parallelizing on components (which is deactivated with -t 1)
vg construct -r small/xy.fa -v small/xy.vcf.gz > xy.vg
vg construct -r small/xy.fa -v small/xy.vcf.gz -R x > x.vg
vg construct -r small/xy.fa -v small/xy.vcf.gz -R y > y.vg
vg snarls x.vg > xy.snarls
vg snarls y.vg >> xy.snarls
is $(vg snarls xy.vg | vg view -R - | wc -l) 35 "correct number of snarls when parallelizing on compoents"
is $(vg snarls xy.vg | vg view -R - | wc -l) $(vg view -R xy.snarls | wc -l) "same number of snarls when parallelizing on components"
rm -f xy.vg xy.snarls x.vg y.vg

# Find trivial snarls from a GBZ
vg gbwt -g graph.gbz --gbz-format -G graphs/components_walks.gfa
vg snarls -T graph.gbz > graph.snarls
is $? 0 "GBZ graphs can be used for finding snarls"
is $(vg view -R graph.snarls | wc -l) 5 "correct number of snarls in the GFA W-line example"
rm -f graph.gbz graph.snarls

# Check snarl output order in a more complex case.
# Snarls need to come out in order along chains, recursing down and then coming back up to the right place in the chain.
vg view -J -v snarls/nested.json > nested.vg
vg snarls -t 1 -T nested.vg > nested.snarls
# Could go either way first in the first snarl
echo "1,2,5,5,7,3,9," > possibilities.txt
echo "1,3,2,5,5,7,9," >> possibilities.txt
GOT="$(vg view -R nested.snarls | jq -r '[(.start.node_id | tonumber), (.end.node_id | tonumber)] | min' | tr '\n' ',')"
is $(echo $GOT | wc -c) 15 "vg snarls made the right snarls when dealing with nested chains"
grep $GOT possibilities.txt >/dev/null
is $? 0 "vg snarls made snarls in the right order when dealing with nested chains"
rm -f nested.vg nested.snarls possibilities.txt

# Find snarls from a GFA in GFA space
is "$(vg snarls --named-coordinates graphs/components_walks_named.gfa | vg view -Rj - | jq -c '[.start.name, .start.backward, .end.name, .end.backward]' | sort)" \
   "$(printf '["cats",null,"goldfish",null]\n["goldfish",null,"sparrows",null]\n["pigs",true,"squirrels",null]\n["squirrels",null,"rabbits",true]\n')" \
   "Correct snarls are found in GFA space"

# Find traversals from a chopped GBZ-ified GFA in node ID space
vg autoindex -p index -g graphs/big_snarl_named.gfa -w giraffe
vg snarls index.giraffe.gbz --traversals traversals.dat >/dev/null
is "$(vg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.node_id // .name)' | wc -l)" "7" "In node ID space traversal visits 7 nodes"
is "$(vg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.snarl)' | wc -l)" "4" "In node ID space traversal visits 4 trivial snarls"

# And in GFA space
vg snarls index.giraffe.gbz --traversals traversals.dat --named-coordinates >/dev/null
is "$(vg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.node_id // .name)' | wc -l)" "3" "In segment name space traversal visits 3 segments"
is "$(vg view -Ej traversals.dat | jq -c 'select(.visit | length > 2) | .visit[] | select(.snarl)' | wc -l)" "0" "In segment name space traversal visits 0 trivial snarls"

rm -f index.* traversals.dat

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 7

vg view -J -v snarls/snarls.json > snarls.vg
is $(vg snarls snarls.vg -r st.pb | vg view -R - | wc -l) 3 "vg snarls made right number of protobuf Snarls"
is $(vg view -E st.pb | wc -l) 6 "vg snarls made right number of protobuf SnarlTraversals"

rm -f st.pb

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

rm -f tiny.vg tiny.exhaustive.trav.sort tiny.exhaustive.trav tiny.vcf.trav.sort tiny.vcf.trav

# vcf alt traversals on an inversion
vg construct -Saf -v sv/x.inv.vcf -r sv/x.fa > x.inv.vg
vg snarls x.inv.vg -a -r x.inv.exhaustive.trav > /dev/null
vg snarls x.inv.vg -a -r x.inv.vcf.trav -v sv/x.inv.vcf -f sv/x.fa > /dev/null
vg view -E x.inv.exhaustive.trav | sort > x.inv.exhaustive.trav.sort
vg view -E x.inv.vcf.trav | sort > x.inv.vcf.trav.sort
diff x.inv.exhaustive.trav.sort x.inv.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for inversion"

rm -f x.inv.vg x.inv.exhaustive.trav.sort x.inv.exhaustive.trav x.inv.vcf.trav.sort x.inv.vcf.trav

# vcf alt traversals on ins_and_del
vg construct -Saf -v ins_and_del/ins_and_del.vcf.gz -r ins_and_del/ins_and_del.fa > ins_and_del.vg
vg snarls ins_and_del.vg -r ins_and_del.exhaustive.trav > /dev/null
vg snarls ins_and_del.vg -r ins_and_del.vcf.trav -v ins_and_del/ins_and_del.vcf.gz > /dev/null
vg view -E ins_and_del.exhaustive.trav | sort > ins_and_del.exhaustive.trav.sort
vg view -E ins_and_del.vcf.trav | sort > ins_and_del.vcf.trav.sort
diff ins_and_del.exhaustive.trav.sort ins_and_del.vcf.trav.sort
is $? 0 "vcf traversals are the same as exhaustive traversals for ins_and_del graph"

rm -f ins_and_del.vg ins_and_del.exhaustive.trav.sort ins_and_del.exhaustive.trav ins_and_del.vcf.trav.sort ins_and_del.vcf.trav


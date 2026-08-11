#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 7

vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > linear.vg
vg circularize -p x linear.vg > circular.vg

is $(vg view -j circular.vg | jq -c '.path[] | select(.is_circular)' | wc -l) 1 "a path may be circularized"

vg index -x circular.xg circular.vg
vg convert circular.xg -v  > extracted.vg

is $(vg view -j extracted.vg | jq -c '.path[] | select(.is_circular)' | wc -l) 1 "a circular path survives a round trip to/from xg"

vg circularize -p y linear.vg
is $? 1 "Not allowed to circularize a nonexistent path (--path)"

echo "y" > paths.txt
vg circularize -P paths.txt linear.vg
is $? 1 "Not allowed to circularize a nonexistent path (--pathfile)"

vg circularize -a 2 -z 1 linear.vg
is $? 1 "Not allowed to have tail ID smaller than head ID"

vg circularize -a 1 linear.vg
is $? 1 "Not allowed to have only a head ID"

vg circularize -z 2 linear.vg
is $? 1 "Not allowed to have only a tail ID"

rm -f circular.vg circular.xg extracted.vg linear.vg paths.txt


#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

plan tests 3

cd ..

make clean
is "$?" "0" "make clean completes without errors"
make
is "$?" "0" "we can make after make clean"
bin/vg help
is "$?" "0" "vg builds successfully"

#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

YEETS="$(find ../src -iname "*.[ch]pp" | xargs grep "yeet" | wc -l)"
is "${YEETS}" 0 "code quality is acceptable"

#!/usr/bin/env bash
#
BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

vgtordf.sh cyclic/all.json | jq .
  ' 

#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 0

## Test that ggsv throws an error if nodes specified in a superbubble file are not found in a graph/index

## Test whether ggsv can depth filter.

## Test the depth filter

## Test whether GGSV successfully generates a valid VCF on a simple graph/index/supbub file.

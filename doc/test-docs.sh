#!/usr/bin/env bash
# test-docs.sh: Test the examples in the vg documentation and wiki with https://github.com/anko/txm
set -e

# Work out where we are.
# See https://stackoverflow.com/a/246128
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Make sure we have the tester
which txm || npm install -g txm

# Go to the test directory, where the tests expect to run
cd "${HERE}/../test"

# Test the readme
echo txm --jobs 1 "${HERE}/../README.md"
txm --jobs 1 "${HERE}/../README.md"

# Run all the wiki tests
find "${HERE}/wiki" -name "*.md" | xargs -n 1 -t txm --jobs 1



#!/usr/bin/env bash

# On CI, we want to cache deps/ and avoid rebuilding unmodified submodules, but
# we also want to make sure every build has the right source code for all the
# submodules.
#
# This script will start with the version of deps/ fetched from the cache, and
# put in place only the differing source files. Differing files will be newer
# than the files from the cache.
#
# This script is meant to run in the root of the repository.
#
# This script will not, on its own, tolerate moving between compiler versions.

set -ex

if [[ -e deps ]]; then
    # Get the cached deps out of the way 
    mv deps deps_cached
fi
mkdir -p deps_cached

# Get the correct code
git submodule update --init --recursive

# Clobber any differing files
rsync -r --links --checksum deps/ deps_cached/

rm -Rf deps

# And move the built dependencies into place
mv deps_cached deps


#!/usr/bin/env bash
# Build vg at one commit in a fresh scratch directory, and print the binary's path last.
#
#   build_vg.sh <commit> [vg repository, default the one this script is in]
#
# The build goes on local disk, because NFS and /private/groups are slow for a source
# tree. It goes in a new directory from mktemp under $TMPDIR or /tmp, so it cannot collide
# with another user's directory. Delete the directory when done.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
COMMIT=${1:?commit}
REPO=${2:-$(git -C "$HERE" rev-parse --show-toplevel)}
JOBS=${JOBS:-64}

BUILD=$(mktemp -d "${TMPDIR:-/tmp}/vgbuild.XXXXXX")
echo "building $COMMIT in $BUILD" >&2
git clone -q "$REPO" "$BUILD"
cd "$BUILD"
git checkout -q "$COMMIT"
# Fetch submodules from GitHub, not from the clone's origin, which is a local directory.
git remote set-url origin https://github.com/vgteam/vg.git
git submodule update --init --recursive --jobs 16 -q
# Commits before caa925e81, such as 202fc35d3, use std::unordered_map in
# linkage_model.hpp without including it, which the current compiler rejects. Including
# it everywhere changes nothing else.
CXXFLAGS="-include unordered_map" make -j"$JOBS" > build.log 2>&1 \
    || { tail -n 30 build.log >&2; echo "build failed; log in $BUILD/build.log" >&2; exit 1; }
# vg links some libraries, such as pixman, from its own lib/ without an rpath, so it runs
# only where the system has a copy. Name any missing library now, rather than at run time.
if ldd bin/vg | grep -q 'not found'; then
    ldd bin/vg | grep 'not found' >&2
    echo "bin/vg needs the libraries above; install them (for pixman: libpixman-1-0)" >&2
    exit 1
fi
echo "$BUILD/bin/vg"

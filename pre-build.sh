#!/usr/bin/env bash
# pre-build.sh: setup commands to run before building anything.
# Invoked by the Makefile automatically.
# If this script fails, the build will fail.
set -e

# TODO: quitting if no protoc doesn't reliably stop the build.
protoc --version >/dev/null 2>/dev/null || (echo "Error: protobuf compiler (protoc) not available!" ; exit 1)
if [ -e include/vg/vg.pb.h ] ; then
    HEADER_VER=$(cat include/vg/vg.pb.h | grep GOOGLE_PROTOBUF_VERSION | sed 's/[^0-9]*\([0-9]*\)[^0-9]*/\1/' | head -n1); \
    WORKDIR=$(pwd); \
    TESTDIR=$(mktemp -d); \
    echo 'syntax = "proto3";' > ${TESTDIR}/empty.proto; \
    protoc ${TESTDIR}/empty.proto --proto_path=${TESTDIR} --cpp_out=${TESTDIR}; \
    PROTOC_VER=$(cat ${TESTDIR}/empty.pb.h | grep GOOGLE_PROTOBUF_VERSION | sed 's/[^0-9]*\([0-9]*\)[^0-9]*/\1/' | head -n1); \
    if [ "${HEADER_VER}" != "${PROTOC_VER}" ] ; then \
        echo "Protobuf version has changed!"; \
        echo "Headers are for ${HEADER_VER} but we make headers for ${PROTOC_VER}"; \
        echo "Need to rebuild libvgio"; \
        rm -f lib/libvgio.a; \
        rm -f include/vg/vg.pb.h; \
    fi; \
    rm ${TESTDIR}/empty.proto ${TESTDIR}/empty.pb.h ${TESTDIR}/empty.pb.cc; \
    rmdir ${TESTDIR}; \
fi;

# A note about Protobuf:
# We have a lot of logic here to make sure that the protoc we have generates headers with exactly the same
# version requirements as the headers we already have.
# If not, we regenerate them.
# Doesn't handle Protobuf 3.12.3 weirdness; just make clean if you change flavors of Protobuf 3.12.3.

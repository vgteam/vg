# Multi-container Dockerfile for build and run containers for vg
FROM ubuntu:18.04 AS base
MAINTAINER vgteam

# Install shared build and run dependencies
RUN apt-get -qq update && apt-get -qq install -y curl wget pigz dstat pv jq samtools tabix parallel fontconfig-config
# One dependency is bwa, for vg tests.
COPY --from=quay.io/ucsc_cgl/bwa:0.7.15--a17c6544342330f6ea7a23a37d23273ab1c52d21 /usr/local/bin/bwa /usr/local/bin/bwa

FROM base AS build
# install build apt dependencies
# note: most vg apt dependencies are installed by "make get-deps" below
RUN apt-get -qq update && apt-get -qq install -y \
    sudo \
    bsdmainutils \
    build-essential \
    make \
    git \
    zlib1g-dev \
    rs \
    gdb \
    time \
    gawk

# Copy vg build tree into place
ADD . /vg
WORKDIR /vg

# To increase portability of the docker image, set the target CPU architecture to
# Ivy Bridge (2012) rather than auto-detecting the build machine's CPU.
# Since we build and run protoc, we need ivybridge or later to actually build the container.
# But we then don't depend on the build host's architecture
RUN sed -i s/march=native/march=ivybridge/ deps/sdsl-lite/CMakeLists.txt
RUN . ./source_me.sh && make get-deps && CXXFLAGS=" -march=ivybridge " make -j$(nproc) && make static

ENV PATH /vg/bin:$PATH

FROM base AS run

COPY --from=build /vg/bin/vg /vg/bin/vg
COPY --from=build /vg/scripts /vg/scripts

WORKDIR /vg
ENV PATH /vg/bin:$PATH

FROM build AS test
RUN make test


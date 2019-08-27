# Multi-container Dockerfile for build and run containers for vg
FROM ubuntu:18.04 AS base
MAINTAINER vgteam

FROM base AS build

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

############################################################################################
FROM build AS test

# The test need BWA
COPY --from=quay.io/ucsc_cgl/bwa:0.7.15--a17c6544342330f6ea7a23a37d23273ab1c52d21 /usr/local/bin/bwa /usr/local/bin/bwa

# The tests need some extra packages.
# TODO: Which of these can we remove?
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    pigz \
    dstat \
    pv \
    jq \
    samtools \
    tabix \
    parallel \
    fontconfig-config \
    && apt-get -qq -y clean

# Fail if any non-portable instructions were used
RUN /bin/bash -e -c 'if objdump -d /vg/bin/vg | grep vperm2i128 ; then exit 1 ; else exit 0 ; fi'
# Run tests in the middle so the final container that gets tagged is the run container.
RUN make test


############################################################################################
FROM base AS run

COPY --from=build /vg/bin/vg /vg/bin/vg
COPY --from=build /vg/scripts /vg/scripts

# Install packages which toil-vg needs to be available inside the image, for pipes
# TODO: which of these can be removed?
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    curl \
    wget \
    pigz \
    dstat \
    pv \
    jq \
    samtools \
    tabix \
    parallel \
    fontconfig-config \
    && apt-get -qq -y clean

WORKDIR /vg
ENV PATH /vg/bin:$PATH




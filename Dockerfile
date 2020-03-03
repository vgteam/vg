# Multi-container Dockerfile for build and run containers for vg
FROM ubuntu:18.04 AS base
MAINTAINER vgteam

RUN echo base > /stage.txt

WORKDIR /vg

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

FROM base AS build

RUN echo build > /stage.txt

# Copy vg build tree into place
COPY . /vg

# Install the base packages needed to let vg install packages.
# Make sure this runs after vg sources are imported so vg will always have an
# up to date package index to get its dependencies.
# We don't need to clean the package index since we don't ship this image and
# don't care about its size.
# We don't want to install too much stuff here, because we want to test vg's
# make get-deps to make sure it isn't missing something
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    make \
    sudo \
    git

# If we're trying to build from a non-recursively-cloned repo, go get the
# submodules.
RUN bash -c "[[ -e deps/sdsl-lite/CMakeLists.txt ]] || git submodule update --init --recursive"

# To increase portability of the docker image, set the target CPU architecture to
# Nehalem (2008) rather than auto-detecting the build machine's CPU.
# This has no AVX1, AVX2, or PCLMUL, but it does have SSE4.2.
# UCSC has a Nehalem machine that we want to support.
RUN sed -i s/march=native/march=nehalem/ deps/sdsl-lite/CMakeLists.txt
RUN make get-deps && . ./source_me.sh && env && make include/vg_git_version.hpp && CXXFLAGS=" -march=nehalem " make -j8 && make static && strip bin/vg

ENV PATH /vg/bin:$PATH

############################################################################################
FROM build AS test

RUN echo test > /stage.txt

# The tests also need some other extra packages.
# TODO: Which of these can we remove?
# No clean necessary since we aren't shipping this
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    bwa \
    pigz \
    dstat \
    pv \
    jq \
    samtools \
    tabix \
    parallel \
    bsdmainutils \
    rs \
    fontconfig-config

# Fail if any non-portable instructions were used
RUN /bin/bash -e -c 'if objdump -d /vg/bin/vg | grep vperm2i128 ; then exit 1 ; else exit 0 ; fi'
# Run tests in the middle so the final container that gets tagged is the run container.
RUN make test


############################################################################################
FROM base AS run

RUN echo run > /stage.txt

COPY --from=build /vg/bin/vg /vg/bin/

COPY --from=build /vg/scripts/* /vg/scripts/
# Make sure we have the flame graph scripts so we can do self-profiling
COPY deps/FlameGraph /vg/deps/FlameGraph

# Install packages which toil-vg needs to be available inside the image, for
# pipes and profiling, and good usability on Kubernetes.
# TODO: which of these can be removed?
# Make sure to clean so we don't ship old apt package indexes in our Docker.
RUN ls -lah /vg && \
    apt-get -qq -y update && \
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
    awscli \
    binutils \
    libssl1.0.0 \
    libpython2.7 \
    libperl-dev \
    libelf1 \
    libdw1 \
    libslang2 \
    libnuma1 \
    numactl \
    bc \
    linux-tools-common \
    linux-tools-generic \
    perl \
    && apt-get -qq -y clean


ENV PATH /vg/bin:$PATH




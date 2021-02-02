# Multi-container Dockerfile for build and run containers for vg
FROM ubuntu:18.04 AS base
MAINTAINER vgteam

RUN echo base > /stage.txt

WORKDIR /vg

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

FROM base AS build
ARG THREADS=8

RUN echo build > /stage.txt

RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install sudo

# Install all vg's dependencies.
# The Makefile will come parse the Dockerfile to get the correct dependencies;
# this is the One True Depencency List.
# We don't need to clean the package index since we don't ship this image and
# don't care about its size.
# We clip out everything between these begin and end markers, except the line
# that starts with RUN, or comments. And we pull out line continuation slashes.
# TODO: can we read them here and in the Makefile from the README instead?
###DEPS_BEGIN###
RUN apt-get -qq -y update && apt-get -qq -y upgrade && apt-get -qq -y install \
    make git build-essential protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev \
    libncurses5-dev automake libtool jq bsdmainutils bc rs parallel npm samtools curl \
    unzip redland-utils librdf-dev cmake pkg-config wget gtk-doc-tools raptor2-utils \
    rasqal-utils bison flex gawk libgoogle-perftools-dev liblz4-dev liblzma-dev \
    libcairo2-dev libpixman-1-dev libffi-dev libcairo-dev libprotobuf-dev libboost-all-dev
###DEPS_END###


# Copy vg build tree into place
COPY . /vg

# Updater again after vg sources are imported so vg will always have an
# up to date package index to get its dependencies.
RUN apt-get -qq -y update

# If we're trying to build from a non-recursively-cloned repo, go get the
# submodules.
RUN bash -c "[[ -e deps/sdsl-lite/CMakeLists.txt ]] || git submodule update --init --recursive"

# To increase portability of the docker image, set the target CPU architecture to
# Nehalem (2008) rather than auto-detecting the build machine's CPU.
# This has no AVX1, AVX2, or PCLMUL, but it does have SSE4.2.
# UCSC has a Nehalem machine that we want to support.
RUN sed -i s/march=native/march=nehalem/ deps/sdsl-lite/CMakeLists.txt
# Do the build. Trim down the resulting binary but make sure to include enough debug info for profiling.
RUN . ./source_me.sh && make include/vg_git_version.hpp && CXXFLAGS=" -march=nehalem " make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) static && strip -d bin/vg

ENV PATH /vg/bin:$PATH

############################################################################################
FROM build AS test
ARG THREADS=8

RUN echo test > /stage.txt

# Fail if any non-portable instructions were used
RUN /bin/bash -e -c 'if objdump -d /vg/bin/vg | grep vperm2i128 ; then exit 1 ; else exit 0 ; fi'
# Run tests in the middle so the final container that gets tagged is the run container.
RUN export OMP_NUM_THREADS=$((THREADS < $(nproc) ? THREADS : $(nproc))) make test


############################################################################################
FROM base AS run

RUN echo run > /stage.txt

# Install packages which toil-vg needs to be available inside the image, for
# pipes and profiling, and good usability on Kubernetes.
# TODO: which of these can be removed?
# Make sure to clean so we don't ship old apt package indexes in our Docker.
RUN ls -lah /vg && \
    apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install --no-upgrade \
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
    
COPY --from=build /vg/bin/vg /vg/bin/

COPY --from=build /vg/scripts/* /vg/scripts/
# Make sure we have the flame graph scripts so we can do self-profiling
COPY deps/FlameGraph /vg/deps/FlameGraph

ENV PATH /vg/bin:$PATH




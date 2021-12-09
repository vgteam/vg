# Multi-container Dockerfile for build and run containers for vg

# Use Google's non-rate-limited mirror of Docker Hub to get our base image.
# This helps automated Quay builds because Quay hasn't built a caching system
# and exposes pull rate limits to users.
FROM mirror.gcr.io/library/ubuntu:20.04 AS base
MAINTAINER vgteam

RUN echo base > /stage.txt

WORKDIR /vg

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

FROM base AS build
ARG THREADS=8
ARG TARGETARCH

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
    libcairo2-dev libpixman-1-dev libffi-dev libcairo-dev libprotobuf-dev libboost-all-dev \
    tabix bcftools
###DEPS_END###

# Prepare to build submodule dependencies
COPY source_me.sh /vg/source_me.sh
COPY deps /vg/deps
# To increase portability of the docker image, when building for amd64, set the
# target CPU architecture to Nehalem (2008) rather than auto-detecting the
# build machine's CPU. This has no AVX1, AVX2, or PCLMUL, but it does have
# SSE4.2. UCSC has a Nehalem machine that we want to support.
RUN if [ -z "${TARGETARCH}" ] || [ "${TARGETARCH}" = "amd64" ] ; then sed -i s/march=native/march=nehalem/ deps/sdsl-lite/CMakeLists.txt; fi
# Clear any CMake caches in case we are building from someone's checkout
RUN find . -name CMakeCache.txt | xargs rm -f
# Build the dependencies
COPY Makefile /vg/Makefile
RUN . ./source_me.sh && CXXFLAGS="$(if [ -z "${TARGETARCH}" ] || [ "${TARGETARCH}" = "amd64" ] ; then echo " -march=nehalem "; fi)" make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) deps

# Bring in the sources, which we need in order to build
COPY src /vg/src

# Build all the object files for vg, but don't link.
# Also pass the arch here
RUN . ./source_me.sh && CXXFLAGS="$(if [ -z "${TARGETARCH}" ] || [ "${TARGETARCH}" = "amd64" ] ; then echo " -march=nehalem "; fi)" make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) objs

# Bring in any includes we pre-made, like the git version, if present
COPY include /vg/include

# Make sure version introspection is up to date
RUN rm -f obj/version.o && . ./source_me.sh && CXXFLAGS="$(if [ -z "${TARGETARCH}" ] || [ "${TARGETARCH}" = "amd64" ] ; then echo " -march=nehalem "; fi)" make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) obj/version.o

# Announce the version file, which must exist by now
RUN ls /vg/include && cat /vg/include/vg_git_version.hpp

# Do the final build and link, knowing the version. Trim down the resulting binary but make sure to include enough debug info for profiling.
RUN . ./source_me.sh && CXXFLAGS="$(if [ -z "${TARGETARCH}" ] || [ "${TARGETARCH}" = "amd64" ] ; then echo " -march=nehalem "; fi)" make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) static && strip -d bin/vg

# Ship the scripts
COPY scripts /vg/scripts

ENV PATH /vg/bin:$PATH

############################################################################################
FROM build AS test
ARG THREADS=8

RUN echo test > /stage.txt

RUN curl -sL https://deb.nodesource.com/setup_16.x | bash - && apt-get -qq -y install nodejs && npm install -g txm@7.4.5

# Fail if any non-portable instructions were used
RUN /bin/bash -e -c 'if objdump -d /vg/bin/vg | grep vperm2i128 ; then exit 1 ; else exit 0 ; fi'

# Bring in the tests and docs, which have doctests
COPY test /vg/test
COPY doc /vg/doc
# We test the README so bring it along.
COPY README.md /vg/

# Run tests in the middle so the final container that gets tagged is the run container.
# Tests may not actually be run by smart builders like buildkit.
RUN /bin/bash -e -c "export OMP_NUM_THREADS=$((THREADS < $(nproc) ? THREADS : $(nproc))); make test"


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
COPY --from=build /vg/deps/FlameGraph /vg/deps/FlameGraph

ENV PATH /vg/bin:$PATH




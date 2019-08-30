# Multi-container Dockerfile for build and run containers for vg
FROM ubuntu:18.04 AS base
MAINTAINER vgteam

RUN echo base > /stage.txt

WORKDIR /vg

RUN ls -lah /vg || echo "No vg directory exists yet"

FROM base AS build

RUN echo build > /stage.txt

RUN ls -lah /vg || echo "No vg directory exists yet"

# Copy vg build tree into place
COPY . /vg

RUN ls -lah /vg || echo "No vg directory exists yet"

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
    sudo

# To increase portability of the docker image, set the target CPU architecture to
# Ivy Bridge (2012) rather than auto-detecting the build machine's CPU.
# Since we build and run protoc, we need ivybridge or later to actually build the container.
# But we then don't depend on the build host's architecture
RUN sed -i s/march=native/march=ivybridge/ deps/sdsl-lite/CMakeLists.txt
RUN make get-deps && . ./source_me.sh && env && make include/vg_git_version.hpp && CXXFLAGS=" -march=ivybridge " make -j$(nproc) && make static && strip bin/vg

ENV PATH /vg/bin:$PATH

############################################################################################
FROM build AS test

RUN echo test > /stage.txt

# The test need BWA
COPY --from=quay.io/ucsc_cgl/bwa:0.7.15--a17c6544342330f6ea7a23a37d23273ab1c52d21 /usr/local/bin/bwa /usr/local/bin/bwa

# The tests need some extra packages.
# TODO: Which of these can we remove?
# No clean necessary since we aren't shipping this
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

RUN ls -lah /vg || echo "No vg directory exists yet"

COPY --from=build /vg/bin/vg /vg/bin/

RUN ls -lah /vg || echo "No vg directory exists yet"

COPY --from=build /vg/scripts/* /vg/scripts/

RUN ls -lah /vg || echo "No vg directory exists yet"

# Install packages which toil-vg needs to be available inside the image, for pipes
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
    && apt-get -qq -y clean


ENV PATH /vg/bin:$PATH




# Dockerfile for a full vg build from source
# (derived from vgteam/vg_docker)

FROM ubuntu:18.04
MAINTAINER vgteam
ARG vg_git_revision=master

# Make sure the en_US.UTF-8 locale exists, since we need it for tests
#RUN locale-gen en_US en_US.UTF-8 && DEBIAN_FRONTEND=noninteractive dpkg-reconfigure locales

# install basic apt dependencies
# note: most vg apt dependencies are installed by "make get-deps" below
RUN apt-get -qq update && apt-get -qq install -y \
    sudo \
    pv \
    pigz \
    bsdmainutils \
    build-essential \
    make \
    git \
    zlib1g-dev \
    rs \
    libffi-dev

ADD https://github.com/vgteam/vg_docker/raw/master/deps/bwa_0.7.15-5_amd64.deb /tmp/bwa.deb
RUN dpkg -i /tmp/bwa.deb

# copy over current directory to docker
ADD . /vg

# set our working directory
WORKDIR /vg

# Build
RUN . ./source_me.sh && make get-deps && make -j$(nproc) && make static

ENV PATH /vg/bin:$PATH

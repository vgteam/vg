FROM ubuntu:16.04

MAINTAINER Erik Garrison <erik.garrison@gmail.com>

# Make sure the en_US.UTF-8 locale exists, since we need it for tests
RUN locale-gen en_US en_US.UTF-8 && dpkg-reconfigure locales

# Set up for make get-deps
RUN mkdir /app
WORKDIR /app
COPY Makefile /app/Makefile

RUN sed -i "s/sudo//g" /app/Makefile

# Install vg dependencies and clear the package index
RUN \
    echo "deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse" | tee -a /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
        build-essential \
        libprotoc-dev \
        gcc-5-base \
        libgcc-5-dev \
        pkg-config \
        jq/trusty-backports \
        zlib1g-dev \
        && \
    make get-deps
    
# Move in all the other files
COPY . /app 

ENV LD_LIBRARY_PATH=/app/lib
ENV LIBRARY_PATH /app/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH /app/lib:$LD_LIBRARY_PATH
ENV LD_INCLUDE_PATH /app/include:$LD_INCLUDE_PATH
ENV C_INCLUDE_PATH /app/include:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH /app/include:$CPLUS_INCLUDE_PATH
ENV INCLUDE_PATH /app/include:$INCLUDE_PATH

RUN cd /app && . ./source_me.sh && make && make test
ENTRYPOINT ["/app/bin/vg"]


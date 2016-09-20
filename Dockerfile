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
        gcc-5-base \
        libgcc-5-dev \
        pkg-config \
        jq/trusty-backports \
        && \
    make get-deps && \
    rm -rf /var/lib/apt/lists/*
    
# Move in all the other files
COPY . /app
    
# Build vg
RUN cd /app && . ./source_me.sh && make -j8

# Make tests. We can't do it in parallel since it cleans up the test binary
RUN cd /app && . ./source_me.sh make test

ENV LD_LIBRARY_PATH=/app/lib
ENV LIBRARY_PATH /app/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH /app/lib:$LD_LIBRARY_PATH
ENV LD_INCLUDE_PATH /app/include:$LD_INCLUDE_PATH
ENV C_INCLUDE_PATH /app/include:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH /app/include:$CPLUS_INCLUDE_PATH
ENV INCLUDE_PATH /app/include:$INCLUDE_PATH

#ENTRYPOINT ["/app/bin/vg"]
RUN cp /app/bin/vg /usr/bin


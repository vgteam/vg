FROM ubuntu:14.04

MAINTAINER Erik Garrison <erik.garrison@gmail.com>

# Make sure the en_US.UTF-8 locale exists, since we need it for tests
RUN locale-gen en_US en_US.UTF-8 && dpkg-reconfigure locales

# Set up for make get-deps
RUN mkdir /app
WORKDIR /app
COPY Makefile /app/Makefile

# Install vg dependencies and clear the package index
RUN \
    apt-get update && \
    apt-get install -y \
        build-essential \
        pkg-config
        sudo && \
    make get-deps && \
    rm -rf /var/lib/apt/lists/*
    
# Move in all the other files
COPY . /app
    
# Build vg
RUN make -j8

# Make tests. We can't do it in parallel since it cleans up the test binary
RUN make test

ENTRYPOINT ["/app/vg"]


#!/usr/bin/env bash
#
# giraffe-speed.sh: evaluate Giraffe mapping speed.
#
# Although CPU instruction counting doesn't work in this configuration for some
# reason, you can run on an AWS r5.4xlarge, with 100 GB disk, the
# "vg-data-user" IAM role and Ubuntu 20.04:
#
# git clone https://github.com/vgteam/vg.git
# cd vg
# git checkout <commit to test>
# git submodule update --init --recursive
# sudo apt update
# sudo apt install -y build-essential awscli
# make get-deps
# make -j16
# ./scripts/giraffe-speed.sh
#
# In this configuration, a mapping speed of 3520.61 reads per second per
# thread, or 3527.4 reads per CPU-second (including output), is typical.
# Significant reductions may indicate a performance regression.
#
# Or, run in Docker or on Kubernetes using the Docker build of the version of
# vg you want to evaluate. Allocate at least 50 GB memory, 16 hyperthreads, and
# 100 GB of disk. The script dependencies should already be available.
#
# In this configuration, speed will vary with the machine type, and the load on
# the host, but a mapping speed of 3836.3 reads per CPU-second (including
# output) has been observed in a heavily loaded Kubernetes pod, and 5333.43
# reads per CPU-second (including output) in a lightly loaded one, on the GI
# Kubernetes cluster. In all cases, 0.918424 M instructions per read is a
# typical value for GI Kubernetes's CPUs; significant changes from that
# indicate that more or less work is now required to map each read.
#
# We expect to be root, or to be able to sudo, so Giraffe can evaluate its
# instruction execution speed.
#
# Partially based on
# https://github.com/vgteam/giraffe-sv-paper/blob/566573b708878d8854acb088a0a8f7c920b120eb/scripts/giraffe/speed_giraffe.sh
set -e

THREAD_COUNT=16

READS=s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz

# The HGSVC graph is smaller and faster to load so we use that here.
GRAPH=hgsvc
GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
# sampled.64 isn't much bigger than the base GBWT for HGSVC
GBWT="sampled.64"

function fetch() {
    # Download a URL to a file, if it isn't there already
    if [ ! -e "${2}" ] ; then
        # For now we only do S3 URLs
        aws s3 cp "${1}" "${2}"
    fi
}

fetch "${READS}" "./reads.fq.gz"
fetch "${GRAPH_BASE}.xg" "./${GRAPH}.xg"
fetch "${GRAPH_BASE}.dist" "./${GRAPH}.dist"
fetch "${GRAPH_BASE}.${GBWT}.gbwt" "./${GRAPH}.${GBWT}.gbwt"
fetch "${GRAPH_BASE}.${GBWT}.gg" "./${GRAPH}.${GBWT}.gg"
fetch "${GRAPH_BASE}.${GBWT}.min" "./${GRAPH}.${GBWT}.min"

# Build a bigger reads file (10m) so we can run for more than like 16 seconds
cat ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz ./reads.fq.gz >./many-reads.fq.gz

SUDO=$(which sudo 2>/dev/null || true)
VG=$(which vg 2>/dev/null || echo "bin/vg")

timeout -k1 1h bash -c "${SUDO} ${VG} giraffe -x ${GRAPH}.xg -H ${GRAPH}.${GBWT}.gbwt -g ${GRAPH}.${GBWT}.gg -m ${GRAPH}.${GBWT}.min -d ${GRAPH}.dist -f ./many-reads.fq.gz -i -t ${THREAD_COUNT} -p 2>log.txt >mapped.gam" || true

cat log.txt


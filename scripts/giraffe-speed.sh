#!/usr/bin/env bash
# giraffe-speed.sh: evaluate Giraffe mapping speed.
# Suggested to run on an AWS
# Partially based on https://github.com/vgteam/giraffe-sv-paper/blob/566573b708878d8854acb088a0a8f7c920b120eb/scripts/giraffe/speed_giraffe.sh
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

${SUDO} timeout -k1 1h bash -c "vg giraffe -x '${GRAPH}.xg' -H '${GRAPH}.${GBWT}.gbwt' -g '${GRAPH}.${GBWT}.gg' -m '${GRAPH}.${GBWT}.min' -d '${GRAPH}.dist' -f './many-reads.fq.gz' -i -t '${THREAD_COUNT}' -p 2>log.txt >mapped.gam" || true


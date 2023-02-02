#!/usr/bin/env bash
# lr_benchmark.sh: Run a benchmark for vg long read mapping
# Meant to be run on UCSC Courtyard/Plaza

set -e
set -o pipefail

# Here we use : and := to set variables to default values if not present in the environment.
# You can set these in the environment to override them and I don't have to write a CLI option parser.
# See https://stackoverflow.com/a/28085062

# Where should output go?
: "${OUT_DIR:="./lr_benchmark"}"
echo "Writing to ${OUT_DIR}"
mkdir -p "${OUT_DIR}"

# Adam Novak's simulated reads, loosely following Stephen Hwang's method.
# Annotated with GRCh38.chr1 style path names.
# Also available in 100 and 10000 read versions
# On GRCh38, sample HG00741:
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/r9.5-acc0.95/HG00741/HG00741-sim-r9.5-acc0.95-1000.gam
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/hifi/HG00741/HG00741-sim-hifi-1000.gam
# On CHM13, sample HG002:
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/r9.5-acc0.95/HG002/HG002-sim-r9.5-acc0.95-1000.gam
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/hifi/HG002/HG002-sim-hifi-1000.gam
# Note that not all of these reads have true positions on CHM13 due to a lack of alignment for some HG002 regions!
: "${INPUT_READ_PATH:=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/r9.5-acc0.95/HG002/HG002-sim-r9.5-acc0.95-1000.gam}"

# For GRCh38 mapping:
# An HPRC graph, linked to /public/groups/cgl/graph-genomes/xhchang/hprc_graph/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.giraffe.gbz
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.giraffe.gbz
# Its indexes
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.dist
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.min
# /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.xg

# For CHM13 mapping:
# An HPRC graph, with .min, .dist, .gbwt, and .xg files at least
: "${INPUT_GRAPH_BASE_URL:=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13}"
# Where it is/will be stored on disk locally, along with a .gbz
: "${INPUT_GRAPH_BASE_PATH:=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/hprc-v1.0-mc-chm13}"

# Find the indexes
: "${INPUT_DIST_PATH:=${INPUT_GRAPH_BASE_PATH}.dist}"
: "${INPUT_MIN_PATH:=${INPUT_GRAPH_BASE_PATH}.min}"
: "${INPUT_XG_PATH:=${INPUT_GRAPH_BASE_PATH}.xg}"
: "${INPUT_GBZ_PATH:=${INPUT_GRAPH_BASE_PATH}.gbz}"

# Download the files that are probably there.
if [[ ! -e "${INPUT_DIST_PATH}" ]] ; then
    wget "${INPUT_GRAPH_BASE_URL}.dist" -O "${INPUT_DIST_PATH}.tmp"
    mv "${INPUT_DIST_PATH}.tmp" "${INPUT_DIST_PATH}"
fi
if [[ ! -e "${INPUT_MIN_PATH}" ]] ; then
    wget "${INPUT_GRAPH_BASE_URL}.min" -O "${INPUT_MIN_PATH}.tmp"
    mv "${INPUT_MIN_PATH}.tmp" "${INPUT_MIN_PATH}"
fi
if [[ ! -e "${INPUT_XG_PATH}" ]] ; then
    wget "${INPUT_GRAPH_BASE_URL}.xg" -O "${INPUT_XG_PATH}.tmp"
    mv "${INPUT_XG_PATH}.tmp" "${INPUT_XG_PATH}"
fi

if [[ ! -e "${INPUT_GBZ_PATH}" ]] ; then
    # Make GBZ from other files if not there
    
    # For which we need the GBWT
    : "${INPUT_GBWT_PATH:=${INPUT_GRAPH_BASE_PATH}.gbwt}"
    if [[ ! -e "${INPUT_GBWT_PATH}" ]] ; then
        wget "${INPUT_GRAPH_BASE_URL}.gbwt" -O "${INPUT_GBWT_PATH}.tmp"
        mv "${INPUT_GBWT_PATH}.tmp" "${INPUT_GBWT_PATH}"
    fi
    
    time vg gbwt -x "${INPUT_XG_PATH}" "${INPUT_GBWT_PATH}" --gbz-format -g "${INPUT_GBZ_PATH}.tmp"
    mv "${INPUT_GBZ_PATH}.tmp" "${INPUT_GBZ_PATH}"
fi

if [[ "${WORK_DIR}" == "" ]] ; then
    # Make a work directory
    WORK_DIR="$(mktemp -d)"
    CLEAN_WORK_DIR=1
else
    # Let the user send one in in the environment.
    mkdir -p "${WORK_DIR}"
    CLEAN_WORK_DIR=0
fi

echo "Working in ${WORK_DIR}"

if [[ ! -e "${WORK_DIR}/possible.txt" ]] ; then
    # Get the list of all reads with ref positions set, and which are thus possible to get right
    vg view -aj "${INPUT_READ_PATH}" | jq -r 'select(.refpos) | .name' > "${WORK_DIR}/possible.txt.tmp"
    mv "${WORK_DIR}/possible.txt.tmp" "${WORK_DIR}/possible.txt"
fi

if [[ ! -e "${WORK_DIR}/annotated.gam" ]] ; then
    # Map reads using correctness tracking.
    # Make sure to apply multi-position annotation which Giraffe won't do.
    vg giraffe -G "${INPUT_READ_PATH}" -t 16 -B 8 --align-from-chains -Z "${INPUT_GBZ_PATH}" -d "${INPUT_DIST_PATH}" -m "${INPUT_MIN_PATH}" -x "${INPUT_XG_PATH}" --track-provenance --track-correctness --progress | vg annotate -x "${INPUT_XG_PATH}" -a - --multi-position -l 100 >"${WORK_DIR}/annotated.gam.tmp"
    mv "${WORK_DIR}/annotated.gam.tmp" "${WORK_DIR}/annotated.gam"
fi

# Compute general stats
vg stats -a "${WORK_DIR}/annotated.gam" >"${OUT_DIR}/stats.txt"

if [[ ! -e "${WORK_DIR}/benchmark.tsv" ]] ; then
    # See if reads get close enough to be correct
    # TODO: vg gamcompare announces a correctness count, which we should save
    vg gamcompare -r 200 "${WORK_DIR}/annotated.gam" "${INPUT_READ_PATH}" --aligner lrgiraffe --tsv >"${WORK_DIR}/benchmark.tsv.tmp"
    mv "${WORK_DIR}/benchmark.tsv.tmp" "${WORK_DIR}/benchmark.tsv"
fi

# Make a QQ plot
scripts/plot-qq.R "${WORK_DIR}/benchmark.tsv" "${OUT_DIR}/qq.png"

# Compute a correctness rate
TOTAL_READS="$(cat "${WORK_DIR}/benchmark.tsv" | tail -n +2 | wc -l)"
POSSIBLE_READS="$(cat "${WORK_DIR}/possible.txt" | wc -l)"
CORRECT_READS="$(cat "${WORK_DIR}/benchmark.tsv" | tail -n +2 | grep "^1" | wc -l)"
CORRECT_FRACTION_TOTAL="$(echo "${CORRECT_READS}/${TOTAL_READS}" | bc -l)"
CORRECT_FRACTION_POSSIBLE="$(echo "${CORRECT_READS}/${POSSIBLE_READS}" | bc -l)"
echo "Correct reads: ${CORRECT_READS}" >"${OUT_DIR}/results.txt"
echo "Total reads: ${TOTAL_READS}" >>"${OUT_DIR}/results.txt"
echo "Reads with truth positions: ${POSSIBLE_READS}" >>"${OUT_DIR}/results.txt"
echo "Correct fraction: ${CORRECT_FRACTION_TOTAL} of all reads, ${CORRECT_FRACTION_POSSIBLE} of reads with truth positions" >>"${OUT_DIR}/results.txt"
cat "${OUT_DIR}/results.txt"

if [[ "${CLEAN_WORK_DIR}" == "1" ]] ; then
    # Clean up the work directory
    rm -Rf "${WORK_DIR}"
fi


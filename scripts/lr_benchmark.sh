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
# Also available: /public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/r9.5-acc0.95/HG00741/HG00741-sim-r9.5-acc0.95-1000.gam
: "${INPUT_READ_PATH:=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim/hifi/HG00741/HG00741-sim-hifi-1000.gam}"

# An HPRC graph, linked to /public/groups/cgl/graph-genomes/xhchang/hprc_graph/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.giraffe.gbz
: "${INPUT_GBZ_PATH:=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.giraffe.gbz}"
# Its indexes
: "${INPUT_DIST_PATH:=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.dist}"
: "${INPUT_MIN_PATH:=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.min}"
: "${INPUT_XG_PATH=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.xg}"

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

if [[ ! -e "${WORK_DIR}/annotated.gam" ]] ; then
    # Map reads using correctness tracking.
    # Make sure to apply multi-position annotation which Giraffe won't do.
    vg giraffe -G "${INPUT_READ_PATH}" -t 16 -B 8 --align-from-chains -Z "${INPUT_GBZ_PATH}" -d "${INPUT_DIST_PATH}" -m "${INPUT_MIN_PATH}" -x "${INPUT_XG_PATH}" --track-provenance --track-correctness --progress | vg annotate -x "${INPUT_XG_PATH}" -a - --multi-position -l 100 >"${WORK_DIR}/annotated.gam"
fi

# Compute general stats
vg stats -a "${WORK_DIR}/annotated.gam" >"${OUT_DIR}/stats.txt"

if [[ ! -e "${WORK_DIR}/benchmark.tsv" ]] ; then
    # See if reads get close enough to be correct
    # TODO: vg gamcompare announces a correctness count, which we should save
    vg gamcompare -r 200 "${WORK_DIR}/annotated.gam" "${INPUT_READ_PATH}" --aligner lrgiraffe --tsv >"${WORK_DIR}/benchmark.tsv" 
fi

# Make a QQ plot
scripts/plot-qq.R "${WORK_DIR}/benchmark.tsv" "${OUT_DIR}/qq.png"

# Compute a correctness rate
TOTAL_READS="$(cat "${WORK_DIR}/benchmark.tsv" | tail -n +2 | wc -l)"
CORRECT_READS="$(cat "${WORK_DIR}/benchmark.tsv" | tail -n +2 | grep "^1" | wc -l)"
CORRECT_FRACTION="$(echo "${CORRECT_READS}/${TOTAL_READS}" | bc -l)"
echo "Correct fraction: ${CORRECT_FRACTION}" >"${OUT_DIR}/results.txt"
cat "${OUT_DIR}/results.txt"

if [[ "${CLEAN_WORK_DIR}" == "1" ]] ; then
    # Clean up the work directory
    rm -Rf "${WORK_DIR}"
fi


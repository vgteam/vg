#!/usr/bin/env bash
# lr_benchmark.sh: Run a benchmark for vg long read mapping
# Meant to be run on UCSC Courtyard/Plaza

set -e

# Where should output go?
OUT_DIR="./lr_benchmark"
echo "Writing to ${OUT_DIR}"
mkdir -p "${OUT_DIR}"

# Stephen Hwang's simulated reads, from /public/groups/vg/sjhwang/data/reads/sim_HiFi/sim_10k_reads.gam.
# Assumed to be annotated already.
INPUT_READ_PATH=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/reads/sim_HiFi/sim_10k_reads.gam
# An HPRC graph, linked to /public/groups/cgl/graph-genomes/xhchang/hprc_graph/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.giraffe.gbz
INPUT_GBZ_PATH=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.giraffe.gbz
# Its indexes
INPUT_DIST_PATH=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.dist
INPUT_MIN_PATH=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.min
INPUT_XG_PATH=/public/groups/cgl/graph-genomes/anovak/data/hprc-lrgiraffe/graphs/GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000.xg

# Make a work directory
WORKDIR="$(mktemp -d)"
echo "Working in ${WORKDIR}"

# Map reads using correctness tracking.
# Make sure to apply multi-position annotation which Giraffe won't do.
vg giraffe -G "${INPUT_READ_PATH}" -t 16 -B 8 --align-from-chains -Z "${INPUT_GBZ_PATH}" -d "${INPUT_DIST_PATH}" -m "${INPUT_MIN_PATH}" -x "${INPUT_XG_PATH}" --track-provenance --track_correctness --progress | vg annotate -x "${INPUT_XG_PATH}" -a - --multi-position >"${WORK_DIR}/annotated.gam"
# Compute general stats
vg stats -a "${WORK_DIR}/annotated.gam" >"${OUT_DIR}/stats.txt"
# See if they get close enough to be correct
vg gamcompare -r 200 "${WORK_DIR}/annotated.gam" "${INPUT_READ_PATH}" --aligner lrgiraffe --tsv "${OUT_DIR}/benchmark.tsv"
# Compute a correctness rate
TOTAL_READS="(cat "${OUT_DIR}/benchmark.tsv" | grep -v "^#" | wc -l)"
CORRECT_READS="(cat "${OUT_DIR}/benchmark.tsv" | grep -v "^#" | grep "^1" | wc -l)"
CORRECT_FRACTION="$(echo "${CORRECT_READS}/${TOTAL_READS}" | bc -l)"
echo "Correct fraction: ${CORRECT_FRACTION}" >"${OUT_DIR}/results.txt"
cat "${OUT_DIR}/results.txt"
# Make a QQ plot
scripts/plot-qq.R "${OUT_DIR}/benchmark.tsv" "${OUT_DIR}/qq.png"
# Clean up the work directory
rm -Rf "${WORKDIR}"


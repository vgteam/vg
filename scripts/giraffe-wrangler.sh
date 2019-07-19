#!/usr/bin/env bash

# giraffe-wrangler.sh: Run and profile vg gaffe and analyze the results.

set -e

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] FASTA XG_INDEX GBWT_INDEX MINIMIZER_INDEX DISTANCE_INDEX SIM_GAM REAL_FASTQ \n"
    printf "Options:\n\n"
    exit 1
}

while getopts "" o; do
    case "${o}" in
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "7" ]]; then
    # Too few arguments
    usage
fi

FASTA="${1}"
shift
XG_INDEX="${1}"
shift
GBWT_INDEX="${1}"
shift
MINIMIZER_INDEX="${1}"
shift
DISTANCE_INDEX="${1}"
shift
SIM_GAM="${1}"
shift
REAL_FASTQ="${1}"
shift

# Define the Giraffe parameters
GIRAFFE_OPTS=(-s75 -u 0.1 -v 25 -w 5 -C 600)
# And the thread count for everyone
THREAD_COUNT=32

# Define a work directory
# TODO: this requires GNU mptemp
WORK="$(mktemp -d)"

# Run simulated reads
vg gaffe -x "${XG_INDEX}" -m "${MINIMIZER_INDEX}" -H "${GBWT_INDEX}" -d "${DISTANCE_INDEX}" -G "${SIM_GAM}" -t "${THREAD_COUNT}" "${GIRAFFE_OPTS[@]}" >"${WORK}/mapped.gam"

# Annotate and compare against truth
vg annotate -p -x "${XG_INDEX}" -a "${WORK}/mapped.gam" >"${WORK}/annotated.gam"

# GAM compare against truth. Use gamcompare to count correct reads to save a JSON scan.
CORRECT_COUNT="$(vg gamcompare -r 100 "${WORK}/annotated.gam" "${SIM_GAM}" 2>&1 >/dev/null | sed 's/[^0-9]//g')"

# Compute identity
MEAN_IDENTITY="$(vg view -aj "${WORK}/mapped.gam" | jq -c '.identity' | awk '{x+=$1} END {print x/NR}')"

# TODO: Compute loss stages

# Now do the real reads

# Get RPS
vg gaffe -p -x "${XG_INDEX}" -m "${MINIMIZER_INDEX}" -H "${GBWT_INDEX}" -d "${DISTANCE_INDEX}" -f "${REAL_FASTQ}" -t "${THREAD_COUNT}" "${GIRAFFE_OPTS[@]}" >"${WORK}/real.gam" 2>"${WORK}/log.txt"

GIRAFFE_RPS="$(cat "${WORK}/log.txt" | grep "reads per second" | sed 's/[^0-9.]//g')"

# Get RPS for bwa-mem
REAL_READ_COUNT="$(cat "${REAL_FASTQ}" | wc -l)"
((REAL_READ_COUNT /= 4))

bwa mem -t "${THREAD_COUNT}" "${FASTA}" "${REAL_FASTQ}" >"${WORK}/mapped.bam" 2>"${WORK}/bwa-log.txt"
cat "${REAL_FASTQ}" "${REAL_FASTQ}" >"${WORK}/double.fq"
bwa mem -t "${THREAD_COUNT}" "${FASTA}" "${WORK}/double.fq" >"${WORK}/mapped-double.bam" 2>"${WORK}/bwa-log-double.txt"

BWA_TIME="$(cat "${WORK}/bwa-log.txt" | grep "Real time:" | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')"
BWA_DOUBLE_TIME="$(cat "${WORK}/bwa-log-double.txt" | grep "Real time:" | sed 's/[^0-9.]*\([0-9.]*\).*/\1/')"

BWA_RPS="$(echo "${REAL_READ_COUNT} / (${BWA_DOUBLE_TIME} - ${BWA_TIME})" | bc -l)"

# Print the report
echo "Giraffe got ${CORRECT_COUNT} simulated reads correct with ${MEAN_IDENTITY} average identity"
echo "Giraffe aligned real reads at ${GIRAFFE_RPS} reads/second vs. bwa-mem's ${BWA_RPS} reads/second"

rm -Rf "${WORK}"





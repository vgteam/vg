#!/usr/bin/env bash
# make_pbsim_reads.sh: script to simulate reads with pbsim2.
# Mostly theoretical; records commands that would have worked better than what was actually run
# Intended to run on UCSC Courtyard/Plaza systems
# You may also need to CFLAGS=-fPIC pip3 install --user bioconvert

set -ex

# Here we use : and := to set variables to default values if not present in the environment.
# You can set these in the environment ot override them and I don't have to write a CLI option parser.
# See https://stackoverflow.com/a/28085062

# Graph to simulate from
: "${GRAPH_XG_URL:=s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.xg}"
: "${GRAPH_GBWT_URL:=s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gbwt}"
# Name to use for graph when downloaded
: "${GRAPH_NAME:=hprc-v1.0-mc-grch38}"
# Sample to simulate from
: "${SAMPLE_NAME:=HG00741}"
# Technology name to use in output filenames
: "${TECH_NAME:=hifi}"
# FASTQ to use as a template
: "${SAMPLE_FASTQ:=/public/groups/vg/sjhwang/data/reads/real_HiFi/tmp/HiFi_reads_100k_real.fq}"
# This needs to be the pbsim2 command, which isn't assumed to be in $PATH
: "${PBSIM:=/public/groups/vg/sjhwang/tools/bin/pbsim}"
# Parameters to use with pbsim for simulating reads for each contig. Parameters are space-separated and internal spaces must be escaped.
: "${PBSIM_PARAMS:=--depth 1 --accuracy-min 0.00 --length-min 10000 --difference-ratio 6:50:54}"
# This needs to be a command line which can execute Stephen's script that adds qualities from a FASTQ back into a SAM that is missing them.
# Arguments are space-separated and internal spaces must be escaped.
# TODO: Put the script up somewhere so it can be run!
: "${ADD_QUALITIES:=python3 /public/groups/vg/sjhwang/vg_scripts/bin/readers/sam_reader.py}"
# Directory to save results in
: "${OUT_DIR:=./reads/sim/${TECH_NAME}/${SAMPLE_NAME}}"
# Number of MAFs to convert at once
: "${MAX_JOBS:=10}"

if [[ "${WORK_DIR}" == "" ]] ; then
    # Make a work directory
    WORK_DIR="$(mktemp -d)"
    CLEAN_WORK_DIR=1
else
    # Let the user send one in in the environment.
    CLEAN_WORK_DIR=0
fi


# Make sure scratch directory exists
mkdir -p "${WORK_DIR}"

# Fetch graph
if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.xg" ]] ; then
    aws s3 cp "${GRAPH_XG_URL}" "${WORK_DIR}/${GRAPH_NAME}.xg.tmp"
    mv "${WORK_DIR}/${GRAPH_NAME}.xg.tmp" "${WORK_DIR}/${GRAPH_NAME}.xg"
fi
if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.gbwt" ]] ; then
    aws s3 cp "${GRAPH_GBWT_URL}" "${WORK_DIR}/${GRAPH_NAME}.gbwt.tmp"
    mv "${WORK_DIR}/${GRAPH_NAME}.gbwt.tmp" "${WORK_DIR}/${GRAPH_NAME}.gbwt"
fi

if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.gbz" ]] ; then
    # Make it one file
    time vg gbwt -x "${WORK_DIR}/${GRAPH_NAME}.xg" "${WORK_DIR}/${GRAPH_NAME}.gbwt" --gbz-format -g "${WORK_DIR}/${GRAPH_NAME}.gbz.tmp"
    mv "${WORK_DIR}/${GRAPH_NAME}.gbz.tmp" "${WORK_DIR}/${GRAPH_NAME}.gbz"
fi

if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz" ]] ; then
    # Make it have our sample as the reference
    vg gbwt -Z trash/${GRAPH_NAME}.gbz --set-tag "reference_samples=${SAMPLE_NAME}" --gbz-format -g "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz.tmp"
    mv "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz.tmp" "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz"
fi

if [[ ! -e "${WORK_DIR}/${SAMPLE_NAME}.fa" ]] ; then
    # Extract sample assembly FASTA from graph where sample is the *reference*. If
    # we so it from the one where the sample is haplotypes, we get different path
    # name strings and we can't inject without hacking them up. We leave the code
    # to hack them up anyway though, for reference later.
    vg paths -x "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz" \
        --sample "${SAMPLE_NAME}" \
        --extract-fasta \
      > "${WORK_DIR}/${SAMPLE_NAME}.fa.tmp"
    mv "${WORK_DIR}/${SAMPLE_NAME}.fa.tmp" "${WORK_DIR}/${SAMPLE_NAME}.fa"
fi

if [[ ! -d "${WORK_DIR}/${SAMPLE_NAME}-reads" ]] ; then
    rm -Rf "${WORK_DIR}/${SAMPLE_NAME}-reads.tmp"
    mkdir "${WORK_DIR}/${SAMPLE_NAME}-reads.tmp"
    
    # Simulate reads
    time "${PBSIM}" \
        ${PBSIM_PARAMS} \
       --sample-fastq "${SAMPLE_FASTQ}" \
       --prefix "${WORK_DIR}/${SAMPLE_NAME}-reads.tmp/sim" \
       "${WORK_DIR}/${SAMPLE_NAME}.fa"
    
    mv "${WORK_DIR}/${SAMPLE_NAME}-reads.tmp" "${WORK_DIR}/${SAMPLE_NAME}-reads"
fi

function do_job() {
    # Run this file in a job
    set -e
    
    SAM_NAME="${MAF_NAME%.maf}.sam"
    FASTQ_NAME="${MAF_NAME%.maf}.fastq"
    REF_NAME="${MAF_NAME%.maf}.ref"
    RENAMED_BAM_NAME="${MAF_NAME%.maf}.renamed.bam"
    # Get the contig name in the format it would be as a reference sense path.
    # It may already be a reference sense path.
    # Can't run under pipefail because some of these may not match.
    CONTIG_NAME="$(cat "${REF_NAME}" | head -n1 | sed 's/^>//' | sed 's/ .*//' | sed 's/#\([0-9]*\)$/[\1]/')"
    # Haplotype paths can end in a 0 offset/fragment but reference paths don't include that in the name.
    CONTIG_NAME="${CONTIG_NAME%\[0\]}"
    if [[ ! -e "${RENAMED_BAM_NAME}" ]] ; then
        echo "Making ${RENAMED_BAM_NAME}..."
        if [[ ! -e "${SAM_NAME}" ]] ; then
            echo "Making SAM ${SAM_NAME}..."
            /usr/bin/time -v bioconvert maf2sam "${MAF_NAME}" "${SAM_NAME}.tmp" 
            mv "${SAM_NAME}.tmp" "${SAM_NAME}"
        fi
        set -o pipefail
        ${ADD_QUALITIES} -s "${SAM_NAME}" -f "${FASTQ_NAME}" | sed "s/ref/${CONTIG_NAME}/g" | samtools view -b - > "${RENAMED_BAM_NAME}.tmp"
        set +o pipefail
        mv "${RENAMED_BAM_NAME}.tmp" "${RENAMED_BAM_NAME}"
    else
        echo "Already have ${RENAMED_BAM_NAME}..."
    fi
}

# Convert all the reads to BAM in the space of the sample as a primary reference
for MAF_NAME in "${WORK_DIR}"/${SAMPLE_NAME}-reads/sim_*.maf ; do
    if [[ "${MAX_JOBS}" == "1" ]] ; then
        # Serial mode
        do_job
    else
        # Parallel mode
        while [[ "$(jobs -p | wc -l)" == "${MAX_JOBS}" ]] ; do
            # Don't do too much in parallel
            # Fake wait on any job without wait -n
            sleep 0.5
        done
        (
            do_job
        ) &
        ((RUNNING_JOBS += 1))
    fi
done
# Wait on all jobs
wait

if [[ ! -e "${WORK_DIR}/${SAMPLE_NAME}-reads/injected.gam" ]] ; then
    # Move reads into graph space
    time vg inject -x "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz" "${WORK_DIR}/${SAMPLE_NAME}-reads/merged.bam" -t 16 >"${WORK_DIR}/${SAMPLE_NAME}-reads/injected.gam.tmp"
    mv "${WORK_DIR}/${SAMPLE_NAME}-reads/injected.gam.tmp" "${WORK_DIR}/${SAMPLE_NAME}-reads/injected.gam"
fi

if [[ ! -e "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam" ]] ; then
    # Annotate reads with linear reference positions
    time vg annotate -x "${WORK_DIR}/${GRAPH_NAME}.gbz" -a "${WORK_DIR}/${SAMPLE_NAME}-reads/injected.gam" --multi-position -l 100 -t 16 >"${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam.tmp"
    mv "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam.tmp" "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam"
fi

# Subset to manageable sizes (always)
vg filter -d 1.0004 "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam" | vg view -aj - | shuf | head -n100 | vg view -JGa - > "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}-100.gam"
vg filter -d 2.004 "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam" | vg view -aj - | shuf | head -n1000 | vg view -JGa - > "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}-1000.gam"
vg filter -d 3.04 "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam" | vg view -aj - | shuf | head -n10000 | vg view -JGa - > "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}-10000.gam"

# Output them
mkdir -p "${OUT_DIR}"
cp "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}.gam" "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}-100.gam" "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}-1000.gam" "${WORK_DIR}/${SAMPLE_NAME}-reads/${SAMPLE_NAME}-sim-${TECH_NAME}-10000.gam" "${OUT_DIR}/"

if [[ "${CLEAN_WORK_DIR}" == "1" ]] ; then
    # Clean up the work directory
    rm -Rf "${WORK_DIR}"
fi

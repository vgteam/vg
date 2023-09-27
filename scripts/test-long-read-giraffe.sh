#!/usr/bin/env bash

# Script to run Giraffe in long read mose on a set of simulated reads and evaluate its speed and accuracy.

set -ex

: "${DATA_DIR:="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe"}"
: "${GRAPH_BASE:="${DATA_DIR}/graphs/hprc-v1.1-mc-chm13.d9"}"
: "${MINPARAMS:="k31.w50.W"}"
: "${CONDITION:="zip-bugfix"}"
# Our GAM file for writing our mapped reads to
: "${GAM_FILE:="trash/mapped-${CONDITION}.gam"}"
: "${INPUT_READS:="${DATA_DIR}/reads/sim/hifi/HG002/HG002-sim-hifi-1000.gam"}"
: "${GIRAFFE_ARGS:=""}"

# Make absolute paths before changing directories
DATA_DIR="$(abspath "${DATA_DIR}")"
GRAPH_BASE="$(abspath "${GRAPH_BASE}")"
GAM_FILE="$(abspath "${GAM_FILE}")"
INPUT_READS="$(abspath "${INPUT_READS}")"

if which sbatch >/dev/null 2>&1 ; then
    # Slurm is available.
    # Put your Slurm command arguments in a JOB_ARGS array and run do_sbatch or
    # do_srun with your command.

    # Run a command wrapped with sbatch
    function do_sbatch() {
        sbatch "${JOB_ARGS[@]}" --wrap "${1}"
    }

    # Run a command and wait on it with srun
    function do_srun() {
        srun "${JOB_ARGS[@]}" "$@"
    }

    # Wait for Slurm jobs to be done and their changes to be visible on disk
    function swait() {
        QUEUE_LINES=0
        while [[ "${QUEUE_LINES}" != "1" ]] ; do
            # On the first loop, or on subsequent loops when running or pending jobs are visible

            # Wait
            sleep 2
            # Check again
            QUEUE_LINES="$(squeue -u $USER | wc -l)"
        done
        # Hope filesystem is no more than this many seconds behind Slurm
        sleep 10
    }

else
    # No Slurm. Run everything locally.

    # Run a quoted command in the backgorund
    function do_sbatch() {
        bash -c "${1}" &
    }

    # Run a command in the foreground
    function do_srun() {
        "$@"
    }

    # Wait on all jobs
    function swait() {
        wait
    }

fi



# Go to the main vg directory
cd "$(dirname -- "$0")"
cd ..

rm -f *.out
JOB_ARGS=(-c16 --mem 400G --job-name zipcode-run)
do_sbatch "time vg giraffe --parameter-preset lr --progress --track-provenance -Z ${GRAPH_BASE}.gbz -d ${GRAPH_BASE}.dist -m ${GRAPH_BASE}.${MINPARAMS}.withzip.min -z ${GRAPH_BASE}.${MINPARAMS}.zipcodes -G ${INPUT_READS} -t16 ${GIRAFFE_ARGS} >${GAM_FILE}"

swait

EXP_DIR="trash/${CONDITION}"
OUT_DIR="${EXP_DIR}/hifi-${CONDITION}"
rm -Rf "${OUT_DIR}"
rm -Rf "${EXP_DIR}"
mkdir -p "${OUT_DIR}"

JOB_ARGS=(-c 3 --mem 10G)
for STAGE in minimizer seed tree fragment chain align winner ; do
    [[ -e "${OUT_DIR}/read-time-${STAGE}.tsv" ]] || do_sbatch "set -e; vg view -aj "${GAM_FILE}" | jq -r '.annotation.stage_'${STAGE}'_time' >${OUT_DIR}/read-time-${STAGE}.tsv"
done
[[ -e "${OUT_DIR}/read-time-to-chain.tsv" ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.annotation.stage_minimizer_time + .annotation.stage_seed_time + .annotation.stage_bucket_time + .annotation.stage_fragment_time + .annotation.stage_chain_time' >${OUT_DIR}/read-time-to-chain.tsv"

    

[[ -e "${OUT_DIR}"/read-best-chain-coverage.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.annotation.best_chain_coverage' > ${OUT_DIR}/read-best-chain-coverage.tsv"
[[ -e "${OUT_DIR}"/read-best-chain-longest-jump.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.annotation.best_chain_longest_jump' > ${OUT_DIR}/read-best-chain-longest-jump.tsv"
[[ -e "${OUT_DIR}"/read-best-chain-average-jump.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.annotation.best_chain_average_jump' > ${OUT_DIR}/read-best-chain-average-jump.tsv"
[[ -e "${OUT_DIR}"/read-best-chain-anchors.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.annotation.best_chain_anchors' > ${OUT_DIR}/read-best-chain-anchors.tsv"
[[ -e "${OUT_DIR}"/read-best-chain-anchor-length.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.annotation.best_chain_anchor_length' > ${OUT_DIR}/read-best-chain-anchor-length.tsv"
[[ -e "${OUT_DIR}"/read-score.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '.score // 0' > ${OUT_DIR}/read-score.tsv"
[[ -e "${OUT_DIR}"/read-unclipped.tsv ]] || do_sbatch "set -e; vg view -aj ${GAM_FILE} | jq -r '1.0 - (([[.path.mapping[0].edit[0], .path.mapping[-1].edit[-1]][] | select(.from_length // 0 == 0) | select(.sequence) | .to_length] + [0] | add) / (.sequence | length))' > ${OUT_DIR}/read-unclipped.tsv"

swait

PLOT_DIR="${EXP_DIR}/plots"
mkdir -p "${PLOT_DIR}"

do_sbatch "set -e; histogram.py ${OUT_DIR}/read-best-chain-coverage.tsv --bins 100 --title '${CONDITION} Fraction Covered' --y_label 'Items' --x_label 'Coverage' --no_n --save ${PLOT_DIR}/read-best-chain-coverage-${CONDITION}.png"
do_sbatch "set -e; histogram.py ${OUT_DIR}/read-best-chain-longest-jump.tsv --bins 100 --title '${CONDITION} Longest Jump' --y_label 'Items' --x_label 'Jump (bp)' --no_n --save ${PLOT_DIR}/read-best-chain-longest-jump-${CONDITION}.png"
do_sbatch "set -e; histogram.py ${OUT_DIR}/read-best-chain-average-jump.tsv --bins 100 --title '${CONDITION} Average Jump' --y_label 'Items' --x_label 'Jump (bp)' --no_n --save ${PLOT_DIR}/read-best-chain-average-jump-${CONDITION}.png"
do_sbatch "set -e; histogram.py ${OUT_DIR}/read-best-chain-anchors.tsv --bins 100 --title '${CONDITION} Chained Anchors' --y_max 60 --y_label 'Items' --x_label 'Anchors (count)' --no_n --save ${PLOT_DIR}/read-best-chain-anchors-${CONDITION}.png"
do_sbatch "set -e; histogram.py ${OUT_DIR}/read-best-chain-anchor-length.tsv --bins 100 --title '${CONDITION} Chained Anchor Length' --y_max 60 --y_label 'Items' --x_label 'Anchor Length (bp)' --no_n --save ${PLOT_DIR}/read-best-chain-anchor-length-${CONDITION}.png"
do_sbatch "set -e; histogram.py ${OUT_DIR}/read-score.tsv --bins 100 --title '${CONDITION} Score' --y_label 'Items' --x_label 'Score' --no_n --save ${PLOT_DIR}/read-score-${CONDITION}.png"
do_sbatch "set -e; histogram.py ${OUT_DIR}/read-unclipped.tsv --bins 100 --title '${CONDITION} Portion Unclipped' --y_label 'Items' --x_label 'Portion Unclipped' --no_n --save ${PLOT_DIR}/read-unclipped-${CONDITION}.png"

do_sbatch "set -e; histogram.py ${OUT_DIR}/read-time-to-chain.tsv --bins 100 --title '${CONDITION} Time To Chain' --x_max 5 --y_label 'Items' --x_label 'Time (s)' --no_n --save ${PLOT_DIR}/read-time-to-chain-${CONDITION}.png"

swait

printf "#Condition\tminimizer_time\tseed_time\ttree_time\tfragment_time\tchain_time\talign_time\twinner_time\n" > "${PLOT_DIR}/stats.tsv"

printf "${CONDITION}\t${REPLICATE}\t" >>"${PLOT_DIR}/stats.tsv"

for STAGE in minimizer seed tree fragment chain align winner ; do 
    echo ${OUT_DIR}/read-time-${STAGE}.tsv
    printf "$(cat "${OUT_DIR}/read-time-${STAGE}.tsv" | mean.sh)\t" >>"${PLOT_DIR}/stats.tsv"
done
 printf "\n" >>"${PLOT_DIR}/stats.tsv"

cat "${PLOT_DIR}/stats.tsv"

JOB_ARGS=(-c16 --mem 20G)
do_srun vg annotate -a ${GAM_FILE} -x ${GRAPH_BASE}.gbz -m >${GAM_FILE%.gam}.annotated.gam
do_srun vg gamcompare --range 200 ${GAM_FILE%.gam}.annotated.gam ${INPUT_READS} -T -a "${CONDITION}" -o ${GAM_FILE%.gam}.compared.gam > ${GAM_FILE%.gam}.compared.tsv

# Now make a PR plot stratified by MAPQ
Rscript scripts/plot-pr.R ${GAM_FILE%.gam}.compared.tsv ${GAM_FILE%.gam}.compared.svg






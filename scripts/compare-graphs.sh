#!/usr/bin/env bash
# compare-graphs.sh: compare a set of graph against each other using toil-vg mapeval on AWS

set -e

# Take the base names for the graphs
GRAPH_NAMES=( snp1kg primary )

# And the region to run
REGION_NAME="BRCA1"

# What cluster should we use?
CLUSTER_NAME="novakcluster2"
# What keypair should we set it up with?
KEYPAIR_NAME="anovak"

# Where do we keep our input files
INPUT_STORE="https://cgl-pipeline-inputs.s3.amazonaws.com/vg_cgl/bakeoff"

# Where do we save our results from the various jobs responsible for writing them?
OUTPUT_STORE="aws:us-west-2:cgl-pipeline-inputs/vg_cgl/comparison-script/run"
OUTPUT_STORE_URL="s3://cgl-pipeline-inputs/vg_cgl/comparison-script/run"

# Where do we store our jobs?
JOB_TREE="aws:us-west-2:amntree1"

function get_input_url() {
    # Prints the input URL to download for the given file name
    local BASE_FILENAME="${1}"
    shift
    echo "${INPUT_STORE}/${BASE_FILENAME}"
}

function get_graph_url() {
    # Prints the base URL for the given graph
    local BASE_GRAPHNAME="${1}"
    shift
    get_input_url "${BASE_GRAPHNAME}-${REGION_NAME}"
}

# Convert just graph stems to full base urls
GRAPH_URLS=()
for GRAPH_STEM in "${GRAPH_NAMES[@]}"; do
    GRAPH_URLS+=(`get_graph_url "${GRAPH_STEM}"`)
done

TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:latest toil launch-cluster "${CLUSTER_NAME}" --nodeType=t2.micro -z us-west-2a "--keyPairName=${KEYPAIR_NAME}"

# We need to manually install git to make pip + git work...
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt update
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" apt install git -y

# For hot deployment to work, toil-vg needs to be in a virtualenv that can see the system Toil
toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" virtualenv --system-site-packages venv

toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/pip install git+https://github.com/adamnovak/toil-vg.git@2fae08f2cdb86e689244e33cf8a6705b01f42a6c

# We need the master's IP to make Mesos go
MASTER_IP="$(toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" hostname -i)"

toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg mapeval \
    --fasta `get_input_url "${REGION_NAME}.fa"` \
    --index-bases "${GRAPH_URLS[@]}" \
    --gam-names "${GRAPH_NAMES[@]}" \
    --gam_input_reads `get_input_url "comparison-${REGION_NAME}.gam"` \
    --bwa --bwa-paired --vg-paired \
    --mapeval-threshold 200 \
    --realTimeLogging --logInfo \
    "${JOB_TREE}" \
    "${OUTPUT_STORE}" \
    `get_input_url "comparison-${REGION_NAME}.pos"` \
    --batchSystem mesos --provisioner=aws "--mesosMaster=${MASTER_IP}:5050" --nodeType=t2.large
    
mkdir -p ./out
aws s3 sync "${OUTPUT_STORE_URL}" ./out
aws s3 rm --recursive "${OUTPUT_STORE_URL}"

toil clean "${JOB_TREE}"

toil destroy-cluster "${CLUSTER_NAME}" -z us-west-2a

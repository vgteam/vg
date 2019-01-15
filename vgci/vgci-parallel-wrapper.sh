#!/usr/bin/env bash

# vgci-parallel-wrapper.sh: run a subset of the available tests with vgci.sh
# Usage: vgci-paralell-wrapper.sh vgci/test-list.txt vgci-docker-vg-local.tar.gz $CI_NODE_INDEX $CI_NODE_TOTAL junit
# Finds our share of the tests in the test list, given our 0-based index and the total node count.
# Runs each of them with the given loadable Docker image and drops a unique junit xml for the test in the given output directory.
# If we have to run more than 1 test we will duplicate vgci.sh's setup work creating venvs and so on.
# Meant to be run fromn the vg project root.

# If any test fails, stop
set +e
# Report what we're up to
set -x

# Parse arguments
TEST_LIST="${1}"
shift
DOCKER_ARCHIVE="${1}"
shift
NODE_INDEX="${1}"
shift
NODE_TOTAL="${1}"
shift
OUT_DIR="${1}"
shift

# Count up all the tests we run
TEST_NUMBER=0

# Get every NODE_TOTAL'th line starting at NODE_INDEX
# See https://superuser.com/a/396557
sed -n "${NODE_INDEX}~${NODE_TOTAL}p" "${TEST_LIST}" | while read TEST_SPEC
do
    # And for each
   
    # Run the main vgci script for that test
    # TODO: We're hardcoding the Docker file here...
    vgci/vgci.sh -D "${DOCKER_ARCHIVE}" -t "${TEST_SPEC}" -j "${OUT_DIR}/junit.${NODE_INDEX}.${TEST_NUMBER}.xml" -H

    # Number the next test differently
    ((TEST_NUMBER+=1))
done



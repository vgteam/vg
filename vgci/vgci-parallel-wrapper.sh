#!/usr/bin/env bash

# vgci-parallel-wrapper.sh: run a subset of the available tests with vgci.sh
# Usage: vgci-paralell-wrapper.sh vgci/test-list.txt vgci-docker-vg-local $CI_NODE_INDEX $CI_NODE_TOTAL junit test_output
# Finds our share of the tests in the test list, given our 0-based index and the total node count.
# Runs each of them with the given Docker tag and drops a unique junit xml for the test in the given output 
# directory, while keeping the test output in the other given output directory in per-test folders.
# If we have to run more than 1 test we will duplicate vgci.sh's setup work creating venvs and so on.
# Meant to be run fromn the vg project root.

# We will handle errors ourselves
set -e
# Report what we're up to
set -x

# Parse arguments
TEST_LIST="${1}"
shift
DOCKER_TAG="${1}"
shift
NODE_INDEX="${1}"
shift
NODE_TOTAL="${1}"
shift
JUNIT_OUT_DIR="${1}"
shift
TEST_OUT_DIR="${1}"
shift

# Count up all the tests we run
TEST_NUMBER=0

# Get every NODE_TOTAL'th line starting at NODE_INDEX
# See https://superuser.com/a/396557
for TEST_SPEC in "$(sed -n "${NODE_INDEX}~${NODE_TOTAL}p" "${TEST_LIST}")"
do
    # And for each
   
    # Run the main vgci script for that test
    vgci/vgci.sh -T "${DOCKER_TAG}" -t "${TEST_SPEC}" -j "${JUNIT_OUT_DIR}/junit.${NODE_INDEX}.${TEST_NUMBER}.xml" -w "${TEST_OUT_DIR}" -H
    TEST_EXIT="${?}"
    
    if [ "${TEST_EXIT}" != "0" ]
    then
        echo "Test ${TEST_SPEC} failed!"
        exit "${TEST_EXIT}"
    fi

    # Number the next test differently
    ((TEST_NUMBER+=1))
done



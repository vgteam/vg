#!/usr/bin/env bash

# Script to build a test Docker and run tests on it.
# Tests are run using toil-vg: https://github.com/BD2KGenomics/toil-vg
# Drops test-report.xml with the test results in the current directory.
# This script is hooked into http://vgci.cgcloud.info
# Note: we assume we run this in vg/ ie inside the vg directory we want to test

# Keep going on errors; we make sure to return the right status code.
set +e

# Should we build and run locally, or should we use Docker?
LOCAL_BUILD=0
# What filename should we export our Docker to?
SAVE_DOCKER=""
# What filename should we import our Docker from instead of building.
# Note that this must be exported from this script, so it will have the right tag!
LOAD_DOCKER=""
# What tag will our Docker use?
# We need exclusive control of the Docker daemon for the duration of the test, so nobody moves it!
DOCKER_TAG=""
# If not specified, use this default
DOCKER_TAG_DEFAULT="vgci-docker-vg-local"
# Should we re-use and keep around the same virtualenv?
REUSE_VENV=0
# Should we keep our test output around after uploading the new baseline?
KEEP_OUTPUT=0
# Should we keep all intermediate output (ie --force_outstore in toil-vg)?
KEEP_INTERMEDIATE_FILES=0
# Should we show stdout and stderr from tests? If so, set to "-s".
SHOW_OPT=""
# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/vgteam/toil-vg.git@d16da00b92c491f90433e151cb4f5a89a44395b8"
# What toil should we install?
# Could be something like "toil[aws,mesos]==3.20.0"
# or "git+https://github.com/DataBiosphere/toil.git@3ab74776a3adebd6db75de16985ce9d734f60743#egg=toil[aws,mesos]"
TOIL_PACKAGE="git+https://github.com/DataBiosphere/toil.git@9041c54d9802add8973e31eddd8ba33529218f94#egg=toil[aws,mesos]"
# What tests should we run?
# Should be something like "vgci/vgci.py::VGCITest::test_sim_brca2_snp1kg_mpmap"
# Must have the Python file in it or Pytest can't find the tests.
# May contain multiple space-separated test specifiers.
PYTEST_TEST_SPEC="vgci/vgci.py"
# What scratch directory should we use to run the tests?
# If unset we use vgci_work and don't persist it.
# If set, must start with / or ./ because it is used to make toil-vg IOstores
SAVE_WORK_DIR=
# What test result directory shoudl we load, if any?
LOAD_WORK_DIR=
# Save JUnit test report to this file
SAVE_JUNIT=""
# Import JUnit test report from this file instead of running tests.
LOAD_JUNIT=""
# Should we analyze the junit test report and post our own HTML report?
CREATE_REPORT=1
# What S3 URL does test output go to?
OUTPUT_DESTINATION="s3://vg-data/vg_ci"
# What bucket owner account ID if any should be granted full control of uploaded objects?
OUTPUT_OWNER="b1cf5e10ba0aeeb00e5ec70b3532826f22a979ae96c886d3081d0bdc1f51f67e"
# Should the output be marked public? If we're uploading to an S3 bucket where
# this isn't allowed, we'll get a permission denied error when we try to do the
# upload.
OUTPUT_PUBLIC=0

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] \n"
    printf "Options:\n\n"
    printf "\t-l\t\tBuild vg locally (instead of in Docker) and don't use Docker at all.\n"
    printf "\t\t\tNon-Python dependencies must be installed.\n"
    printf "\t-d FILE\tSave built Docker to the given file.\n"
    printf "\t-D FILE\tLoad a previously built Docker from a file instead of building.\n"
    printf "\t-T TAG\tLoad a previously built Docker from the given tag/specifier instead of building.\n"
    printf "\t-r\t\tRe-use virtualenvs across script invocations. \n"
    printf "\t-k\t\tKeep on-disk output from tests. \n"
    printf "\t-i\t\tKeep intermediate on-disk output from tests. \n"
    printf "\t-s\t\tShow test output and error streams (pass -s to pytest). \n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg.\n"
    printf "\t-t TESTSPEC\tUse the given PyTest test specifier(s), space-separated, to select tests to run, or 'None' for no tests.\n"
    printf "\t-w WORKDIR\tOutput test result data to the given absolute or ./ path (also used for scratch)\n"
    printf "\t-W WORKDIR\tLoad test result data from the given path instead of building or running tests\n"
    printf "\t-j FILE\tSave the JUnit test report XML to the given file (default: test-report.xml)\n"
    printf "\t-J FILE\tLoad the JUnit test report from the given file instead of building or running tests\n"
    printf "\t-H\tSkip generating HTML report based on JUnit report\n"
    exit 1
}

while getopts "ld:D:T:rkisp:t:w:W:j:J:H" o; do
    case "${o}" in
        l)
            LOCAL_BUILD=1
            ;;
        d)
            SAVE_DOCKER="${OPTARG}"
            ;;
        D)
            LOAD_DOCKER="${OPTARG}"
            ;;
        T)
            DOCKER_TAG="${OPTARG}"
            ;;
        r)
            REUSE_VENV=1
            ;;
        k)
            KEEP_OUTPUT=1
            ;;
        i)
            KEEP_INTERMEDIATE_FILES=1
            ;;	  
        s) 
            SHOW_OPT="-s"
            ;;
        p)
            TOIL_VG_PACKAGE="${OPTARG}"
            ;;
        t)
            PYTEST_TEST_SPEC="${OPTARG}"
            ;;
        w)
            SAVE_WORK_DIR="${OPTARG}"
            ;;
        W)
            LOAD_WORK_DIR="${OPTARG}"
            ;;
        j)
            SAVE_JUNIT="${OPTARG}"
            ;;
        J)
            LOAD_JUNIT="${OPTARG}"
            ;;
        H)
            CREATE_REPORT=0
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

#########
# COMMON PREP PHASE
#########

if [ ! -e ~/.aws/credentials ] && [ -z "${CI}" ]
then
    # We're not on cloud CI and we have no AWS credentials.
    >&2 echo "WARNING: No AWS credentials at ~/.aws/credentials; test data may not be able to be downloaded!"
fi

if [ -e /sys/fs/cgroup/cpu/cpu.cfs_quota_us ] && [ -e /sys/fs/cgroup/cpu/cpu.cfs_period_us ]
then
    # If confined to a container, use the container's CPU limit (if >= 1)
    NUM_CORES=$(($(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us) / $(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)))
fi

if [ "${NUM_CORES}" == "0" ]
then
    PLATFORM=`uname -s`
    if [ $PLATFORM == "Darwin" ]
    then
        NUM_CORES=`sysctl -n hw.ncpu`
    else
        NUM_CORES=`cat /proc/cpuinfo | grep "^processor" | wc -l`
    fi
fi

if [ "${NUM_CORES}" == "0" ]
then
    echo "could not determine NUM_CORES, using 2"
	NUM_CORES=2
fi

# Add these arguments to grant ownership of uploads
GRANT_ARGS=()
if [[ ! -z "${OUTPUT_OWNER}" || "${OUTPUT_PUBLIC}" == "1" ]] ; then
    GRANT_ARGS=(--grants)
    if [[ "${OUTPUT_PUBLIC}" == "1" ]] ; then
        GRANT_ARGS+=("read=uri=http://acs.amazonaws.com/groups/global/AllUsers")
    fi
    if [[ ! -z "${OUTPUT_OWNER}" ]] ; then
        GRANT_ARGS+=("full=id=${OUTPUT_OWNER}")
    fi
fi

# We have 3 phases: build, test, and report.
# Each has its own prep work and maybe some associated cleanup.
# We have to figure out which of them need doing.
DO_BUILD=1
DO_TEST=1
DO_REPORT=1

if [ ! -z "${LOAD_DOCKER}" ]
then
    # Skip the build phase
    DO_BUILD=0
fi

if [ ! -z "${DOCKER_TAG}" ]
then
    # Skip build and use this tag
    DO_BUILD=0
else
    # Use the default tag
    DOCKER_TAG="${DOCKER_TAG_DEFAULT}"
fi

if [ "${PYTEST_TEST_SPEC}" == "None" ]
then
    # Skip the test phase because we were asked to
    DO_TEST=0
fi

if [ ! -z "${LOAD_JUNIT}" ] || [ ! -z "${LOAD_WORK_DIR}" ]
then
    # Skip the build phase because we don't even need to run anything
    DO_BUILD=0
    # Skip the test phase because we'll just load stuff
    DO_TEST=0
fi

if [ "${CREATE_REPORT}" == "0" ]
then
    # Skip the report phase
    DO_REPORT=0
fi


# Also each phase can fail with a nonzero status
BUILD_FAIL=0
TEST_FAIL=0
REPORT_FAIL=0

if [ "${DO_BUILD}" != "0" ]
then

    echo "VGCI: Run build"

    #########
    # BUILD PREP PHASE
    #########

    # Make sure we have submodules
    git submodule update --init --recursive

    #########
    # BUILD PHASE
    #########

    if [ "${LOCAL_BUILD}" == "1" ]
    then
        # Just build vg here
        make -j ${NUM_CORES}

        if [ "$?" -ne 0 ]
        then
            echo "vg local build fail"
            BUILD_FAIL=1
        fi
    else

        # Build a docker image locally.  Can be useful when don't
        # have priveleges to easily install dependencies
        
        # Build the git version file first, so the Docker knows its version
        make version

        docker pull mirror.gcr.io/library/ubuntu:20.04
        docker build --no-cache -t "${DOCKER_TAG}" -f Dockerfile .
        if [ "$?" -ne 0 ]
        then
            echo "vg docker build fail"
            BUILD_FAIL=1
        else
        
            if [ ! -z "${SAVE_DOCKER}" ]
            then
                # Save the Docker to a file.
                # Loading it will set the tag.
                docker save "${DOCKER_TAG}" -o "${SAVE_DOCKER}"
                if [ "$?" -ne 0 ]
                then
                    echo "vg docker save fail"
                    BUILD_FAIL=1
                fi
            fi
        fi
    fi
fi

if [ "${BUILD_FAIL}" != "0" ]
then
    # The build failed, so skip the test
    DO_TEST=0
fi

if ([ "${DO_TEST}" != "0" ] || [ "${DO_REPORT}" != "0" ]) && [ ! -z "${LOAD_DOCKER}" ]
then
    # Just load the Docker instead of building.
    # It will set the tag it was saved from.
    # We need it both for testing and reportign because both need the vg version
    docker load -i "${LOAD_DOCKER}"
fi

# Create Toil venv
if [ ! "${REUSE_VENV}" == "1" ]; then
    rm -rf .env
fi
if [ ! -e .env ]; then
    virtualenv --python=python3  .env
fi
. .env/bin/activate

# Upgrade pip so that it can use the wheels for numpy & scipy, so that they
# don't try to build from source
pip3 install --upgrade pip setuptools==45.0.0
pip3 install awscli

if [ "${DO_TEST}" != "0" ]
then

    echo "VGCI: Run test"

    #########
    # TEST PREP PHASE
    #########

    # Prepare directory for temp files (assuming cgcloud file structure)
    # Sometimes the instances have un-deletable files in tmp, so we continue through errors
    if [ -d "/mnt/ephemeral" ]
    then
         TMPDIR=/mnt/ephemeral/tmp
         rm -rf $TMPDIR
         mkdir -p $TMPDIR
         export TMPDIR
    fi

    # Dependencies for running tests.  Need numpy, scipy and sklearn for
    # running toil-vg mapeval, and dateutils and reqests for
    # ./mins_since_last_build.py. Need exactly the right version of requests
    # for the Python Docker API to work (see
    # <https://github.com/docker/docker-py/issues/3256>).
    pip3 install pytest timeout_decorator 'requests==2.31.0' dateutils
    # TODO: To upgrade to Numpy 1.24+, we need to remove usages of `numpy.int`
    # AKA `np.int` from toil-vg. See <https://stackoverflow.com/a/74946903>
    # Make sure to get the giant buildable modules as binaries only to avoid
    # wasting CI time building some slightly nicer version Pip prefers.
    pip3 install 'numpy==1.23.5' scipy scikit-learn --only-binary :all:

    # Install Toil
    echo "Installing toil from ${TOIL_PACKAGE}"
    pip3 install --upgrade "${TOIL_PACKAGE}"
    if [ "$?" -ne 0 ]
    then
        echo "pip3 install toil fail"
        exit 1
    fi

    # Don't manually install boto since toil just installs its preferred version

    # Install toil-vg itself
    echo "Installing toil-vg from ${TOIL_VG_PACKAGE}"
    pip3 install --upgrade "${TOIL_VG_PACKAGE}"
    if [ "$?" -ne 0 ]
    then
        echo "pip3 install toil-vg fail"
        exit 1
    fi
    
    #########
    # TEST PHASE
    #########

    # we pass some parameters through pytest by way of our config file
    # in particular, we set the vg version and cores, and specify
    printf "cores ${NUM_CORES}\n" > vgci_cfg.tsv
    printf "teardown False\n" >> vgci_cfg.tsv
    printf "workdir ${SAVE_WORK_DIR:-./vgci-work}\n" >> vgci_cfg.tsv
    if [ "${KEEP_INTERMEDIATE_FILES}" == "0" ]; then
        printf "force_outstore False\n" >> vgci_cfg.tsv
    else
        printf "force_outstore True\n" >> vgci_cfg.tsv
    fi
    #printf "verify False\n" >> vgci_cfg.tsv
    #printf "baseline ./vgci-baseline\n" >> vgci_cfg.tsv

    if [ "${LOCAL_BUILD}" == "1" ]
    then
        # Test the locally built vg
        VG_VERSION=`vg version -s`
        printf "vg-docker-version None\n" >> vgci_cfg.tsv
        printf "container None\n" >> vgci_cfg.tsv
    else
        # Test the Dockerized vg
        VG_VERSION=`docker run ${DOCKER_TAG} vg version -s`
        printf "vg-docker-version ${DOCKER_TAG}\n" >> vgci_cfg.tsv
        
        # Pull down the docker images, so time costs (and instability) of doing so doesn't affect
        # individual test results (looking at you, rocker/tidyverse:3.4.2)
        # Allow two tries
        for img in $(toil-vg generate-config | grep docker: | grep -v vg | awk '{print $2}' | sed "s/^\([\"']\)\(.*\)\1\$/\2/g"); do docker pull $img ; done
        for img in $(toil-vg generate-config | grep docker: | grep -v vg | awk '{print $2}' | sed "s/^\([\"']\)\(.*\)\1\$/\2/g"); do docker pull $img ; done
    fi

    mkdir -p "${SAVE_WORK_DIR:-./vgci-work}"
    
    # run the tests, output the junit report 
    rm -f test-report.xml
    pytest -vv ${PYTEST_TEST_SPEC} --junitxml=test-report.xml ${SHOW_OPT}
    TEST_FAIL="$?"
    
    if [ ! -z "${SAVE_JUNIT}" ]
    then
        # Copy it to the destination
        cp test-report.xml "${SAVE_JUNIT}" || touch "${SAVE_JUNIT}"
    fi
    
    if [ ! -z "${CI}" ] && [ "${TEST_FAIL}" != 0 ]
    then
        # We are running on cloud CI (and not manually running the tests), so
        # we probably have AWS and Github credentials and can upload stuff to S3.
        # A test faled, so we should make sure we upload its outstore for debugging.
        # TODO: If we get the report job to always run and include this, maybe we don't need individual uploads too.
        
        # Upload the results of this test in particular, as soon as it is done, instead of waiting for the final report job to do it.
        tar czf "test_output.tar.gz" "${SAVE_WORK_DIR}/" test-report.xml
        DEST_URL="${OUTPUT_DESTINATION}/vgci_output_archives/${VG_VERSION}/${CI_PIPELINE_ID}/${CI_JOB_ID}/test_output.tar.gz"
	aws s3 cp --only-show-errors \
            "test_output.tar.gz" "${DEST_URL}" "${GRANT_ARGS[@]}"
            
        
        echo "Test(s) failed. Output is available at ${DEST_URL}"
    fi

    # Load from the work directory we saved to
    LOAD_WORK_DIR="${SAVE_WORK_DIR:-./vgci-work}"
fi

if [ ! -z "${LOAD_JUNIT}" ]
then
    # Load up the input JUnit report
    cp "${LOAD_JUNIT}" test-report.xml
fi

if [ "${DO_REPORT}" != "0" ]
then

    echo "VGCI: Run report"

    #########
    # REPORT PREP PHASE
    #########

    pip3 install pygithub
    
    #########
    # REPORT PHASE
    #########

    # Generate a report in two files: HTML full output, and a Markdown summary.
    # Takes as input the test result XML and the work directory with the
    # test output files.
    python3 vgci/mine-logs.py test-report.xml "${LOAD_WORK_DIR}/" report-html/ summary.md
    if [ "$?" -ne 0 ]
    then
        echo "Log mining fail"
        REPORT_FAIL=1
    fi

    if [ ! -z "${CI}" ]
    then
        # We are running on cloud CI (and not manually running the tests), so
        # we probably have AWS and Github credentials and can upload stuff to S3.
        
        # Put the report on Github for the current pull request or commit.
        vgci/post-report report-html summary.md
        if [ "$?" -ne 0 ]
        then
            echo "Report posting fail"
            REPORT_FAIL=1
        fi
        
        # We need to re-determine the vg version since we may not have run the tests.
        if [ "${LOCAL_BUILD}" == "1" ]
        then
            # Report on the locally built vg
            VG_VERSION=`vg version -s`
        else
            # Test the Dockerized vg
            VG_VERSION=`docker run ${DOCKER_TAG} vg version -s`
        fi
        
        # we publish the results to the archive
        tar czf "${VG_VERSION}_output.tar.gz" "${LOAD_WORK_DIR}/" test-report.xml vgci/vgci.py vgci/vgci.sh vgci_cfg.tsv
        aws s3 cp --only-show-errors \
            "${VG_VERSION}_output.tar.gz" "${OUTPUT_DESTINATION}/vgci_output_archives/" "${GRANT_ARGS[@]}"
        if [ "$?" -ne 0 ]
        then
            echo "Archive upload fail"
            REPORT_FAIL=1
        fi

        # if we're merging the PR (and not just testing it), we publish results to the baseline
        if [ -z "${CI_MERGE_REQUEST_IID}" ] && [ "${CI_COMMIT_REF_NAME}" == "master" ]
        then
            echo "Updating baseline"
            aws s3 sync --only-show-errors --delete \
                "${LOAD_WORK_DIR}/" "${OUTPUT_DESTINATION}/vgci_regression_baseline" "${GRANT_ARGS[@]}"
            if [ "$?" -ne 0 ]
            then
                echo "Baseline upload fail"
                REPORT_FAIL=1
            fi
        
            printf "${VG_VERSION}\n" > "vg_version_${VG_VERSION}.txt"
            printf "${CI_COMMIT_TITLE}" >> "vg_version_${VG_VERSION}.txt"
            aws s3 cp --only-show-errors \
                "vg_version_${VG_VERSION}.txt" "${OUTPUT_DESTINATION}/vgci_regression_baseline/" "${GRANT_ARGS[@]}"
            if [ "$?" -ne 0 ]
            then
                echo "Version upload fail"
                REPORT_FAIL=1
            fi
        fi
    fi    
fi

# General cleanup of test stuff we may have had to keep for the report
if ([ "${LOCAL_BUILD}" == "0" ] || [ "${TEST_FAIL}" == 0 ]) && \
    [ ! "${KEEP_OUTPUT}" == "1" ] && \
    [ -z "${SAVE_WORK_DIR}" ] && \
    [ "${DO_TEST}" != "0" ] ; then
    # On anything other than a failed local run, and if we haven't been told not to,
    # and if we ran tests, clean up the test output.
    # We only ever clean the default work path
    rm -rf "./vgci_work"
fi
if [ ! "${REUSE_VENV}" == "1" ]; then
    # If we aren't re-using the virtualenv, clean it up
    rm -rf .env
fi

if [ -d "/mnt/ephemeral" ]
then
    rm -rf $TMPDIR
fi

# Decide an exit status: use the first failing stage
if [ "${BUILD_FAIL}" != "0" ]
then
    echo "Build phase has failed"
    exit "${BUILD_FAIL}"
fi

if [ "${TEST_FAIL}" != "0" ]
then
    echo "Test phase has failed"
    exit "${TEST_FAIL}"
fi

if [ "${REPORT_FAIL}" != "0" ]
then
    echo "Report phase has failed"
    exit "${REPORT_FAIL}"
fi

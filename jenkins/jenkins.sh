#!/bin/bash

# Run some CI tests on vg using toil-vg

# https://github.com/BD2KGenomics/toil-vg

# This script is hooked into 

# http://jenkins.cgcloud.info

# Most of the setup here is cribbed from other cgcloud jenkins projects such as toil-vg
# itself

# Note: we assume we run this in vg/ ie inside the vg directory we want to test

#!/bin/bash

# Should we build and run locally, or should we use Docker?
LOCAL_BUILD=0
# Should we re-use and keep around the same virtualenv?
REUSE_VENV=0
# Should we keep our test output around after uploading the new baseline?
KEEP_OUTPUT=0
# Should we show stdout and stderr from tests? If so, set to "-s".
SHOW_OPT=""
# What toil-vg should we install?
TOIL_VG_PACKAGE="git+https://github.com/vgteam/toil-vg.git@5f6756485ecd3cc36d4f91538b56f75be26f50eb"
# What tests should we run?
# Should be something like "jenkins/vgci.py::VGCITest::test_sim_brca2_snp1kg"
PYTEST_TEST_SPEC="jenkins/vgci.py"

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] \n"
    printf "Options:\n\n"
    printf "\t-l\t\tBuild vg locally (instead of in Docker) and don't use Docker at all.\n"
    printf "\t\t\tNon-Python dependencies must be installed.\n"
    printf "\t-r\t\tRe-use virtualenvs across script invocations. \n"
    printf "\t-k\t\tKeep on-disk output. \n"
    printf "\t-s\t\tShow test output and error streams (pass -s to pytest). \n"
    printf "\t-p PACKAGE\tUse the given Python package specifier to install toil-vg.\n"
    printf "\t-t TESTSPEC\tUse the given PyTest test specifier to select tests to run.\n"
    exit 1
}

while getopts "lrksp:t:" o; do
    case "${o}" in
        l)
            LOCAL_BUILD=1
            ;;
        r)
            REUSE_VENV=1
            ;;
        k)
            KEEP_OUTPUT=1
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
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [ ! -e ~/.aws/credentials ]; then
    >&2 echo "WARNING: No AWS credentials at ~/.aws/credentials; test data may not be able to be downloaded!"
fi

# Most of the script, we want to die on error
set -e

# Maximum number of minutes that can have passed since new vg docker image built
NUM_CORES=`cat /proc/cpuinfo | grep "^processor" | wc -l`
# Create Toil venv
if [ ! "${REUSE_VENV}" == "1" ]; then
    rm -rf .env
fi
if [ ! -e .env ]; then
    virtualenv  .env
fi
. .env/bin/activate

# Prepare directory for temp files (assuming cgcloud file structure)
if [ -d "/mnt/ephemeral" ]
then
     TMPDIR=/mnt/ephemeral/tmp
     rm -rf $TMPDIR
     mkdir $TMPDIR
     export TMPDIR
fi

# Upgrade pip so that it can use the wheels for numpy & scipy, so that they
# don't try to build from source
pip install --upgrade pip

# Create s3am venv
if [ ! "${REUSE_VENV}" == "1" ]; then
    rm -rf s3am
fi
if [ ! -e s3am ]; then
    virtualenv --never-download s3am && s3am/bin/pip install s3am==2.0
fi
mkdir -p bin
# Expose binaries to the PATH
ln -snf ${PWD}/s3am/bin/s3am bin/
export PATH=$PATH:${PWD}/bin

# Create awscli venv
if [ ! "${REUSE_VENV}" == "1" ]; then
    rm -rf awscli
fi
if [ ! -e awscli ]; then
    virtualenv --never-download awscli && awscli/bin/pip install awscli
fi
# Expose binaries to the PATH
ln -snf ${PWD}/awscli/bin/aws bin/
export PATH=$PATH:${PWD}/bin

# Dependencies for running tests.  Need numpy, scipy and sklearn
# for running toil-vg mapeval, and dateutils and reqests for ./mins_since_last_build.py
pip install numpy
pip install scipy
pip install sklearn
pip install dateutils
pip install requests
pip install timeout_decorator
pip install pytest
pip install pygithub
pip install toil[aws,mesos]
# Don't manually install boto since toil just installs its preferred version

# Install toil-vg itself
echo "Installing toil-vg from ${TOIL_VG_PACKAGE}"
pip install --upgrade "${TOIL_VG_PACKAGE}"
if [ "$?" -ne 0 ]
then
    echo "pip install toil-vg fail"
    exit 1
fi

# Make sure we have submodules
git submodule update --init --recursive

# we pass some parameters through pytest by way of our config file
# in particular, we set the vg version and cores, and specify
# that we want to keep all the results in vgci-work/
printf "cores ${NUM_CORES}\n" > vgci_cfg.tsv
printf "teardown False\n" >> vgci_cfg.tsv
printf "workdir ./vgci-work\n" >> vgci_cfg.tsv
#printf "verify False\n" >> vgci_cfg.tsv
#printf "baseline ./vgci-baseline\n" >> vgci_cfg.tsv

rm -rf vgci-work
mkdir vgci-work

if [ "${LOCAL_BUILD}" == "1" ]
then
    # Just build vg here
    . ./source_me.sh
    make -j ${NUM_CORES}

    if [ "$?" -ne 0 ]
    then
        echo "vg local build fail"
        exit 1
    fi
    VG_VERSION=`vg version`
    printf "vg-docker-version None\n" >> vgci_cfg.tsv
    printf "container None\n" >> vgci_cfg.tsv
else
    # Build a docker image locally.  Can be useful when don't
    # have priveleges to easily install dependencies

    # we actually want to throw git in our local image so we can get
    # a proper version
    rm -f .dockerignore

    docker pull ubuntu:16.04
    DOCKER_TAG="jenkins-docker-vg-local"
    docker build --no-cache -t "jenkins-docker-vg-local" -f jenkins/Dockerfile.jenkins .
    if [ "$?" -ne 0 ]
    then
        echo "vg docker build fail"
        exit 1
    fi
    VG_VERSION=`docker run jenkins-docker-vg-local vg version`
    printf "vg-docker-version jenkins-docker-vg-local\n" >> vgci_cfg.tsv
fi

# For the actual test and the cleanup, continue on error
set +e

# run the tests, output the junit report for Jenkins
pytest -vv "${PYTEST_TEST_SPEC}" --junitxml=test-report.xml ${SHOW_OPT}
PYRET="$?"

# Generate a report in two files: HTML full output, and a Markdown summary.
# Takes as input the Jenkins test result XML and the work directory with the
# test output files.
jenkins/mine-logs.py test-report.xml vgci-work/ report-html/ summary.md

# Put the report on Github for the current pull request or commit.
jenkins/post-report report-html summary.md


if [ ! -z "${BUILD_NUMBER}" ]
then
    # We are running on Jenkins (and not manually running the Jenkins tests), so
    # we probably have AWS credentials and can upload stuff to S3.

    # we publish the results to the archive
    tar czf "${VG_VERSION}_output.tar.gz" vgci-work test-report.xml jenkins/vgci.py jenkins/jenkins.sh vgci_cfg.tsv
    aws s3 cp --acl public-read "${VG_VERSION}_output.tar.gz" s3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_output_archives/

    # if success and we're merging the PR (and not just testing it), we publish results to the baseline
    if [ "$PYRET" -eq 0 ] && [ -z ${ghprbActualCommit} ]
    then
        echo "Tests passed. Updating baseline"
        aws s3 sync --acl public-read ./vgci-work/ s3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_regression_baseline
        printf "${VG_VERSION}\n" > vg_version_${VG_VERSION}.txt
        printf "${ghprbActualCommitAuthor}\n${ghprbPullTitle}\n${ghprbPullLink}\n" >> vg_version_${VG_VERSION}.txt
        aws s3 cp --acl public-read vg_version_${VG_VERSION}.txt s3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_regression_baseline/
    fi
fi
    

# clean up changes to bin
# Don't disturb bin/protoc or vg will want to rebuild protobuf needlessly 
rm bin/aws bin/s3am

if [ ! "${REUSE_VENV}" == "1" ]; then
    rm -rf awscli s3am
fi

if ([ "${LOCAL_BUILD}" == "0" ] || [ "${PYRET}" == 0 ]) && [ ! "${KEEP_OUTPUT}" == "1" ]; then
    # On anything other than a failed local run, and if we haven't been told not to, clean up the test output.
    rm -rf vgci-work
fi
if [ ! "${REUSE_VENV}" == "1" ]; then
    # If we aren't re-using the virtualenv, clean it up
    rm -rf .env
fi

if [ -d "/mnt/ephemeral" ]
then
    rm -rf $TMPDIR
fi

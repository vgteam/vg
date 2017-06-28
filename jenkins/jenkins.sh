#!/bin/sh

# Run some CI tests on vg using toil-vg

# https://github.com/BD2KGenomics/toil-vg

# This script is hooked into 

# http://jenkins.cgcloud.info

# Most of the setup here is cribbed from other cgcloud jenkins projects such as toil-vg
# itself

# Note: we assume we run this in vg/ ie inside the vg directory we want to test

#!/bin/bash

usage() { printf "Usage: $0 [Options] \nOptions:\n\t-l\t Build vg locally (instead in Docker). Non-python dependencies must be installed\n" 1>&2; exit 1; }

while getopts "l" o; do
    case "${o}" in
        l)
            LOCAL_BUILD=1
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

# Maximum number of minutes that can have passed since new vg docker image built
NUM_CORES=`cat /proc/cpuinfo | grep "^processor" | wc -l`
# Create Toil venv
rm -rf .env
virtualenv  .env
. .env/bin/activate

# Prepare directory for temp files (assuming cgcloud file structure)
if [ -d "/mnt/ephemeral" ]
then
     TMPDIR=/mnt/ephemeral/tmp
     rm -rf $TMPDIR
     mkdir $TMPDIR
     export TMPDIR
fi

# Create s3am venv
rm -rf s3am
virtualenv --never-download s3am && s3am/bin/pip install s3am==2.0
mkdir -p bin
# Expose binaries to the PATH
ln -snf ${PWD}/s3am/bin/s3am bin/
export PATH=$PATH:${PWD}/bin

# Create awscli venv
rm -rf awscli
virtualenv --never-download awscli && awscli/bin/pip install awscli
# Expose binaries to the PATH
ln -snf ${PWD}/awscli/bin/aws bin/
export PATH=$PATH:${PWD}/bin

# Dependencies for running tests.  Need numpy, scipy and sklearn
# for running toil-vg mapeval, and dateutils and reqests for ./mins_since_last_build.py
pip install numpy scipy sklearn dateutils requests timeout_decorator pytest boto

# Install toil-vg itself
pip install toil[aws,mesos] "toil-vg==1.2.1a1.dev409"
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

if [ -n "${LOCAL_BUILD}" ]
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

# run the tests, output the junit report for Jenkins
pytest -vv jenkins/vgci.py --junitxml=test-report.xml
PYRET="$?"

# we publish the results to the archive
tar czf "${VG_VERSION}_output.tar.gz" vgci-work test-report.xml jenkins/vgci.py jenkins/jenkins.sh vgci_cfg.tsv
aws s3 cp "${VG_VERSION}_output.tar.gz" s3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_output_archives/

# if success and we're merging the PR, we publish results to the baseline
if [ "$PYRET" -eq 0 ] && [[ ${ghprbActualCommit+x} ]]
then
    echo "Tests passed. Updating baseline"
    aws s3 sync ./vgci-work/ s3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_regression_baseline
    printf "${VG_VERSION}\n" > vg_version_${VG_VERSION}.txt
    printf "${ghprbActualCommitAuthor}\n${ghprbPullTitle}\n${ghprbPullLink}\n" >> vg_version_${VG_VERSION}.txt
    aws s3 cp vg_version_${VG_VERSION}.txt s3://cgl-pipeline-inputs/vg_cgl/vg_ci/jenkins_regression_baseline/
fi

# clean working copy to satisfy corresponding check in Makefile
rm -rf bin awscli s3am

rm -rf .env vgci-work
if [ -d "/mnt/ephemeral" ]
then
    rm -rf $TMPDIR
fi


# Run some CI tests on the latest vg docker image

# https://quay.io/repository/vgteam/vg?tag=latest&tab=tags

# using toil-vg

# https://github.com/BD2KGenomics/toil-vg

# This script is hooked into 

# http://jenkins.cgcloud.info

# and is run after each push to this repository (vg_docker), after a time delay which
# allows Travis to do its thing and build the docker image first.  A cleaner alternative
# would be to have Travis trigger Jenkins at the end of ./build.sh but it seems easier
# now to be able to piggy back on the GitHub hooks where credentials are already set up.

# Most of the setup here is cribbed from other cgcloud jenkins projects such as toil-vg
# itself

#!/bin/bash

usage() { printf "Usage: $0 [Options] \nOptions:\n\t-b <B>\t Build vg for branch B locally\n\t-r <R>\t Build docker image for vg revision R locally\n" 1>&2; exit 1; }

while getopts "b:r:" o; do
    case "${o}" in
        b)
            BRANCH=$OPTARG
            ;;
        r)
            REVISION=$OPTARG
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

# Maximum number of minutes that can have passed since new vg docker image built
MAX_MINUTES_ELAPSED=6000000
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
pip install toil[aws,mesos] toil-vg

# we pass some parameters through pytest by way of our config file
# in particular, we set the vg version and cores, and specify
# that we want to keep all the results in vgci-work/
printf "cores ${NUM_CORES}\n" > vgci_cfg.tsv
printf "teardown False\n" >> vgci_cfg.tsv
printf "workdir ./vgci-work\n" >> vgci_cfg.tsv
#printf "verify False\n" >> vgci_cfg.tsv
#printf "baseline ./vgci-baseline\n" >> vgci_cfg.tsv

# Build a local docker image for given revision. 
if [ -n "${REVISION}" ]
then
    git clone https://github.com/vgteam/vg_docker.git ${TMPDIR}/vg_docker
    docker pull ubuntu:16.04
    # hack to work on jenkins image
    grep -v locale ${TMPDIR}/vg_docker/Dockerfile.build > Dockerfile.nl
    docker build --no-cache --build-arg "vg_git_revision=${REVISION}" -t "jenkins-local-${REVISION}-build" - < Dockerfile.nl
    if [ "$?" -ne 0 ]
    then
	echo "vg local docker build fail"
	exit 1
    fi

    printf "vg-docker-version jenkins-local-${REVISION}-build\n" >> vgci_cfg.tsv
    
    
# Build a local vg for given branch.  Note this requires gcc 4.9 which docker-jenkins-image
# doesn't have.  It's hard to run this from a normal ssh (which only gives sudo apt-get),
# not much hope until image upgraded.
elif [ -n "${BRANCH}" ]
then
    # Make sure we're running gcc 4.9 by def
    # https://gist.github.com/ibogun/ec0a4005c25df57a1b9d
    #sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    #sudo apt-get -qq update
    #sudo apt-get install -y gcc-4.9 g++-4.9
    #sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 50
    #sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 50
    
    # (copied from vg docker.  may be slight overkill here.  only do on jenkins instances)
    sudo apt-get -qq install -y pkg-config sudo curl pv wget pigz unzip bsdmainutils \
build-essential make automake cmake libtool bison flex git zlib1g-dev libbz2-dev libncurses5-dev \
libgoogle-perftools-dev libjansson-dev librdf-dev jq bc rs redland-utils raptor2-utils \
rasqal-utils samtools

    git clone https://github.com/vgteam/vg.git --branch $BRANCH --recursive ${TMPDIR}/vg.local
    pushd ${TMPDIR}/vg.local
    . source_me.sh	 
    make -j ${NUM_CORES} ; make
    if [ "$?" -ne 0 ]
    then
	echo "vg make fail"
	exit 1
    fi
    VG_VERSION=`vg version`
    popd	 
    # disable docker (which is on bydefault)
    printf "vg-docker-version None\n" >> vgci_cfg.tsv

# if no branch specified, we look for a new docker image on quay
else
	# We only proceed if we have a new docker image to use
	QUAY_TAG=`python ./quay_tag_info.py vgteam/vg --max-age $MAX_MINUTES_ELAPSED`
	if [ "$?" -eq 0 ] && [ "${#QUAY_TAG}" -ge 10 ]
	then
		 VG_VERSION=`docker run ${QUAY_TAG} vg version`
		 printf "vg-docker-version ${QUAY_TAG}\n" >> vgci_cfg.tsv
	else
		 echo "Could not find vg docker image younger than ${MAX_MINUTES_ELAPSED} minutes"
		 exit 1
	fi
fi

# run the tests, output the junit report for Jenkins
pytest -vv ./vgci.py --junitxml=test-report.xml
PYRET="$?"

# we publish the results to the archive
tar czf "${VG_VERSION}_output.tar.gz" ./vgci-work ./test-report.xml ./vgci.py ./jenkins.sh ./vgci_cfg.tsv
aws s3 cp "${VG_VERSION}_output.tar.gz" s3://glennhickey-vgci-output/

# if success, we publish results to the baseline
if [ "$PYRET" -eq 0 ]
then
    aws s3 sync ./vgci-work/ s3://glennhickey-vg-regression-baseline
fi


# clean working copy to satisfy corresponding check in Makefile
rm -rf bin awscli s3am

rm -rf .env vgci-work
if [ -d "/mnt/ephemeral" ]
then
	 rm -rf $TMPDIR
fi

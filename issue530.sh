#!/bin/bash
set -e -o pipefail
HERE="${BASH_SOURCE[0]}"
while [ -h "$HERE" ] ; do HERE="$(readlink "$HERE")"; done
HERE="$( cd -P "$( dirname "$HERE" )" && pwd )"

if [ ! -e "/tmp/TSNano_1lane_L008_13801_NA12878_R1_001.18.gam" ]; then
    wget -P /tmp https://dl.dnanex.us/F/D/YBvZXJ5k2B6jJZzb31y2KyggK440Qb3zY3Gg4Zk3/TSNano_1lane_L008_13801_NA12878_R1_001.18.gam
fi
for t in 0 1 2 4 8 16; do
    echo -n "t = $t: "
    /usr/bin/time $HERE/bin/vg count -t $t /tmp/TSNano_1lane_L008_13801_NA12878_R1_001.18.gam
done

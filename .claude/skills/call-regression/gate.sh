#!/usr/bin/env bash
# Run the whole gate for one vg binary. Call both arms at once, then score them and
# compare them with the baseline.
#
#   gate.sh <vg binary> <label>
#
# Runs go to $REGRESS_DIR/runs/<label>/<arm>. Exits 1 if any output or score differs
# from the baseline.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/config.sh"
VG=${1:?vg binary}
LABEL=${2:?label for this run, such as a commit}

for arm in short ont; do
    "$HERE/call.sh" "$VG" "$arm" "$REGRESS_DIR/runs/$LABEL/$arm" &
done
wait
for arm in short ont; do
    "$HERE/score.sh" "$REGRESS_DIR/runs/$LABEL/$arm" 2> "$REGRESS_DIR/runs/$LABEL/$arm/score.log" &
done
wait
status=0
for arm in short ont; do
    "$HERE/compare.sh" "$REGRESS_DIR/runs/$LABEL/$arm" "$arm" || status=1
done
exit $status

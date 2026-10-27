#!/usr/bin/env bash
# Check a run's output files against the baseline's, by md5.
#
#   compare.sh <run directory from call.sh> <arm: short|ont>
#
# The anchors file names the vg version that wrote it, so its #vg-version line is left
# out of its hash. Exits 1 if any file differs.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=${1:?run directory}
ARM=${2:?arm: short or ont}

status=0
while read -r arm file want; do
    [ "$arm" = "$ARM" ] || continue
    if [ "$file" = anchors.tsv ]; then
        got=$(grep -v '^#vg-version' "$RUN/$file" | md5sum | cut -d' ' -f1)
    else
        got=$(md5sum < "$RUN/$file" | cut -d' ' -f1)
    fi
    if [ "$got" = "$want" ]; then
        echo "$ARM $file: identical to the baseline"
    else
        echo "$ARM $file: DIFFERS from the baseline ($got, expected $want)"
        status=1
    fi
done < <(grep -v '^#' "$HERE/expected.tsv")

# Compare the scores too, if score.sh has run. Identical calls give identical scores, so
# this matters when the calls differ. It shows whether accuracy moved, and which way.
if [ -s "$RUN/score/summary.json" ]; then
    if cmp -s "$RUN/score/summary.json" "$HERE/expected-score-$ARM.json"; then
        echo "$ARM scores: identical to the baseline"
    else
        echo "$ARM scores: DIFFER from the baseline (baseline <, this run >):"
        diff <(python3 -I -m json.tool "$HERE/expected-score-$ARM.json") \
             <(python3 -I -m json.tool "$RUN/score/summary.json") || true
        status=1
    fi
fi
exit $status

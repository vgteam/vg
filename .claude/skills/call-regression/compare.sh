#!/usr/bin/env bash
# Check a run's output files against the baseline's, by md5.
#
#   compare.sh <run directory from call.sh> <arm: short|ont>
#
# The baseline is per contig, and CONTIG (default chr20) says which to compare with. The
# mosaic and anchors headers name the vg version that wrote them and the paths of the graph and
# reads, so those header lines are left out of their hashes. Exits 1 if any file differs.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=${1:?run directory}
ARM=${2:?arm: short or ont}
CONTIG=${CONTIG:-chr20}
EXPECTED_SCORE="$HERE/expected-score-$CONTIG-$ARM.json"
grep -q "^$CONTIG	$ARM	" "$HERE/expected.tsv" \
    || { echo "no baseline for $CONTIG $ARM in expected.tsv" >&2; exit 1; }

status=0
while read -r contig arm file want; do
    [ "$contig" = "$CONTIG" ] && [ "$arm" = "$ARM" ] || continue
    if [ "$file" = anchors.tsv ] || [ "$file" = mosaic.tsv ]; then
        got=$(grep -v -E '^#(vg-version|graph|reads)'$'\t' "$RUN/$file" | md5sum | cut -d' ' -f1)
    else
        got=$(md5sum < "$RUN/$file" | cut -d' ' -f1)
    fi
    if [ "$got" = "$want" ]; then
        echo "$CONTIG $ARM $file: identical to the baseline"
    else
        echo "$CONTIG $ARM $file: DIFFERS from the baseline ($got, expected $want)"
        status=1
    fi
done < <(grep -v '^#' "$HERE/expected.tsv")

# Compare the scores too, if score.sh has run. Identical calls give identical scores, so
# this matters when the calls differ. It shows whether accuracy moved, and which way.
if [ -s "$RUN/score/summary.json" ]; then
    if cmp -s "$RUN/score/summary.json" "$EXPECTED_SCORE"; then
        echo "$CONTIG $ARM scores: identical to the baseline"
    else
        echo "$CONTIG $ARM scores: DIFFER from the baseline (baseline <, this run >):"
        diff <(python3 -I -m json.tool "$EXPECTED_SCORE") \
             <(python3 -I -m json.tool "$RUN/score/summary.json") || true
        status=1
    fi
fi
exit $status

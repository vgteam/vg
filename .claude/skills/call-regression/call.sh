#!/usr/bin/env bash
# Run one arm of the gate with one vg binary.
#
#   call.sh <vg binary> <arm: short|ont> <output directory>
#
# Writes calls.vcf, mosaic.tsv, anchors.tsv, stderr.log (vg's stderr followed by
# /usr/bin/time -v's report) and version.txt to the output directory.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/config.sh"

VG=$(realpath "${1:?vg binary}")
ARM=${2:?arm: short or ont}
OUT=${3:?output directory}
case $ARM in short|ont) ;; *) echo "unknown arm $ARM" >&2; exit 1 ;; esac

SUB="$REGRESS_DIR/data/$ARM/${CONTIG}_0_${CONTIG}.gbz"
GAF="$REGRESS_DIR/data/$ARM/$CONTIG.gaf.gz"
[ -s "$GAF.tbi" ] || { echo "missing $GAF.tbi: run prep.sh first" >&2; exit 1; }

mkdir -p "$OUT"
"$VG" version > "$OUT/version.txt"
# shellcheck disable=SC2046
/usr/bin/time -v "$VG" call "$SUB" -p "$REF_PATH" -s "$SAMPLE" -t "$THREADS" \
    $(args_for "$ARM") \
    --gaf-reads "$GAF" --gaf-index "$GAF.tbi" \
    --anchors-out "$OUT/anchors.tsv" --mosaic-out "$OUT/mosaic.tsv" \
    > "$OUT/calls.vcf" 2> "$OUT/stderr.log"
echo "$ARM with $VG: $(grep -vc '^#' "$OUT/calls.vcf") records," \
    "$(grep 'Elapsed (wall' "$OUT/stderr.log" | awk '{print $NF}') wall," \
    "$(grep 'Maximum resident' "$OUT/stderr.log" | awk '{printf "%.1f GiB", $NF/2^20}') peak"

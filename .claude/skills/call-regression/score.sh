#!/usr/bin/env bash
# Score one run's calls against GIAB HG002 v5.0q, inside the truth's benchmark regions.
#
#   score.sh <run directory from call.sh>
#
# This scores the way the vgcall-wg-benchmark-wf workflow of user shnegi does. We split multiallelic
# records and keep only PASS records. aardvark scores the records under 50 bp against the
# small-variant truth. truvari scores all records against the structural truth, with
# bench and then refine. whatshap compares the phasing with the truth's.
#
# Writes score/ in the run directory, and prints a summary.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/config.sh"
RUN=$(realpath "${1:?run directory}")
S="$RUN/score"
T="$REGRESS_DIR/truth/$CONTIG"
REF="$REGRESS_DIR/data/short/$CONTIG.fa"
mkdir -p "$S"

bgzip -c "$RUN/calls.vcf" > "$S/calls.vcf.gz"
bcftools norm -m -any -f "$REF" -Ou "$S/calls.vcf.gz" \
    | bcftools view -f PASS,. -Ou \
    | bcftools sort -Oz -o "$S/split.vcf.gz" 2> /dev/null
tabix -f -p vcf "$S/split.vcf.gz"
bcftools view -i 'abs(strlen(REF)-strlen(ALT))<50' -Oz -o "$S/small.vcf.gz" "$S/split.vcf.gz"
tabix -f -p vcf "$S/small.vcf.gz"

rm -rf "$S/aardvark"
"$AARDVARK" compare --reference "$REF" --truth-vcf "$T/smvar.vcf.gz" \
    --query-vcf "$S/small.vcf.gz" --regions "$T/smvar.bed" --output-dir "$S/aardvark" \
    --threads 16 --truth-sample HG002 --query-sample "$SAMPLE" > "$S/aardvark.log" 2>&1

if [ -n "$TRUVARI" ]; then
    rm -rf "$S/truvari"
    "$TRUVARI" bench -b "$T/stvar.vcf.gz" -c "$S/split.vcf.gz" -f "$REF" -o "$S/truvari" \
        --bSample HG002 --cSample "$SAMPLE" --pick ac --passonly --includebed "$T/stvar.bed" \
        > "$S/truvari.log" 2>&1
    if [ -s "$S/truvari/candidate.refine.bed" ]; then
        # mafft leaves temporary directories in the working directory.
        (cd "$S" && "$TRUVARI" refine -f "$REF" -t 16 -u -a mafft "$S/truvari" >> "$S/truvari.log" 2>&1)
    fi
fi

if [ -n "$WHATSHAP" ]; then
    "$WHATSHAP" compare --names truth,calls --tsv-pairwise "$S/whatshap.tsv" \
        --switch-error-bed "$S/switch-errors.bed" \
        "$T/smvar.vcf.gz" "$S/calls.vcf.gz" > "$S/whatshap.log" 2>&1 || true
fi

python3 -I "$HERE/summarize.py" "$RUN"

#!/usr/bin/env bash
# Build the gate's inputs. For each arm, these are the contig's subgraph, its reads and its
# reference FASTA. Both arms share the truth sets, sliced to the contig. A rerun skips the
# steps that are done.
#
#   prep.sh <vg binary>
#
# Any vg build can cut the subgraphs. The gate's results do not depend on which.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/config.sh"
VG=$(realpath "${1:?vg binary}")

# The input files the expected results were made from, by md5. A different file gives
# different calls, so a mismatch stops here rather than at an unexplained diff. Checking
# the GAFs themselves takes several minutes, so by default we check only their indexes.
# Set CHECK_GAFS=1 to check the GAFs too.
check_md5() {
    local file=$1 want=$2
    local got
    got=$(md5sum "$file" | cut -d' ' -f1)
    [ "$got" = "$want" ] || { echo "md5 mismatch for $file: $got, expected $want" >&2; exit 1; }
}
check_md5 "$SHORT_GBZ" ddcd39170224632325fbb92d585b4f9e
check_md5 "$SHORT_GAF.tbi" c6705e9df8d03189cad0de28221e85f2
check_md5 "$ONT_GBZ" 95fa448761c4025730a406ec8bd4f1c0
check_md5 "$ONT_GAF.tbi" a4ec37d5394ee5dbd741853215468eb5
if [ "${CHECK_GAFS:-0}" = 1 ]; then
    check_md5 "$SHORT_GAF" db7028f3c1f434f75b5a236d00e7b717
    check_md5 "$ONT_GAF" 2f3596942072ad41effbc45824c70dc7
fi

# Truth files, checked against GIAB's checksum.md5 for the release.
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q
mkdir -p "$TRUTH_DIR"
for kind in smvar stvar; do
    for ext in benchmark.bed vcf.gz vcf.gz.tbi; do
        f="$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_$kind.$ext"
        [ -s "$f" ] || curl -sSf -o "$f" "$GIAB/$(basename "$f")"
    done
done
check_md5 "$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_smvar.vcf.gz" 0c8011ca31092b2bba0c6d583f0275d3
check_md5 "$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_smvar.benchmark.bed" bd87a42f92e01b2042eee9909534746f
check_md5 "$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz" 9471dbe3ab576179462e3575de68476d
check_md5 "$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed" ac99fa1b0dc787b0aa5c14a567518702

T="$REGRESS_DIR/truth/$CONTIG"
mkdir -p "$T"
for kind in smvar stvar; do
    if [ ! -s "$T/$kind.vcf.gz.tbi" ]; then
        bcftools view -r "$CONTIG" -Oz -o "$T/$kind.vcf.gz" \
            "$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_$kind.vcf.gz"
        bcftools index -t -f "$T/$kind.vcf.gz"
        awk -v c="$CONTIG" '$1 == c' "$TRUTH_DIR/HG002_CHM13v2.0_v5.0q_$kind.benchmark.bed" \
            > "$T/$kind.bed"
    fi
done

# Cut the component that holds the contig, and keep the reads whose paths stay inside it.
for arm in short ont; do
    D="$REGRESS_DIR/data/$arm"
    mkdir -p "$D"
    if [ ! -s "$D/${CONTIG}_0_${CONTIG}.gbz" ]; then
        (cd "$D" && "$VG" chunk -x "$(gbz_for $arm)" --gbz --contig "$CONTIG" -b "$CONTIG" -t 16)
    fi
    if [ ! -s "$D/$CONTIG.fa.fai" ]; then
        "$VG" paths -x "$D/${CONTIG}_0_${CONTIG}.gbz" -F -Q "$REF_PATH" \
            | sed "1s/.*/>$CONTIG/" > "$D/$CONTIG.fa"
        samtools faidx "$D/$CONTIG.fa"
    fi
    if [ ! -s "$D/$CONTIG.gaf.gz.tbi" ]; then
        "$VG" convert -f "$D/${CONTIG}_0_${CONTIG}.gbz" | awk '$1 == "S" {print $2}' \
            > "$D/$CONTIG.nodes.txt"
        RANGE=$("$VG" stats -r "$D/${CONTIG}_0_${CONTIG}.gbz" | cut -f2 | tr : -)
        tabix "$(gaf_for $arm)" "{}:$RANGE" \
            | python3 -I "$HERE/extract_reads.py" "$D/$CONTIG.nodes.txt" \
            | bgzip -@ 8 > "$D/$CONTIG.gaf.gz"
        tabix -p gaf "$D/$CONTIG.gaf.gz"
    fi
done
# Both graphs carry the same CHM13, so the two FASTAs must agree, or the truth's
# coordinates would not fit one of them.
cmp -s "$REGRESS_DIR/data/short/$CONTIG.fa" "$REGRESS_DIR/data/ont/$CONTIG.fa" \
    || { echo "the two graphs' $CONTIG reference sequences differ" >&2; exit 1; }
echo "prep done in $REGRESS_DIR"

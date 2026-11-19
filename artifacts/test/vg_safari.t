#!/usr/bin/env bash

# Command-line tests of vg rymer and vg safari, the SAFARI ancient DNA mapper. Run from test/:
#
#     cd test && prove -v ../artifacts/test/vg_safari.t

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 15

WORK="$(mktemp -d)"

vg construct -a -r small/x.fa -v small/x.vcf.gz > "${WORK}/x.vg"
vg index -x "${WORK}/x.xg" "${WORK}/x.vg"
vg gbwt -o "${WORK}/x.gbwt" -x "${WORK}/x.vg" -v small/x.vcf.gz
vg gbwt -x "${WORK}/x.xg" -g "${WORK}/x.gbz" "${WORK}/x.gbwt"
vg index -j "${WORK}/x.dist" "${WORK}/x.vg"
vg minimizer -k 29 -w 11 -d "${WORK}/x.dist" -o "${WORK}/x.min" -z "${WORK}/x.zipcodes" "${WORK}/x.gbz"

# Profiles in which C reads as T, and G as A, often at the ends of a fragment.
HEADER="A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n"
printf "${HEADER}" > "${WORK}/x_5p.prof"
printf "${HEADER}" > "${WORK}/x_3p.prof"
for RATE in 0.30 0.20 0.12 0.08 0.05 0.03 0.02 0.01 ; do
    printf "0.001\t0.001\t0.001\t0.001\t0.001\t${RATE}\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\n" >> "${WORK}/x_5p.prof"
    printf "0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t${RATE}\t0.001\t0.001\t0.001\t0.001\t0.001\n" >> "${WORK}/x_3p.prof"
done

vg rymer -k 29 -w 11 -o "${WORK}/x.ry" "${WORK}/x.gbz"
is "${?}" "0" "vg rymer builds an RYmer index"

vg rymer 2>&1 | grep -q "10.1093/nargab/lqaf170"
is "${?}" "0" "vg rymer help asks to cite the SAFARI paper"

vg safari 2>&1 | grep -q "@article{rubin2025safari"
is "${?}" "0" "vg safari help gives the BibTeX entry of the SAFARI paper"

vg safari -h 2>&1 | grep -q -- "--rymer-name"
is "${?}" "0" "vg safari full help describes the SAFARI options"

vg safari -h 2>&1 | grep -q -- "--hit-cap"
is "${?}" "0" "vg safari full help describes the Giraffe options"

SAFARI_INDEXES=(-Z "${WORK}/x.gbz" -d "${WORK}/x.dist" -m "${WORK}/x.min" -z "${WORK}/x.zipcodes")
PROFILES=(--deam-5p "${WORK}/x_5p.prof" --deam-3p "${WORK}/x_3p.prof")

vg safari "${SAFARI_INDEXES[@]}" -q "${WORK}/x.ry" "${PROFILES[@]}" -f reads/small.middle.ref.fq > "${WORK}/mapped.gam"
is "${?}" "0" "vg safari maps reads"
is "$(vg view -aj "${WORK}/mapped.gam" | wc -l)" "$(( $(wc -l < reads/small.middle.ref.fq) / 4 ))" "vg safari writes one alignment per read"

vg safari "${SAFARI_INDEXES[@]}" "${PROFILES[@]}" -f reads/small.middle.ref.fq > /dev/null 2>&1
isnt "${?}" "0" "vg safari needs an RYmer index"

vg safari "${SAFARI_INDEXES[@]}" -q "${WORK}/x.ry" -f reads/small.middle.ref.fq > /dev/null 2>&1
isnt "${?}" "0" "vg safari needs deamination profiles"

vg safari "${SAFARI_INDEXES[@]}" -q "${WORK}/x.min" "${PROFILES[@]}" -f reads/small.middle.ref.fq > /dev/null 2>&1
isnt "${?}" "0" "vg safari refuses an ordinary minimizer index as the RYmer index"

vg safari "${SAFARI_INDEXES[@]}" -q "${WORK}/x.ry" "${PROFILES[@]}" -j 2 -f reads/small.middle.ref.fq > /dev/null 2>&1
isnt "${?}" "0" "vg safari refuses a posterior threshold above 1"

# Short reads with C read as T at their starts, and G as A at their ends.
vg sim -s 7 -n 300 -l 45 -x "${WORK}/x.xg" -a > "${WORK}/sim.gam"
vg view -X "${WORK}/sim.gam" | awk 'NR % 4 == 2 {
        n = split($0, b, "")
        for (i = 1; i <= n; i++) {
            if (i <= 6 && b[i] == "C") { b[i] = "T" }
            if (i > n - 6 && b[i] == "G") { b[i] = "A" }
            if (i % 9 == 0 && b[i] == "C") { b[i] = "T" }
        }
        s = ""; for (i = 1; i <= n; i++) { s = s b[i] }
        print s; next
    } { print }' > "${WORK}/damaged.fq"

vg giraffe "${SAFARI_INDEXES[@]}" -f "${WORK}/damaged.fq" > "${WORK}/giraffe.gam"
vg safari "${SAFARI_INDEXES[@]}" -q "${WORK}/x.ry" "${PROFILES[@]}" -f "${WORK}/damaged.fq" > "${WORK}/safari.gam"
is "${?}" "0" "vg safari maps damaged reads"
GIRAFFE_ALIGNED="$(vg stats -a "${WORK}/giraffe.gam" | grep '^Total aligned' | cut -f2 -d':' | tr -d ' ')"
SAFARI_ALIGNED="$(vg stats -a "${WORK}/safari.gam" | grep '^Total aligned' | cut -f2 -d':' | tr -d ' ')"
is "$(( SAFARI_ALIGNED >= GIRAFFE_ALIGNED ))" "1" "vg safari aligns at least as many damaged reads as vg giraffe (${SAFARI_ALIGNED} vs ${GIRAFFE_ALIGNED})"

vg safari "${SAFARI_INDEXES[@]}" -q "${WORK}/x.ry" "${PROFILES[@]}" -f "${WORK}/damaged.fq" -f "${WORK}/damaged.fq" 2>&1 > /dev/null | grep -q "only added to single-end reads"
is "${?}" "0" "vg safari warns that paired reads get no RYmer seeds"

vg giraffe "${SAFARI_INDEXES[@]}" -f reads/small.middle.ref.fq > "${WORK}/giraffe_again.gam"
is "$(vg view -aj "${WORK}/giraffe_again.gam" | wc -l)" "$(( $(wc -l < reads/small.middle.ref.fq) / 4 ))" "vg giraffe still maps reads as before"

rm -rf "${WORK}"

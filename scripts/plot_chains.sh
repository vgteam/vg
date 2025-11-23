#!/usr/bin/env bash
# plot_chains.sh: plot read vs. path dotplots for Giraffe chaining dumps

# To use this script, run vg giraffe with --track-correctness or --set-refpos,
# and --show-work.
#
# It will dump a bunch of expalantion_<READNAME> firectories with chaining
# information in the current directory.
#
# Then run this script, which will make plots in those directories of every
# chain (numbered in decreasing order of score) against every path (including
# haplotype paths).
#
# Then you can view them with commands like:
#
#   imgcat explanation_S40_22229/plot-*CHM13*.png
# 
#   imgcat explanation_S40_22229/plot-chain0-*.png

for READ_DIR in explanation_* ; do
    READ_NAME="$(echo "${READ_DIR}" | cut -f2- -d'_')"
    for CHAIN_FILE in explanation_"${READ_NAME}"/chain*-dotplot*.tsv ; do
        CHAIN_NAME="$(echo "${CHAIN_FILE}" | rev | cut -f1 -d'/' | rev | grep -o "chain[0-9]*")"
        for CHAIN_CONTIG in $(cat "${CHAIN_FILE}" | cut -f1 | cut -f1 -d'-'| sort | uniq) ; do
            echo "Plot for ${CHAIN_NAME} of ${READ_NAME} on ${CHAIN_CONTIG} from ${CHAIN_FILE}"
            "$(dirname -- "${BASH_SOURCE[0]}")"/scatter.py --category_regex --categories "chain" "frag" "${CHAIN_CONTIG}" --types "line" "line" "point" --markers "o" "." "." --colors "g" "r" "b" --no_legend --x_label "${CHAIN_CONTIG}" --y_label "${READ_NAME}" --title "${READ_NAME} ${CHAIN_NAME} on ${CHAIN_CONTIG}" --save "explanation_${READ_NAME}/plot-${CHAIN_NAME}-on-${CHAIN_CONTIG}.png" <(cat "${CHAIN_FILE}" | grep "${CHAIN_CONTIG}") &
        done
    done
    wait
done

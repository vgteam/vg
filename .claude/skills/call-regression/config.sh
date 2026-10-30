# Inputs for the vg call regression gate. Sourced by the other scripts.
#
# Set any of these in the environment to override them. Only REGRESS_DIR has no default.

# Where the subgraphs, reads, truth slices, scoring tools and runs go. It needs about 5 GB
# for the data and a few hundred MB per run.
REGRESS_DIR=${REGRESS_DIR:?set REGRESS_DIR to a directory for the gate\'s data and runs}

# The contig the gate calls on.
export CONTIG=${CONTIG:-chr20}
REF_PATH="CHM13#0#${CONTIG}"
SAMPLE=HG002
# Threads per vg call. The output does not depend on it, but run times do, so keep it
# the same between runs whose times you compare.
THREADS=${THREADS:-32}

# The sampled graphs and their reads. Each GAF was mapped to its own graph, so a GAF works
# only with the GBZ of its read set. The defaults are where `shnegi`'s workflow left them
# when the gate was made. prep.sh checks the files by md5, wherever they are.
WF=/private/groups/migalab/shnegi/vg_anchors_project/vgcall-wg-benchmark-wf
SHORT_GBZ=${SHORT_GBZ:-$WF/results/sampled/novaseq-pcrfree-30x/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gbz}
SHORT_GAF=${SHORT_GAF:-$WF/results/aln/novaseq-pcrfree-30x/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gir-default.sorted.gaf.gz}
ONT_GBZ=${ONT_GBZ:-$WF/results/sampled/E821/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gbz}
ONT_GAF=${ONT_GAF:-$WF/results/aln/E821/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gir-r10.sorted.gaf.gz}

# GIAB HG002 v5.0q on CHM13v2.0, from
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/
TRUTH_DIR=${TRUTH_DIR:-$REGRESS_DIR/truth}

# Scoring tools. aardvark is required for small-variant scores; truvari and whatshap are
# optional, and their scores are skipped when they are missing.
# The baseline was scored with aardvark 1.0.0. The fallback is where `shnegi` kept a copy
# of that release's binary.
AARDVARK=${AARDVARK:-$(command -v aardvark \
    || echo /private/groups/migalab/shnegi/apps/aardvark/v1.0.0/aardvark-v1.0.0-x86_64-unknown-linux-gnu/aardvark)}
# setup_tools.sh installs them in TOOLS_VENV; they are also looked for on PATH.
TOOLS_VENV=${TOOLS_VENV:-$REGRESS_DIR/tools/venv}
TRUVARI=${TRUVARI:-$(PATH="$TOOLS_VENV/bin:$PATH" command -v truvari || true)}
WHATSHAP=${WHATSHAP:-$(PATH="$TOOLS_VENV/bin:$PATH" command -v whatshap || true)}

# The arguments of each arm, after the graph, reference path, sample and threads.
# ont: what the vgcall-wg-benchmark-wf workflow runs on E821 reads, as one contig.
# short: Benedict Paten's Illumina gate arm for PR 4990.
ARM_ont="--read-likelihood --preset ont"
ARM_short="--read-likelihood --mismap-max 0.95 --read-phasing --regenotype"

gbz_for() { case $1 in short) echo "$SHORT_GBZ" ;; ont) echo "$ONT_GBZ" ;; esac; }
gaf_for() { case $1 in short) echo "$SHORT_GAF" ;; ont) echo "$ONT_GAF" ;; esac; }
args_for() { case $1 in short) echo "$ARM_short" ;; ont) echo "$ARM_ont" ;; esac; }

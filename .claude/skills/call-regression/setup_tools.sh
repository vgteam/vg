#!/usr/bin/env bash
# Install the tools the gate needs besides vg and aardvark. truvari and whatshap go in a
# venv at $TOOLS_VENV. mafft, which truvari refine runs, and /usr/bin/time come from apt.
# A rerun skips what is installed.
#
#   setup_tools.sh
#
# The venv lives under $REGRESS_DIR by default, beside the data, so it lasts as long as
# the data does and is shared by everyone who uses that directory.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/config.sh"

if [ ! -x /usr/bin/time ] || ! command -v mafft > /dev/null; then
    if command -v apt-get > /dev/null; then
        sudo DEBIAN_FRONTEND=noninteractive apt-get install -qq -y time mafft
    else
        echo "install GNU time (as /usr/bin/time) and mafft, then rerun" >&2
        exit 1
    fi
fi
if [ ! -x "$TOOLS_VENV/bin/truvari" ] || [ ! -x "$TOOLS_VENV/bin/whatshap" ]; then
    python3 -m venv "$TOOLS_VENV"
    # The versions the baseline's scores were made with.
    "$TOOLS_VENV/bin/pip" install -q truvari==5.5.0 whatshap==2.8
fi
"$TOOLS_VENV/bin/truvari" version
"$TOOLS_VENV/bin/whatshap" --version
"$AARDVARK" --version

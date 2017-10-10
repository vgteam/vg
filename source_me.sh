# Work out where this script is.
# See <https://stackoverflow.com/a/246128> and
# <https://unix.stackexchange.com/q/4650>
VG_ROOT="$(dirname "$(readlink -f "$_")")"

export LIBRARY_PATH="${VG_ROOT}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${VG_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_INCLUDE_PATH="${VG_ROOT}/include:${LD_INCLUDE_PATH}"
export C_INCLUDE_PATH="${VG_ROOT}/include:${C_INCLUDE_PATH}"
export CPLUS_INCLUDE_PATH="${VG_ROOT}/include:${CPLUS_INCLUDE_PATH}"
export INCLUDE_PATH="${VG_ROOT}/include:${INCLUDE_PATH}"
export PATH="${VG_ROOT}/bin:${VG_ROOT}/scripts:${PATH}"
export CC=$(which gcc)
export CXX=$(which g++)

#
#  disable until file arguments work as in normal bash :(
#
# add bash autocompletion
#if test -n "$BASH_VERSION"
#then
#
#	 . ./autocomp.bash
#fi

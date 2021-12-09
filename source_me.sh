export LIBRARY_PATH=`pwd`/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=`pwd`/lib:$DYLD_LIBRARY_PATH
export LD_INCLUDE_PATH=`pwd`/include:$LD_INCLUDE_PATH
# Setting include directories via C_INCLUDE_PATH/CPLUS_INCLUDE_PATH will
# automatically get them demoted to the end of the search list even if a -I
# option is passed to try and bump them up earlier, before other -I options.
# We leave the Makefile in charge of finding all the include directories.
export CFLAGS="-I $(pwd)/include ${CFLAGS}"
export CXXFLAGS="-I $(pwd)/include -I$(pwd)/include/dynamic ${CXXFLAGS}"
export PATH=`pwd`/bin:`pwd`/scripts:"$PATH"
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

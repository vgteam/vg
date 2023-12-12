DEP_DIR:=./deps
SRC_DIR:=src
ALGORITHMS_SRC_DIR:=$(SRC_DIR)/algorithms
CONFIG_SRC_DIR:=$(SRC_DIR)/config
IO_SRC_DIR:=$(SRC_DIR)/io
SUBCOMMAND_SRC_DIR:=$(SRC_DIR)/subcommand
UNITTEST_SRC_DIR:=$(SRC_DIR)/unittest
UNITTEST_SUPPORT_SRC_DIR:=$(SRC_DIR)/unittest/support
BIN_DIR:=bin
UNITTEST_BIN_DIR:=$(BIN_DIR)/unittest
OBJ_DIR:=obj
SHARED_OBJ_DIR:=obj/pic
ALGORITHMS_OBJ_DIR:=$(OBJ_DIR)/algorithms
ALGORITHMS_SHARED_OBJ_DIR:=$(SHARED_OBJ_DIR)/algorithms
CONFIG_OBJ_DIR:=$(OBJ_DIR)/config
IO_OBJ_DIR:=$(OBJ_DIR)/io
IO_SHARED_OBJ_DIR:=$(SHARED_OBJ_DIR)/io
SUBCOMMAND_OBJ_DIR:=$(OBJ_DIR)/subcommand
UNITTEST_OBJ_DIR:=$(OBJ_DIR)/unittest
UNITTEST_SUPPORT_OBJ_DIR:=$(OBJ_DIR)/unittest/support
LIB_DIR:=lib
# INC_DIR must be a relative path
INC_DIR:=include
CWD:=$(shell pwd)
CXX ?= g++
PKG_CONFIG ?= pkg-config

SFX :=
EXE:=vg$(SFX)

all: $(BIN_DIR)/$(EXE)

# Magic dependencies (see <http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr>)
include $(wildcard $(OBJ_DIR)/*.d)
include $(wildcard $(SHARED_OBJ_DIR)/*.d)
include $(wildcard $(ALGORITHMS_OBJ_DIR)/*.d)
include $(wildcard $(ALGORITHMS_SHARED_OBJ_DIR)/*.d)
include $(wildcard $(CONFIG_OBJ_DIR)/*.d)
include $(wildcard $(IO_OBJ_DIR)/*.d)
include $(wildcard $(IO_SHARED_OBJ_DIR)/*.d)
include $(wildcard $(SUBCOMMAND_OBJ_DIR)/*.d)
include $(wildcard $(UNITTEST_OBJ_DIR)/*.d)
include $(wildcard $(UNITTEST_BIN_DIR)/*.d)

# What pkg-config-controlled system dependencies should we use compile and link flags from?
# Use PKG_CONFIG_PATH to point the build system at the right versions of these, if they aren't picked up automatically.
# We can't do this for our bundled, pkg-config-supporting dependencies (like htslib) because they won't be built yet.
PKG_CONFIG_DEPS := cairo libzstd
# These are like PKG_CONFIG_DEPS but we try to always link them statically, if possible.
# Note that we then must *always* link anything *else* that uses them statically.
# Jansson has to be in here because it has to come after libvgio, which is in the static deps.
PKG_CONFIG_STATIC_DEPS := protobuf jansson

# We don't ask for -fopenmp here because how we get it can depend on the compiler.
# We don't ask for automatic Make dependency file (*.d) generation here because
# the options we pass can interfere with similar options in dependency project.
CXXFLAGS := -O3 -Werror=return-type -ggdb -g $(CXXFLAGS)
# Keep dependency generation flags for just our own sources
DEPGEN_FLAGS := -MMD -MP

# Set include flags. All -I options need to go in here, so the first directory
# listed is genuinely searched first.
# We make our dependency install directory -isystem; this might not be
# necessary on all platforms and suppresses warnings.
# Also, pkg-config flags need to be made -isystem if our dependency install
# directory is, or they might put a system HTSlib before ours.
# Also, Protobuf produces an absurd number of these now, so we deduplicate them
# even though that's not *always* safe. See
# <https://stackoverflow.com/a/11532197> and
# <https://github.com/protocolbuffers/protobuf/issues/12998>
INCLUDE_FLAGS :=-I$(CWD)/$(INC_DIR) -isystem $(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR) -I$(CWD)/$(UNITTEST_SRC_DIR) -I$(CWD)/$(UNITTEST_SUPPORT_SRC_DIR) -I$(CWD)/$(SUBCOMMAND_SRC_DIR) -I$(CWD)/$(INC_DIR)/dynamic $(shell $(PKG_CONFIG) --cflags $(PKG_CONFIG_DEPS) $(PKG_CONFIG_STATIC_DEPS) | tr ' ' '\n' | awk '!x[$$0]++' | tr '\n' ' ' | sed 's/ -I/ -isystem /g')

# Define libraries to link vg against.
LD_LIB_DIR_FLAGS := -L$(CWD)/$(LIB_DIR)
LD_LIB_FLAGS := -lvcflib -ltabixpp -lgssw -lssw -lsublinearLS -lpthread -lncurses -lgcsa2 -lgbwtgraph -lgbwt -lkff -ldivsufsort -ldivsufsort64 -lvcfh -lraptor2 -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -lstructures -lbdsg -lxg -lsdsl -lzstd -lhandlegraph
# We omit Boost Program Options for now; we find it in a platform-dependent way.
# By default it has no suffix
BOOST_SUFFIX=""
# We define some more libraries to link against at the end, in static linking mode if possible, so we can use faster non-PIC code.
LD_STATIC_LIB_FLAGS := -lvgio $(CWD)/$(LIB_DIR)/libtabixpp.a $(CWD)/$(LIB_DIR)/libhts.a $(CWD)/$(LIB_DIR)/libdeflate.a -lz -lbz2 -llzma
# Some of our static libraries depend on libraries that may not always be avilable in static form.
LD_STATIC_LIB_DEPS := -lpthread -lm
# Use pkg-config to find dependencies.
# Always use --static so that we have the -l flags for transitive dependencies, in case we're doing a full static build.
# But only force static linking of the dependencies we want to use non-PIC code for, for speed.
LD_LIB_FLAGS += $(shell $(PKG_CONFIG) --libs --static $(PKG_CONFIG_DEPS))
LD_STATIC_LIB_FLAGS += $(shell $(PKG_CONFIG) --libs --static $(PKG_CONFIG_STATIC_DEPS))
# We also use plain LDFLAGS to point at system library directories that we want
# to propagate through to dependencies' builds.

# CMake builds that need to find OpenMP might not know about all the prefixes it could be installed into.
# So we make a list of prefixes to search for it.
OMP_PREFIXES:=/

COMPILER_ID=$(strip $(shell $(CXX) --version 2>&1))
ifeq ($(shell uname -s),Darwin)
    $(info OS is Mac)
    # Don't try and set an rpath on any dependency utilities because that's not
    # a thing and install names will work.
    LD_UTIL_RPATH_FLAGS=""

    # Homebrew installs a Protobuf that uses an Abseil that is built with C++17, so we need to build with at least C++17
    CXX_STANDARD=17

    # We may need libraries from Macports
    ifeq ($(shell if [ -d /opt/local/lib ];then echo 1;else echo 0;fi), 1)
        # Use /opt/local/lib if present when building dependencies
        LDFLAGS += -L/opt/local/lib
    endif
    ifeq ($(shell if [ -d /usr/local/lib ];then echo 1;else echo 0;fi), 1)
        # Use /usr/local/lib if present.
        LDFLAGS += -L/usr/local/lib
    endif
    ifeq ($(shell if [ -d /usr/local/include ];then echo 1;else echo 0;fi), 1)
        # Use /usr/local/include to the end of the include search path.
        # Make sure it is system level only so it comes after other -I paths.
        INCLUDE_FLAGS += -isystem /usr/local/include

        ifeq ($(shell if [ -d /usr/local/include/cairo ];then echo 1;else echo 0;fi), 1)	
            # pkg-config is not always smart enough to find Cairo's include path for us.
            # We make sure to grab its directory manually if we see it.
            INCLUDE_FLAGS += -isystem /usr/local/include/cairo
            LD_LIB_FLAGS += -lcairo
        endif
    endif

    ifndef HOMEBREW_PREFIX
        BREW_PATH=$(shell which brew 2>/dev/null)
        ifneq ($(BREW_PATH),)
            # Get prefix from Homebrew instead of environment
            HOMEBREW_PREFIX=$(shell brew --prefix)
        endif
    endif

    ifdef HOMEBREW_PREFIX
        # We need Bison from Homebrew instead of Apple's old Bison, and GNU coreutils
        export PATH:=$(HOMEBREW_PREFIX)/opt/bison/bin:$(HOMEBREW_PREFIX)/opt/coreutils/libexec/gnubin:$(PATH)
        # If we have homebrew, use Homebrew in general
        CXXFLAGS += -I$(HOMEBREW_PREFIX)/include
        LDFLAGS += -L$(HOMEBREW_PREFIX)/lib
    endif

	# We need to find Boost Program Options. It is usually
	# -lboost_program_options, except for on Macports installs of Boost where
	# it is -lboost_program_options-mt. If we were a real build system we would
	# try things until it worked. Instead, we guess.
    ifeq ($(shell if [ -f /opt/local/lib/libboost_program_options-mt.dylib ];then echo 1;else echo 0;fi), 1)
        # This is where Macports puts it, so use that name
		BOOST_SUFFIX="-mt"
    endif

    # Our compiler might be Apple clang, which doesn't have -fopenmp.
    ifneq ($(strip $(shell echo "$(COMPILER_ID)" | grep -i clang | wc -l)), 0)
        # This is Clang.
        $(info Compiler $(CXX) is Clang)

        # We need to use the hard way of getting OpenMP not bundled with the compiler.
        # The compiler only needs to do the preprocessing
        CXXFLAGS += -Xpreprocessor -fopenmp

        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            # libomp used to be globally installed in Homebrew
            $(info OMP source is Homebrew libomp global install)
            OMP_PREFIXES:=$(OMP_PREFIXES);$(HOMEBREW_PREFIX)
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            # libomp moved to these directories, recently, because it is now keg-only to not fight GCC
            $(info OMP source is Homebrew libomop keg)
            CXXFLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LDFLAGS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
            OMP_PREFIXES:=$(OMP_PREFIXES);$(HOMEBREW_PREFIX)/opt/libomp
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            # Macports installs libomp to /opt/local/lib/libomp
            $(info OMP source Macports)
            CXXFLAGS += -I/opt/local/include/libomp
            LDFLAGS += -L/opt/local/lib/libomp
            OMP_PREFIXES:=$(OMP_PREFIXES);/opt/local
        else
            $(error OMP is not available from either Homebrew or Macports)
        endif

        # We also need to link it
        LD_LIB_FLAGS += -lomp
    else
        $(info Compiler $(CXX) is GCC)
        # The compiler is (probably?) GNU GCC
        # On Mac, we need to make sure to configure it to use libc++ like
        # Clang, and not GNU libstdc++.
        # Otherwise, we won't be able to use any C++ system libraries from
        # Homebrew or Macports, which will be built against libc++.

        # See https://stackoverflow.com/q/22228208

        CXXFLAGS += -fopenmp

        # Find includes using Clang
        LIBCXX_INCLUDES := $(shell clang++ -print-search-dirs | perl -ne 's{^libraries: =(.*)}{$$1/../../../} && print')
        # Use them and libc++ and not the normal standard library
        CXXFLAGS := -isystem $(LIBCXX_INCLUDES)/include/c++/v1 -nostdinc++ -nodefaultlibs -lc -lc++ -lc++abi -lgcc_s.1 -Wl,-no_compact_unwind $(CXXFLAGS)

        # Make sure to use the right libgomp to go with libomp
        LD_LIB_FLAGS += -lomp -lgomp.1
    endif

    # We care about building only for the current machine. If we do something
    # more restrictive we can have trouble inlining parts of the standard
    # library that were built for something less restrictive. However,
    # Apple Clang does not recognize -march=native on ARM.
	ifeq ($(shell uname -m), x86_64)
		CXXFLAGS += -march=native
	endif

    # Note shared libraries are dylibs
    SHARED_SUFFIX = dylib
    # Define options to start static linking of libraries.
    # We don't actually do any static linking on Mac, so we leave this empty.
    START_STATIC =
    END_STATIC =
else
    # We are not running on OS X
    $(info OS is Linux)
    $(info Compiler $(CXX) is assumed to be GCC)

    # Linux can have some old compilers so we want to work back to C++14
    CXX_STANDARD=14

    # Set an rpath for vg and dependency utils to find installed libraries
    LD_UTIL_RPATH_FLAGS="-Wl,-rpath,$(CWD)/$(LIB_DIR)"
    LD_LIB_FLAGS += $(LD_UTIL_RPATH_FLAGS)
    # Make sure to allow backtrace access to all our symbols, even those which are not exported.
    # Absolutely no help in a static build.
    LD_LIB_FLAGS += -rdynamic

    # We want to link against the elfutils libraries
    LD_LIB_FLAGS += -ldwfl -ldw -ldwelf -lelf -lebl

    # We want to link against libatomic which the GNU C++ standard library needs.
    # See <https://github.com/nodejs/node/issues/30093> and <https://stackoverflow.com/q/30591313>
    LD_LIB_FLAGS += -latomic

    # We get OpenMP the normal way, using whatever the compiler knows about
    CXXFLAGS += -fopenmp

    ifeq ($(shell arch), x86_64)
        # We care about building for SSE4.2 only and not AVX, to have vaguely portable binaries
        CXXFLAGS += -msse4.2
    endif

    # Note shared libraries are so files
    SHARED_SUFFIX = so
    # Define options to start static linking of libraries on GNU ld.
    START_STATIC = -Wl,-Bstatic
    # Note that END_STATIC is only safe to use in a mostly-dynamic build, and has to appear or we will try to statically link secret trailing libraries.
    END_STATIC = -Wl,-Bdynamic


endif

# Set the C++ standard we are using
CXXFLAGS := -std=c++$(CXX_STANDARD) $(CXXFLAGS)

# Propagate CXXFLAGS and LDFLAGS to child makes and other build processes
export CXXFLAGS
$(info CXXFLAGS are $(CXXFLAGS))
export LDFLAGS
$(info LDFLAGS are $(LDFLAGS))

OMP_MISSING=$(strip $(shell echo \\\#include \<omp.h\> | $(CXX) $(CXXFLAGS) -x c++ -E /dev/stdin -o /dev/null 2>&1 | head -n1 | grep error | wc -l))
ifeq ($(OMP_MISSING), 1)
    $(warning OpenMP header omp.h is not available! vg will not be able to build!)
endif

# Actually set the Boost library option, with the determined suffix
LD_LIB_FLAGS += "-lboost_program_options$(BOOST_SUFFIX)"

# These libs need to come after libdw if used, because libdw depends on them
LD_LIB_FLAGS += -ldl -llzma -lbz2 -lzstd

# Sometimes we need to filter the assembler output. The assembler can run during
# ./configure scripts, compiler calls, or $(MAKE) calls (other than $(MAKE)
# install). So we just stick $(FILTER) at the end of all such commands.
ifeq ($(shell uname -s),Darwin)
    # We need to apply a filter to all our build command output. This discards
    # all the assembler warnings which can overwhelm Travis log storage.
    FILTER=2>&1 | python3 $(CWD)/scripts/filter-noisy-assembler-warnings.py
    # For the filter to work and not just swallow errors we also need to turn on
    # pipefail in the shell
    SHELL=/bin/bash -o pipefail
else
    # No filter
    FILTER=
endif

# When building statically, we need to tell the linker not to bail if it sees multiple definitions.
# libc on e.g. our Jenkins host does not define malloc as weak, so other mallocs can't override it in a static build.
# TODO: Why did this problem only begin to happen when libvw was added?
STATIC_FLAGS=-static -static-libstdc++ -static-libgcc -Wl,--allow-multiple-definition

# These are put into libvg. Grab everything except main
OBJ = $(filter-out $(OBJ_DIR)/main.o,$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp)))
SHARED_OBJ = $(patsubst $(OBJ_DIR)/%.o,$(SHARED_OBJ_DIR)/%.o,$(OBJ))

# And all the algorithms
ALGORITHMS_OBJ = $(patsubst $(ALGORITHMS_SRC_DIR)/%.cpp,$(ALGORITHMS_OBJ_DIR)/%.o,$(wildcard $(ALGORITHMS_SRC_DIR)/*.cpp))
ALGORITHMS_SHARED_OBJ = $(patsubst $(ALGORITHMS_OBJ_DIR)/%.o,$(ALGORITHMS_SHARED_OBJ_DIR)/%.o,$(ALGORITHMS_OBJ))

# These aren't put into libvg. They are linked into vg itself to communicate
# things about the platform.
# Config objects are built individually and conditionally; that's the point.
CONFIG_OBJ =

# But always build all the IO logic
IO_OBJ = $(patsubst $(IO_SRC_DIR)/%.cpp,$(IO_OBJ_DIR)/%.o,$(wildcard $(IO_SRC_DIR)/*.cpp))
IO_SHARED_OBJ = $(patsubst $(IO_OBJ_DIR)/%.o,$(IO_SHARED_OBJ_DIR)/%.o,$(IO_OBJ))

# These aren't put into libvg, but they provide subcommand implementations for the vg bianry
SUBCOMMAND_OBJ = $(patsubst $(SUBCOMMAND_SRC_DIR)/%.cpp,$(SUBCOMMAND_OBJ_DIR)/%.o,$(wildcard $(SUBCOMMAND_SRC_DIR)/*.cpp))

# These aren't put into libvg. But they do go into the main vg binary to power its self-test.
UNITTEST_OBJ = $(patsubst $(UNITTEST_SRC_DIR)/%.cpp,$(UNITTEST_OBJ_DIR)/%.o,$(wildcard $(UNITTEST_SRC_DIR)/*.cpp))

# These support the tests. Some should go into the main vg binary but some should only go into test-suite binaries.
UNITTEST_SUPPORT_OBJ = $(patsubst $(UNITTEST_SUPPORT_SRC_DIR)/%.cpp,$(UNITTEST_SUPPORT_OBJ_DIR)/%.o,$(wildcard $(UNITTEST_SUPPORT_SRC_DIR)/*.cpp))

# These are per-test-suite binaries we can build faster
UNITTEST_EXE = $(patsubst $(UNITTEST_SRC_DIR)/%.cpp,$(UNITTEST_BIN_DIR)/%,$(wildcard $(UNITTEST_SRC_DIR)/*.cpp))


RAPTOR_DIR:=deps/raptor
JEMALLOC_DIR:=deps/jemalloc
LOCKFREE_MALLOC_DIR:=deps/lockfree-malloc
SDSL_DIR:=deps/sdsl-lite
SNAPPY_DIR:=deps/snappy
GCSA2_DIR:=deps/gcsa2
GBWT_DIR:=deps/gbwt
GBWTGRAPH_DIR=deps/gbwtgraph
KFF_DIR=deps/kff-cpp-api
PROGRESS_BAR_DIR:=deps/progress_bar
FASTAHACK_DIR:=deps/fastahack
FERMI_DIR:=deps/fermi-lite
VCFLIB_DIR:=deps/vcflib
TABIXPP_DIR:=deps/tabixpp
HTSLIB_DIR:=deps/htslib
GSSW_DIR:=deps/gssw
SPARSEHASH_DIR:=deps/sparsehash
SPARSEPP_DIR:=deps/sparsepp
SHA1_DIR:=deps/sha1
DYNAMIC_DIR:=deps/DYNAMIC
SSW_DIR:=deps/ssw/src
LINLS_DIR:=deps/sublinear-Li-Stephens
STRUCTURES_DIR:=deps/structures
BACKWARD_CPP_DIR:=deps/backward-cpp
DOZEU_DIR:=deps/dozeu
ELFUTILS_DIR:=deps/elfutils
LIBDEFLATE_DIR:=deps/libdeflate
LIBVGIO_DIR:=deps/libvgio
LIBHANDLEGRAPH_DIR:=deps/libhandlegraph
LIBBDSG_DIR:=deps/libbdsg
XG_DIR:=deps/xg
MMMULTIMAP_DIR=deps/mmmultimap
IPS4O_DIR=deps/ips4o
BBHASH_DIR=deps/BBHash
MIO_DIR=deps/mio
ATOMIC_QUEUE_DIR=deps/atomic_queue

# Dependencies that go into libvg's archive
# These go in libvg but come from dependencies
DEP_OBJ =
DEP_OBJ += $(OBJ_DIR)/progress_bar.o
DEP_OBJ += $(OBJ_DIR)/sha1.o
DEP_OBJ += $(OBJ_DIR)/Fasta.o
DEP_SHARED_OBJ = $(patsubst $(OBJ_DIR)/%.o,$(SHARED_OBJ_DIR)/%.o,$(DEP_OBJ))

# These are libraries that we need to build before we link vg.
# It would be nice to dump their contents into libvg to make it stand-alone.
# But that requires fancy ar scripting.
# If you just pass them to ar it puts the library *file* in libvg where nothing can read it.
LIB_DEPS =
LIB_DEPS += $(LIB_DIR)/libsdsl.a
LIB_DEPS += $(LIB_DIR)/libssw.a
LIB_DEPS += $(LIB_DIR)/libsnappy.a
LIB_DEPS += $(LIB_DIR)/libgcsa2.a
LIB_DEPS += $(LIB_DIR)/libgbwt.a
LIB_DEPS += $(LIB_DIR)/libgbwtgraph.a
LIB_DEPS += $(LIB_DIR)/libkff.a
LIB_DEPS += $(LIB_DIR)/libhts.a
LIB_DEPS += $(LIB_DIR)/libtabixpp.a
LIB_DEPS += $(LIB_DIR)/libvcflib.a
LIB_DEPS += $(LIB_DIR)/libgssw.a
LIB_DEPS += $(LIB_DIR)/libvcfh.a
LIB_DEPS += $(LIB_DIR)/libsonlib.a
LIB_DEPS += $(LIB_DIR)/libpinchesandcacti.a
LIB_DEPS += $(LIB_DIR)/libraptor2.a
LIB_DEPS += $(LIB_DIR)/libfml.a
LIB_DEPS += $(LIB_DIR)/libsublinearLS.a
LIB_DEPS += $(LIB_DIR)/libstructures.a
LIB_DEPS += $(LIB_DIR)/libdeflate.a
LIB_DEPS += $(LIB_DIR)/libvgio.a
LIB_DEPS += $(LIB_DIR)/libhandlegraph.a
LIB_DEPS += $(LIB_DIR)/libbdsg.a
LIB_DEPS += $(LIB_DIR)/libxg.a
ifneq ($(shell uname -s),Darwin)
    # On non-Mac (i.e. Linux), where ELF binaries are used, pull in libdw which
    # backward-cpp will use.
    LIB_DEPS += $(LIB_DIR)/libdw.a
    LIB_DEPS += $(LIB_DIR)/libdwfl.a
    LIB_DEPS += $(LIB_DIR)/libdwelf.a
    LIB_DEPS += $(LIB_DIR)/libebl.a
    LIB_DEPS += $(LIB_DIR)/libelf.a
endif

# Control varialbe for address sanitizer
# Like valgrind but fast!
# You can `make clean && make jemalloc=off asan=on` to build with it.
asan = off
ifeq ($(asan),on)
	CXXFLAGS += -fsantitze=address
endif

# Control variable for allocator
# On the command line, you can `make jemalloc=off` if you definitely don't want jemalloc.
# Or you can `make jemalloc=debug` to use a version that tries to find memory errors.
jemalloc = on
ifeq ($(shell uname -s),Darwin)
	jemalloc = off
endif

# Only depend on these files for the final linking stage.	
# These libraries provide no headers to affect the vg build.	
LINK_DEPS =

ifeq ($(jemalloc),on)
    # Use jemalloc at link time
	LINK_DEPS += $(LIB_DIR)/libjemalloc.a
    # We have to use it statically or we can't get at its secret symbols.
	LD_LIB_FLAGS += $(LIB_DIR)/libjemalloc.a
	# Use the config object for jemalloc
    CONFIG_OBJ += $(CONFIG_OBJ_DIR)/allocator_config_jemalloc.o
else ifeq ($(jemalloc),debug)
    # Use jemalloc at link time
	LINK_DEPS += $(LIB_DIR)/libjemalloc_debug.a
    # We have to use it statically or we can't get at its secret symbols.
	LD_LIB_FLAGS += $(LIB_DIR)/libjemalloc_debug.a
	# Use the config object for jemalloc
    CONFIG_OBJ += $(CONFIG_OBJ_DIR)/allocator_config_jemalloc_debug.o
else
	# Use the config object for the normal allocator
    CONFIG_OBJ += $(CONFIG_OBJ_DIR)/allocator_config_system.o
endif

# common dependencies to build before all vg src files
DEPS = $(LIB_DEPS)
DEPS += $(INC_DIR)/gcsa/gcsa.h
DEPS += $(INC_DIR)/gbwt/dynamic_gbwt.h
DEPS += $(INC_DIR)/gbwtgraph/gbwtgraph.h
DEPS += $(INC_DIR)/kff_io.hpp
DEPS += $(INC_DIR)/lru_cache.h
DEPS += $(INC_DIR)/dynamic/dynamic.hpp
DEPS += $(INC_DIR)/sparsehash/sparse_hash_map
DEPS += $(INC_DIR)/sparsepp/spp.h
DEPS += $(INC_DIR)/sha1.hpp
DEPS += $(INC_DIR)/progress_bar.hpp
DEPS += $(INC_DIR)/backward.hpp
DEPS += $(INC_DIR)/dozeu/dozeu.h
DEPS += $(INC_DIR)/mmmultimap.hpp
DEPS += $(INC_DIR)/ips4o.hpp
DEPS += $(INC_DIR)/raptor2/raptor2.h
DEPS += $(INC_DIR)/BooPHF.h
DEPS += $(INC_DIR)/mio/mmap.hpp
DEPS += $(INC_DIR)/atomic_queue.h

.PHONY: clean clean-tests get-deps deps test set-path objs static static-docker docs man .pre-build .check-environment .check-git .no-git

# Aggregate all libvg deps, and exe deps other than libvg
LIBVG_DEPS = $(OBJ) $(ALGORITHMS_OBJ) $(IO_OBJ) $(DEP_OBJ) $(DEPS)
LIBVG_SHARED_DEPS = $(SHARED_OBJ) $(ALGORITHMS_SHARED_OBJ) $(IO_SHARED_OBJ) $(DEP_SHARED_OBJ) $(DEPS)
EXE_DEPS = $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(CONFIG_OBJ) $(DEPS) $(LINK_DEPS)

# We have a target we can build to do everything but link the library and executable
objs: $(LIBVG_DEPS) $(EXE_DEPS)

$(LIB_DIR)/libvg.a: $(LIBVG_DEPS)
	rm -f $@
	ar rs $@ $(OBJ) $(ALGORITHMS_OBJ) $(IO_OBJ) $(DEP_OBJ)
	
$(LIB_DIR)/libvg.$(SHARED_SUFFIX): $(LIBVG_SHARED_DEPS)
	rm -f $@
	$(CXX) -shared -o $@ $(SHARED_OBJ) $(ALGORITHMS_SHARED_OBJ) $(IO_SHARED_OBJ) $(DEP_SHARED_OBJ) $(LDFLAGS) $(LD_LIB_DIR_FLAGS) $(LD_LIB_FLAGS) $(LD_STATIC_LIB_FLAGS) $(LD_STATIC_LIB_DEPS) 

# Each test set can have its own binary, and not link everything static
$(UNITTEST_EXE): $(UNITTEST_BIN_DIR)/%: $(UNITTEST_OBJ_DIR)/%.o $(UNITTEST_SUPPORT_OBJ) $(CONFIG_OBJ) $(LIB_DIR)/libvg.$(SHARED_SUFFIX)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(UNITTEST_SUPPORT_OBJ) $(CONFIG_OBJ) $(LIB_DIR)/libvg.$(SHARED_SUFFIX) $(LDFLAGS) $(LD_LIB_DIR_FLAGS) $(LD_LIB_FLAGS) $(LD_STATIC_LIB_FLAGS) $(LD_STATIC_LIB_DEPS)

# For a normal dynamic build we remove the static build marker
$(BIN_DIR)/$(EXE): $(LIB_DIR)/libvg.a $(EXE_DEPS)
	-rm -f $(LIB_DIR)/vg_is_static
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$(EXE) $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(CONFIG_OBJ) $(LDFLAGS) $(LIB_DIR)/libvg.a $(LD_LIB_DIR_FLAGS) $(LD_LIB_FLAGS) $(START_STATIC) $(LD_STATIC_LIB_FLAGS) $(END_STATIC) $(LD_STATIC_LIB_DEPS)
# We keep a file that we touch on the last static build.
# If the vg linkables are newer than the last static build, we do a build
$(LIB_DIR)/vg_is_static: $(INC_DIR)/vg_environment_version.hpp $(OBJ_DIR)/main.o $(LIB_DIR)/libvg.a $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(CONFIG_OBJ) $(DEPS) $(LINK_DEPS)
	$(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$(EXE) $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(CONFIG_OBJ) $(LDFLAGS) $(LIB_DIR)/libvg.a $(STATIC_FLAGS) $(LD_LIB_DIR_FLAGS) $(LD_LIB_FLAGS) $(LD_STATIC_LIB_FLAGS) $(LD_STATIC_LIB_DEPS)
	-touch $(LIB_DIR)/vg_is_static

# We don't want to always rebuild the static vg if no files have changed.
# But we do need to rebuild it if files have changed.
# TODO: is there a way to query the mtimes of all the files and rebuild if they changed *or* vg isn't static?
# For now we link dynamically and then link statically, if we actually need to rebuild anything.
static: $(LIB_DIR)/vg_is_static

# Make sure to strip out the symbols that make the binary 300 MB, but leave the
# symbols perf needs for profiling.
static-docker: static scripts/*
	strip -d $(BIN_DIR)/$(EXE)
	DOCKER_BUILDKIT=1 docker build . -f Dockerfile.static -t vg

# We have system-level deps to install
# We want the One True Place for them to be in the Dockerfile.
get-deps:
	sudo apt-get install -qq -y --no-upgrade $(shell cat Dockerfile | sed -n '/^###DEPS_BEGIN###/,$${p;/^###DEPS_END###/q}' | grep -v '^ *#' | grep -v "^RUN" | tr '\n' ' ' | tr -d '\\')

# And we have submodule deps to build
deps: $(DEPS)

test: $(BIN_DIR)/$(EXE) $(LIB_DIR)/libvg.a test/build_graph $(BIN_DIR)/shuf $(BIN_DIR)/vcf2tsv $(FASTAHACK_DIR)/fastahack $(BIN_DIR)/rapper
	. ./source_me.sh && cd test && prove -v t
	. ./source_me.sh && doc/test-docs.sh

# Somebody has been polluting the test directory with temporary files that are not deleted after the tests.
# To make git status more useful, we delete everything that looks like a temporary file.
clean-test:
	cd test && rm -rf tmp && mkdir tmp && mv 2_2.mat build_graph.cpp default.mat tmp && rm -f *.* && mv tmp/* . && rmdir tmp

docs: $(SRC_DIR)/*.cpp $(SRC_DIR)/*.hpp $(ALGORITHMS_SRC_DIR)/*.cpp $(ALGORITHMS_SRC_DIR)/*.hpp $(SUBCOMMAND_SRC_DIR)/*.cpp $(SUBCOMMAND_SRC_DIR)/*.hpp $(UNITTEST_SRC_DIR)/*.cpp $(UNITTEST_SRC_DIR)/*.hpp $(UNITTEST_SUPPORT_SRC_DIR)/*.cpp
	doxygen
	echo "View documentation at: file://$(PWD)/doc/doxygen/index.html"
	
man: $(patsubst doc/asciidoc/man/%.adoc,doc/man/%.1,$(wildcard doc/asciidoc/man/*.adoc))

doc/man/%.1: doc/asciidoc/man/%.adoc
	asciidoctor -b manpage -d manpage -o $@ $<

# Hack to use gshuf or shuf as appropriate to the platform when testing
$(BIN_DIR)/shuf:
ifeq ($(shell uname -s),Darwin)
	ln -s `which gshuf` $(BIN_DIR)/shuf
else
	ln -s `which shuf` $(BIN_DIR)/shuf
endif

test/build_graph: test/build_graph.cpp $(LIB_DIR)/libvg.a $(SRC_DIR)/vg.hpp
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o test/build_graph test/build_graph.cpp $(LDFLAGS) $(LIB_DIR)/libvg.a $(LD_LIB_DIR_FLAGS) $(LD_LIB_FLAGS) $(START_STATIC) $(LD_STATIC_LIB_FLAGS) $(END_STATIC) $(FILTER)

$(LIB_DIR)/libjemalloc.a: $(JEMALLOC_DIR)/src/*.c
	+. ./source_me.sh && rm -f $(LIB_DIR)/libjemalloc_debug.* $(LIB_DIR)/libjemalloc.* $(LIB_DIR)/libjemalloc_pic.* && rm -Rf $(CWD)/$(INC_DIR)/jemalloc && cd $(JEMALLOC_DIR) && ./autogen.sh && ./configure --enable-prof --disable-libdl --prefix=`pwd` $(FILTER) && $(MAKE) clean && $(MAKE) $(FILTER) && cp -r lib/libjemalloc.a $(CWD)/$(LIB_DIR)/ && cp -r include/* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libjemalloc_debug.a: $(JEMALLOC_DIR)/src/*.c
	+. ./source_me.sh && rm -f $(LIB_DIR)/libjemalloc_debug.* $(LIB_DIR)/libjemalloc.* $(LIB_DIR)/libjemalloc_pic.* && rm -Rf $(CWD)/$(INC_DIR)/jemalloc && cd $(JEMALLOC_DIR) && ./autogen.sh && ./configure --enable-prof --disable-libdl --enable-debug --enable-fill --prefix=`pwd` $(FILTER) && $(MAKE) clean && $(MAKE) $(FILTER) && cp -r lib/libjemalloc.a $(CWD)/$(LIB_DIR)/libjemalloc_debug.a && cp -r include/* $(CWD)/$(INC_DIR)/

# Use fake patterns to tell Make that this rule generates all these files when run once.
# Here % should always match "lib" which is a common substring.
# See https://stackoverflow.com/a/19822767
$(LIB_DIR)/%sdsl.a $(LIB_DIR)/%divsufsort.a $(LIB_DIR)/%divsufsort64.a : $(SDSL_DIR)/lib/*.cpp $(SDSL_DIR)/include/sdsl/*.hpp
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cd $(SDSL_DIR) && AS_INTEGRATED_ASSEMBLER=1 BUILD_PORTABLE=1 CXXFLAGS="$(CPPFLAGS) $(CXXFLAGS)" ./install.sh $(CWD) $(FILTER)
else
	+. ./source_me.sh && cd $(SDSL_DIR) && BUILD_PORTABLE=1 CXXFLAGS="$(CPPFLAGS) $(CXXFLAGS)" ./install.sh $(CWD) $(FILTER)
endif

$(LIB_DIR)/libssw.a: $(SSW_DIR)/*.c $(SSW_DIR)/*.cpp $(SSW_DIR)/*.h
	+. ./source_me.sh && cd $(SSW_DIR) && $(MAKE) $(FILTER) && ar rs $(CWD)/$(LIB_DIR)/libssw.a ssw.o ssw_cpp.o && cp ssw_cpp.h ssw.h $(CWD)/$(INC_DIR)

# We need to hide -Xpreprocessor -fopenmp from Snappy, at least on Mac, because
# it will drop the -Xpreprocessor and keep the -fopenmp and upset Clang.
$(LIB_DIR)/libsnappy.a: $(SNAPPY_DIR)/*.cc $(SNAPPY_DIR)/*.h
	+. ./source_me.sh && cd $(SNAPPY_DIR) && ./autogen.sh && CXXFLAGS="$(filter-out -Xpreprocessor -fopenmp,$(CXXFLAGS))" ./configure --prefix=$(CWD) $(FILTER) && CXXFLAGS="$(filter-out -Xpreprocessor -fopenmp,$(CXXFLAGS))" $(MAKE) libsnappy.la $(FILTER) && cp .libs/libsnappy.a $(CWD)/lib/ && cp snappy-c.h snappy-sinksource.h snappy-stubs-public.h snappy.h $(CWD)/include/

$(INC_DIR)/gcsa/gcsa.h: $(LIB_DIR)/libgcsa2.a

$(LIB_DIR)/libgcsa2.a: $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(wildcard $(GCSA2_DIR)/*.cpp) $(wildcard $(GCSA2_DIR)/include/gcsa/*.h)
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cp -r $(GCSA2_DIR)/include/gcsa $(CWD)/$(INC_DIR)/ && cd $(GCSA2_DIR) && make directories && AS_INTEGRATED_ASSEMBLER=1 $(MAKE) lib/libgcsa2.a $(FILTER) && mv lib/libgcsa2.a $(CWD)/$(LIB_DIR)
else
	+. ./source_me.sh && cp -r $(GCSA2_DIR)/include/gcsa $(CWD)/$(INC_DIR)/ && cd $(GCSA2_DIR) && make directories && $(MAKE) lib/libgcsa2.a $(FILTER) && mv lib/libgcsa2.a $(CWD)/$(LIB_DIR)
endif

$(INC_DIR)/gbwt/dynamic_gbwt.h: $(LIB_DIR)/libgbwt.a

$(LIB_DIR)/libgbwt.a: $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(wildcard $(GBWT_DIR)/src/*.cpp) $(wildcard $(GBWT_DIR)/include/gbwt/*.h)
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cp -r $(GBWT_DIR)/include/gbwt $(CWD)/$(INC_DIR)/ && cd $(GBWT_DIR) && $(MAKE) clean && AS_INTEGRATED_ASSEMBLER=1 $(MAKE) $(FILTER) && mv lib/libgbwt.a $(CWD)/$(LIB_DIR)
else
	+. ./source_me.sh && cp -r $(GBWT_DIR)/include/gbwt $(CWD)/$(INC_DIR)/ && cd $(GBWT_DIR) && $(MAKE) clean && $(MAKE) $(FILTER) && mv lib/libgbwt.a $(CWD)/$(LIB_DIR)
endif

$(INC_DIR)/gbwtgraph/gbwtgraph.h: $(LIB_DIR)/libgbwtgraph.a

$(LIB_DIR)/libgbwtgraph.a: $(LIB_DIR)/libgbwt.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(LIB_DIR)/libhandlegraph.a $(wildcard $(GBWTGRAPH_DIR)/src/*.cpp) $(wildcard $(GBWTGRAPH_DIR)/include/gbwtgraph/*.h)
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cp -r $(GBWTGRAPH_DIR)/include/gbwtgraph $(CWD)/$(INC_DIR)/ && cd $(GBWTGRAPH_DIR) && $(MAKE) clean && AS_INTEGRATED_ASSEMBLER=1 $(MAKE) $(FILTER) && mv lib/libgbwtgraph.a $(CWD)/$(LIB_DIR)
else
	+. ./source_me.sh && cp -r $(GBWTGRAPH_DIR)/include/gbwtgraph $(CWD)/$(INC_DIR)/ && cd $(GBWTGRAPH_DIR) && $(MAKE) clean && $(MAKE) $(FILTER) && mv lib/libgbwtgraph.a $(CWD)/$(LIB_DIR)
endif

$(INC_DIR)/kff_io.hpp: $(LIB_DIR)/libkff.a

$(LIB_DIR)/libkff.a: $(KFF_DIR)/kff_io.cpp $(KFF_DIR)/kff_io.hpp.in
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cd $(KFF_DIR) && rm -Rf build && mkdir build && cd build && cmake .. && AS_INTEGRATED_ASSEMBLER=1 $(MAKE) $(FILTER) && cp kff_io.hpp $(CWD)/$(INC_DIR) && mv libkff.a $(CWD)/$(LIB_DIR)
else
	+. ./source_me.sh && cd $(KFF_DIR) && rm -Rf build && mkdir build && cd build && cmake .. && $(MAKE) $(FILTER) && cp kff_io.hpp $(CWD)/$(INC_DIR) && mv libkff.a $(CWD)/$(LIB_DIR)
endif

$(INC_DIR)/BooPHF.h: $(BBHASH_DIR)/BooPHF.h
	+cp $(BBHASH_DIR)/BooPHF.h $(CWD)/$(INC_DIR)

$(INC_DIR)/progress_bar.hpp: $(PROGRESS_BAR_DIR)/progress_bar.hpp
	+cp $(PROGRESS_BAR_DIR)/progress_bar.hpp $(CWD)/$(INC_DIR)

$(OBJ_DIR)/progress_bar.o: $(PROGRESS_BAR_DIR)/progress_bar.cpp $(PROGRESS_BAR_DIR)/*.hpp
	+. ./source_me.sh && $(CXX) -I$(FASTAHACK_DIR) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<
$(SHARED_OBJ_DIR)/progress_bar.o: $(PROGRESS_BAR_DIR)/progress_bar.cpp $(PROGRESS_BAR_DIR)/*.hpp
	+. ./source_me.sh && $(CXX) -I$(FASTAHACK_DIR) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $<

$(INC_DIR)/Fasta.h:  $(FASTAHACK_DIR)/Fasta.h
	+. ./source_me.sh && cd $(FASTAHACK_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

$(OBJ_DIR)/Fasta.o: $(FASTAHACK_DIR)/Fasta.cpp $(INC_DIR)/Fasta.h $(FASTAHACK_DIR)/fastahack
	+. ./source_me.sh && $(CXX) -I$(FASTAHACK_DIR) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $< $(FILTER)
$(SHARED_OBJ_DIR)/Fasta.o: $(FASTAHACK_DIR)/Fasta.cpp $(INC_DIR)/Fasta.h $(FASTAHACK_DIR)/fastahack
	+. ./source_me.sh && $(CXX) -I$(FASTAHACK_DIR) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $< $(FILTER)

# We have this target to clean up the old Protobuf we used to have.
# We can remove it after we no longer care about building properly on a dirty
# build from vg versions that shipped Protobuf themselves.
$(LIB_DIR)/cleaned_old_protobuf_v003: $(wildcard $(LIB_DIR)/libproto*) $(wildcard $(LIB_DIR)/pkgconfig/protobuf*)
	+rm -f $(LIB_DIR)/cleaned_old_protobuf*
	+rm -f $(LIB_DIR)/libproto* $(LIB_DIR)/pkgconfig/protobuf* $(BIN_DIR)/protoc
	+rm -Rf $(INC_DIR)/google/protobuf deps/protobuf
	+touch $(LIB_DIR)/cleaned_old_protobuf_v003

# We used to ship our own version of boost, but now we use the system version instead.
$(LIB_DIR)/cleaned_old_boost: $(wildcard $(LIB_DIR)/libboost_*) $(wildcard $(INC_DIR)/boost/*)
	+rm -f $(LIB_DIR)/libboost_*
	+rm -Rf $(INC_DIR)/boost
	+touch $(LIB_DIR)/cleaned_old_boost

# We used to build elfutils with libdebuginfod, but we now need to build
# without it.
$(LIB_DIR)/cleaned_old_elfutils:
	+rm -f $(LIB_DIR)/libelf.a $(LIB_DIR)/libebl.a $(LIB_DIR)/libdwfl.a  $(LIB_DIR)/libdwelf.a $(LIB_DIR)/libdw.a
	+touch $(LIB_DIR)/cleaned_old_elfutils

$(LIB_DIR)/libvgio.a: $(LIB_DIR)/libhts.a $(LIB_DIR)/libhandlegraph.a $(LIB_DIR)/pkgconfig/htslib.pc $(LIB_DIR)/cleaned_old_protobuf_v003 $(LIBVGIO_DIR)/CMakeLists.txt $(LIBVGIO_DIR)/src/*.cpp $(LIBVGIO_DIR)/include/vg/io/*.hpp $(LIBVGIO_DIR)/deps/vg.proto
	+rm -f $(CWD)/$(INC_DIR)/vg.pb.h $(CWD)/$(INC_DIR)/vg/vg.pb.h
	+rm -Rf $(CWD)/$(INC_DIR)/vg/io/
	+. ./source_me.sh && export CXXFLAGS="$(CPPFLAGS) $(CXXFLAGS)" && export LDFLAGS="$(LDFLAGS) $(LD_LIB_DIR_FLAGS)" && cd $(LIBVGIO_DIR) && rm -Rf CMakeCache.txt CMakeFiles *.cmake install_manifest.txt *.pb.cc *.pb.h *.a && rm -rf build-vg && mkdir build-vg && cd build-vg && PKG_CONFIG_PATH=$(CWD)/$(LIB_DIR)/pkgconfig:$(PKG_CONFIG_PATH) cmake -DCMAKE_CXX_STANDARD=$(CXX_STANDARD) -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_PREFIX_PATH="/usr;$(OMP_PREFIXES)" -DCMAKE_INSTALL_PREFIX=$(CWD) -DCMAKE_INSTALL_LIBDIR=lib .. $(FILTER) && $(MAKE) clean && VERBOSE=1 $(MAKE) $(FILTER) && $(MAKE) install

$(LIB_DIR)/libhandlegraph.a: $(LIBHANDLEGRAPH_DIR)/src/include/handlegraph/*.hpp $(LIBHANDLEGRAPH_DIR)/src/*.cpp
	+. ./source_me.sh && cd $(LIBHANDLEGRAPH_DIR) && rm -Rf build CMakeCache.txt CMakeFiles && mkdir build && cd build && CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS)" cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=$(CWD) -DCMAKE_INSTALL_LIBDIR=lib .. && $(MAKE) $(FILTER) && $(MAKE) install


# On Linux, libdeflate builds a .so.
# On Mac, it *still* builds an so, which is just a dylib with .so extension.
# On Mac we need to make sure to set the install name. We do that by renaming to dylib.
# We don't just leave it as .so because we need to deal with outdated .so files with no paths set.
$(LIB_DIR)/libdeflate.$(SHARED_SUFFIX): $(LIB_DIR)/libdeflate.a
	+cd $(LIBDEFLATE_DIR) && cp libdeflate.so $(CWD)/$(LIB_DIR)
	+touch $(CWD)/$(LIB_DIR)/libdeflate.so
ifeq ($(shell uname -s),Darwin)
	+mv $(LIB_DIR)/libdeflate.so $(LIB_DIR)/libdeflate.$(SHARED_SUFFIX)
	+install_name_tool -id $(CWD)/$(LIB_DIR)/libdeflate.$(SHARED_SUFFIX) $(LIB_DIR)/libdeflate.$(SHARED_SUFFIX)
endif

$(LIB_DIR)/libdeflate.a: $(LIBDEFLATE_DIR)/*.h $(LIBDEFLATE_DIR)/lib/*.h $(LIBDEFLATE_DIR)/lib/*/*.h $(LIBDEFLATE_DIR)/lib/*.c $(LIBDEFLATE_DIR)/lib/*/*.c
	+. ./source_me.sh && cd $(LIBDEFLATE_DIR) && V=1 $(MAKE) $(FILTER) && cp libdeflate.a $(CWD)/$(LIB_DIR) && cp libdeflate.h $(CWD)/$(INC_DIR)

# We build htslib after libdeflate so it can use libdeflate.
# We need to do some wizardry to get it to pick up the right build and target system types on modern autotools.
# We have to do a full build in order to install, to get the pkg-config file so libvgio can link against it.
# We also have to have the shared libdeflate or we will get complaints that the static one is not position independent.
# If we need either the library or the pkg-config file (which we didn't used to ship), run the whole build.
# We use a wildcard match to make sure make understands that both files come from one command run.
# See https://stackoverflow.com/a/3077254
# We also need to make sure that htslib searches itself before system paths, as
# a system path, in case another htslib is installed on the system. Some HTSlib
# headers look for the current HTSlib with <>.
$(LIB_DIR)/libhts%a $(LIB_DIR)/pkgconfig/htslib%pc: $(LIB_DIR)/libdeflate.a $(LIB_DIR)/libdeflate.$(SHARED_SUFFIX) $(HTSLIB_DIR)/*.c $(HTSLIB_DIR)/*.h $(HTSLIB_DIR)/htslib/*.h $(HTSLIB_DIR)/cram/*.c $(HTSLIB_DIR)/cram/*.h
	+. ./source_me.sh && cd $(HTSLIB_DIR) && rm -Rf $(CWD)/$(INC_DIR)/htslib $(CWD)/$(LIB_DIR)/libhts* && autoreconf -i && autoheader && autoconf || true
	+. ./source_me.sh && cd $(HTSLIB_DIR) && (./configure -n 2>&1 || true) | grep "build system type" | rev | cut -f1 -d' ' | rev >systype.txt
	+. ./source_me.sh && cd $(HTSLIB_DIR) && CFLAGS="-I$(CWD)/$(HTSLIB_DIR) -isystem $(CWD)/$(HTSLIB_DIR) -I$(CWD)/$(INC_DIR) $(CFLAGS)" LDFLAGS="$(LDFLAGS) -L$(CWD)/$(LIB_DIR) $(LD_UTIL_RPATH_FLAGS)" ./configure --with-libdeflate --disable-s3 --disable-gcs --disable-libcurl --disable-plugins --prefix=$(CWD) --host=$$(cat systype.txt) $(FILTER) && $(MAKE) clean && $(MAKE) $(FILTER) && $(MAKE) install

# Build and install tabixpp for vcflib.
$(LIB_DIR)/libtabixpp.a: $(LIB_DIR)/libhts.a $(TABIXPP_DIR)/*.cpp $(TABIXPP_DIR)/*.hpp
	+. ./source_me.sh && cd $(TABIXPP_DIR) && rm -f tabix.o libtabixpp.a && INCLUDES="-I$(CWD)/$(INC_DIR)" HTS_HEADERS="" $(MAKE) tabix.o $(FILTER) && ar rcs libtabixpp.a tabix.o
	+cp $(TABIXPP_DIR)/libtabixpp.a $(LIB_DIR) && cp $(TABIXPP_DIR)/tabix.hpp $(INC_DIR)
	+echo "Name: tabixpp" > $(LIB_DIR)/pkgconfig/tabixpp.pc
	+echo "Description: Self-packaged tabixpp" >> $(LIB_DIR)/pkgconfig/tabixpp.pc
	+echo "Version: 1.0" >> $(LIB_DIR)/pkgconfig/tabixpp.pc
	+echo "Cflags: -I$(CWD)/$(INC_DIR)" >> $(LIB_DIR)/pkgconfig/tabixpp.pc
	+echo "Libs: -L$(CWD)/$(LIB_DIR) -ltabixpp" >> $(LIB_DIR)/pkgconfig/tabixpp.pc

# Build vcflib. Install the library and headers but not binaries or man pages.
# We need to build as RelWithDebInfo to avoid vcflib using its own
# -march=native, which would conflict with the -march that comes in through
# CXXFLAGS from the vg Dockerfile.
# We also need to use the magic path hint to let CMake find Mac OpenMP.
# We need to use /usr first for CMake search or Ubuntu 22.04 will decide pybind11 is installed in / when actually it is only fully installed in /usr.
$(LIB_DIR)/libvcflib.a: $(LIB_DIR)/libhts.a $(LIB_DIR)/libtabixpp.a $(VCFLIB_DIR)/src/*.cpp $(VCFLIB_DIR)/src/*.hpp $(VCFLIB_DIR)/contrib/*/*.cpp $(VCFLIB_DIR)/contrib/*/*.h
	+rm -f $(VCFLIB_DIR)/contrib/WFA2-lib/VERSION
	+. ./source_me.sh && cd $(VCFLIB_DIR) && rm -Rf build && mkdir build && cd build && PKG_CONFIG_PATH="$(CWD)/$(LIB_DIR)/pkgconfig:$(PKG_CONFIG_PATH)" cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DZIG=OFF -DCMAKE_C_FLAGS="$(CFLAGS)" -DCMAKE_CXX_FLAGS="$(CXXFLAGS) ${CPPFLAGS}" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=$(CWD) -DCMAKE_PREFIX_PATH="/usr;$(OMP_PREFIXES)" .. && cmake --build .
	+cp $(VCFLIB_DIR)/contrib/filevercmp/*.h* $(INC_DIR)
	+cp $(VCFLIB_DIR)/contrib/fastahack/*.h* $(INC_DIR)
	+cp $(VCFLIB_DIR)/contrib/smithwaterman/*.h* $(INC_DIR)
	+cp $(VCFLIB_DIR)/contrib/intervaltree/*.h* $(INC_DIR)
	+cp $(VCFLIB_DIR)/contrib/multichoose/*.h* $(INC_DIR)
	+cp $(VCFLIB_DIR)/src/*.h* $(INC_DIR)
	+cp $(VCFLIB_DIR)/build/libvcflib.a $(LIB_DIR)

# vcflib binaries are all automatically built. We need this one.
$(BIN_DIR)/vcf2tsv: $(VCFLIB_DIR)/src/*.cpp $(VCFLIB_DIR)/src/*.h $(LIB_DIR)/libvcflib.a
	+cp $(VCFLIB_DIR)/build/vcf2tsv $(BIN_DIR)

$(FASTAHACK_DIR)/fastahack: $(FASTAHACK_DIR)/*.c $(FASTAHACK_DIR)/*.h $(FASTAHACK_DIR)/*.cpp
	+. ./source_me.sh && cd $(FASTAHACK_DIR) && $(MAKE) $(FILTER)

$(LIB_DIR)/libgssw.a: $(GSSW_DIR)/src/gssw.c $(GSSW_DIR)/src/gssw.h
	+. ./source_me.sh && cd $(GSSW_DIR) && $(MAKE) $(FILTER) && cp lib/libgssw.a $(CWD)/$(LIB_DIR)/ && cp src/gssw.h $(CWD)/$(INC_DIR)/

$(INC_DIR)/lru_cache.h: $(DEP_DIR)/lru_cache/*.h $(DEP_DIR)/lru_cache/*.cc
	+cd $(DEP_DIR)/lru_cache && cp *.h* $(CWD)/$(INC_DIR)/

# We moved the Dynamic headers so make sure to clean up the old ones.
$(INC_DIR)/dynamic/dynamic.hpp: $(DYNAMIC_DIR)/include/dynamic/*.hpp $(DYNAMIC_DIR)/include/dynamic/*/*.hpp
	+rm -Rf $(INC_DIR)/dynamic.hpp $(INC_DIR)/dynamic
	# annoyingly doesn't have an install option on the cmake, so we manually move their external dependency headers
	+cd $(CWD)/$(DYNAMIC_DIR) && rm -Rf build && mkdir -p build && cd build && export CXXFLAGS="$(CPPFLAGS) $(CXXFLAGS)" && cmake -DCMAKE_VERBOSE_MAKEFILE=ON .. && make && cp -r $(CWD)/$(DYNAMIC_DIR)/deps/hopscotch_map/include/* $(CWD)/$(INC_DIR)/
	# Do the copy of the main file last so we can tell if this recipe failed and redo it.
	# Otherwise we get dynamic.hpp without its deps
	+mkdir -p $(INC_DIR)/dynamic && cp -r $(CWD)/$(DYNAMIC_DIR)/include/dynamic/* $(INC_DIR)/dynamic/

$(INC_DIR)/sparsehash/sparse_hash_map: $(wildcard $(SPARSEHASH_DIR)/**/*.cc) $(wildcard $(SPARSEHASH_DIR)/**/*.h)
	+. ./source_me.sh && cd $(SPARSEHASH_DIR) && ./autogen.sh && LDFLAGS="$(LDFLAGS) $(LD_LIB_DIR_FLAGS)" ./configure --prefix=$(CWD) $(FILTER) && $(MAKE) $(FILTER) && $(MAKE) install

$(INC_DIR)/sparsepp/spp.h: $(wildcard $(SPARSEHASH_DIR)/sparsepp/*.h)
	+cp -r $(SPARSEPP_DIR)/sparsepp $(INC_DIR)/

#$(INC_DIR)/Variant.h
$(LIB_DIR)/libvcfh.a: $(DEP_DIR)/libVCFH/*.cpp $(DEP_DIR)/libVCFH/*.hpp
	+. ./source_me.sh && cd $(DEP_DIR)/libVCFH && $(MAKE) $(FILTER) && cp libvcfh.a $(CWD)/$(LIB_DIR)/ && cp vcfheader.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libsonlib.a: $(CWD)/$(DEP_DIR)/sonLib/C/inc/*.h $(CWD)/$(DEP_DIR)/sonLib/C/impl/*.c
	+. ./source_me.sh && cd $(DEP_DIR)/sonLib && kyotoTycoonLib="" $(MAKE) $(FILTER) && cp lib/sonLib.a $(CWD)/$(LIB_DIR)/libsonlib.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

$(LIB_DIR)/libpinchesandcacti.a: $(LIB_DIR)/libsonlib.a $(CWD)/$(DEP_DIR)/pinchesAndCacti/inc/*.h $(CWD)/$(DEP_DIR)/pinchesAndCacti/impl/*.c
	+. ./source_me.sh && cd $(DEP_DIR)/pinchesAndCacti && $(MAKE) $(FILTER) && cd $(CWD)/$(DEP_DIR)/sonLib && cp lib/stPinchesAndCacti.a $(CWD)/$(LIB_DIR)/libpinchesandcacti.a && cp lib/3EdgeConnected.a $(CWD)/$(LIB_DIR)/lib3edgeconnected.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

# When building raptor we need to make sure to pre-generate and fix up the lexer
# We also need to clear out its cmake stuff in case it found a wrong Bison and cached it.
$(LIB_DIR)/libraptor2.a: $(RAPTOR_DIR)/src/* $(wildcard $(RAPTOR_DIR)/build/*)
	which bison
	+. ./source_me.sh && cd $(RAPTOR_DIR)/build && rm -Rf CMakeCache.txt CMakeFiles CTestTestfile.cmake Makefile cmake_install.cmake src tests utils && cmake .. && rm -f src/turtle_parser.c && rm -f src/turtle_lexer.c && make turtle_lexer_tgt && make -f src/CMakeFiles/raptor2.dir/build.make src/turtle_lexer.c && sed -i.bak '/yycleanup/d' src/turtle_lexer.c && $(MAKE) $(FILTER) && cp src/libraptor2.a $(CWD)/$(LIB_DIR)
	+touch $(LIB_DIR)/libraptor2.a

# We need rapper from Raptor for the tests
$(BIN_DIR)/rapper: $(LIB_DIR)/libraptor2.a
	+cp $(RAPTOR_DIR)/build/utils/rapper $(BIN_DIR)/

# The Raptor header needs to be newer than the library.
# Mac Travis managed to get an old header with a new binary.
$(INC_DIR)/raptor2/raptor2.h: $(LIB_DIR)/libraptor2.a $(RAPTOR_DIR)/build/*
	+cd $(RAPTOR_DIR)/build && mkdir -p $(CWD)/$(INC_DIR)/raptor2 && cp src/*.h $(CWD)/$(INC_DIR)/raptor2
	+touch $(INC_DIR)/raptor2/raptor2.h

$(LIB_DIR)/libstructures.a: $(STRUCTURES_DIR)/src/include/structures/*.hpp $(STRUCTURES_DIR)/src/*.cpp $(STRUCTURES_DIR)/Makefile
	+. ./source_me.sh && cd $(STRUCTURES_DIR) && $(MAKE) clean && $(MAKE) lib/libstructures.a $(FILTER) && cp lib/libstructures.a $(CWD)/$(LIB_DIR)/ && cp -r src/include/structures $(CWD)/$(INC_DIR)/

$(INC_DIR)/sha1.hpp: $(SHA1_DIR)/sha1.hpp
	+cp $(SHA1_DIR)/*.h* $(CWD)/$(INC_DIR)/

$(INC_DIR)/backward.hpp: $(BACKWARD_CPP_DIR)/backward.hpp
	+cp $(BACKWARD_CPP_DIR)/backward.hpp $(CWD)/$(INC_DIR)/

$(INC_DIR)/simde/x86/sse4.1.h: $(DOZEU_DIR)/simde/*.h $(DOZEU_DIR)/simde/x86/*.h
	+cp -r $(DOZEU_DIR)/simde $(INC_DIR)

$(INC_DIR)/dozeu/dozeu.h: $(DOZEU_DIR)/*.h $(INC_DIR)/simde/x86/sse4.1.h
	+mkdir -p $(CWD)/$(INC_DIR)/dozeu && cp $(DOZEU_DIR)/*.h $(CWD)/$(INC_DIR)/dozeu/

$(LIB_DIR)/libebl.a: $(LIB_DIR)/libelf.a

$(LIB_DIR)/libdw.a: $(LIB_DIR)/libelf.a

$(LIB_DIR)/libdwelf.a: $(LIB_DIR)/libelf.a

$(LIB_DIR)/libdwfl.a: $(LIB_DIR)/libelf.a

# We can't build elfutils from Git without "maintainer mode".
# There are some release-only headers or something that it complains it can't find otherwise.
# We also don't do a normal make and make install here because we don't want to build and install all the elfutils binaries and libasm.
# We need to disable libdebuginfod or the static binary will try and load it at
# runtime and pull in incompatible libs it depends on on whatever system it's
# running on.
$(LIB_DIR)/libelf.a: $(ELFUTILS_DIR)/libebl/*.c $(ELFUTILS_DIR)/libebl/*.h $(ELFUTILS_DIR)/libdw/*.c $(ELFUTILS_DIR)/libdw/*.h $(ELFUTILS_DIR)/libelf/*.c $(ELFUTILS_DIR)/libelf/*.h $(ELFUTILS_DIR)/src/*.c $(ELFUTILS_DIR)/src/*.h $(LIB_DIR)/cleaned_old_elfutils
	+cd $(CWD)/$(INC_DIR)/ && rm -Rf elfutils gelf.h libelf.h dwarf.h libdwflP.h libdwfl.h libebl.h libelf.h
	+. ./source_me.sh && cd $(ELFUTILS_DIR) && autoreconf -i -f && ./configure --enable-maintainer-mode --disable-libdebuginfod --disable-debuginfod --prefix=$(CWD) $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/libelf && $(MAKE) clean && $(MAKE) libelf.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/libebl && $(MAKE) clean && $(MAKE) libebl.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/libdwfl && $(MAKE) clean && $(MAKE) libdwfl.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/libdwelf && $(MAKE) clean && $(MAKE) libdwelf.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/lib && $(MAKE) clean && $(MAKE) libeu.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/libcpu && $(MAKE) clean && $(MAKE) libcpu.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/backends && $(MAKE) clean && $(MAKE) libebl_backends.a $(FILTER)
	+. ./source_me.sh && cd $(ELFUTILS_DIR)/libdw && $(MAKE) clean && $(MAKE) libdw.a known-dwarf.h $(FILTER)
	+cd $(ELFUTILS_DIR) && mkdir -p $(CWD)/$(INC_DIR)/elfutils && cp libdw/known-dwarf.h libdw/libdw.h libebl/libebl.h libelf/elf-knowledge.h version.h libdwfl/libdwfl.h libdwelf/libdwelf.h $(CWD)/$(INC_DIR)/elfutils && cp libelf/gelf.h libelf/libelf.h libdw/dwarf.h $(CWD)/$(INC_DIR) && cp libebl/libebl.a libdw/libdw.a libdwfl/libdwfl.a libdwelf/libdwelf.a libelf/libelf.a $(CWD)/$(LIB_DIR)/

$(OBJ_DIR)/sha1.o: $(SHA1_DIR)/sha1.cpp $(SHA1_DIR)/sha1.hpp
	+$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $< $(FILTER)
$(SHARED_OBJ_DIR)/sha1.o: $(SHA1_DIR)/sha1.cpp $(SHA1_DIR)/sha1.hpp
	+$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $< $(FILTER)

$(LIB_DIR)/libfml.a: $(FERMI_DIR)/*.h $(FERMI_DIR)/*.c
	. ./source_me.sh && cd $(FERMI_DIR) && $(MAKE) $(FILTER) && cp *.h $(CWD)/$(INC_DIR)/ && cp libfml.a $(CWD)/$(LIB_DIR)/

# We don't need to hack the build to point at our htslib because sublinearLS gets its htslib from the include flags we set
$(LIB_DIR)/libsublinearLS.a: $(LINLS_DIR)/src/*.cpp $(LINLS_DIR)/src/*.hpp $(LIB_DIR)/libhts.a
	. ./source_me.sh && cd $(LINLS_DIR) && $(MAKE) clean && INCLUDE_FLAGS="-I$(CWD)/$(INC_DIR)" $(MAKE) libs $(FILTER) && cp lib/libsublinearLS.a $(CWD)/$(LIB_DIR)/ && mkdir -p $(CWD)/$(INC_DIR)/sublinearLS && cp src/*.hpp $(CWD)/$(INC_DIR)/sublinearLS/

$(LIB_DIR)/libbdsg.a: $(INC_DIR)/BooPHF.h $(LIBBDSG_DIR)/Makefile $(LIBBDSG_DIR)/bdsg/src/*.cpp $(LIBBDSG_DIR)/bdsg/include/bdsg/*.hpp $(LIBBDSG_DIR)/bdsg/include/bdsg/internal/*.hpp $(LIBBDSG_DIR)/bdsg/include/bdsg/overlays/*.hpp $(LIB_DIR)/libhandlegraph.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(INC_DIR)/sparsepp/spp.h $(INC_DIR)/dynamic/dynamic.hpp $(INC_DIR)/mio/mmap.hpp
	+. ./source_me.sh && rm -Rf $(CWD)/$(INC_DIR)/bdsg && cd $(LIBBDSG_DIR) && $(MAKE) clean && CPLUS_INCLUDE_PATH=$(CWD)/$(INC_DIR):$(CWD)/$(INC_DIR)/dynamic:$(CPLUS_INCLUDE_PATH) CXXFLAGS="$(INCLUDE_FLAGS) $(CXXFLAGS)" $(MAKE) $(FILTER) && cp lib/libbdsg.a $(CWD)/$(LIB_DIR) && cp -r bdsg/include/* $(CWD)/$(INC_DIR)

$(INC_DIR)/mio/mmap.hpp: $(MIO_DIR)/include/mio/*
	+. ./source_me.sh && cp -r $(MIO_DIR)/include/mio $(CWD)/$(INC_DIR)/

# It would be better to copy the atomic_queue directory rather than its contents, but to avoid re-writing mmmultimap...
$(INC_DIR)/atomic_queue.h: $(ATOMIC_QUEUE_DIR)/include/*
	+. ./source_me.sh && cp -r $(ATOMIC_QUEUE_DIR)/include/atomic_queue/* $(CWD)/$(INC_DIR)/

$(INC_DIR)/mmmultiset.hpp: $(MMMULTIMAP_DIR)/src/mmmultiset.hpp $(INC_DIR)/mmmultimap.hpp
$(INC_DIR)/mmmultimap.hpp: $(MMMULTIMAP_DIR)/src/mmmultimap.hpp $(MMMULTIMAP_DIR)/src/mmmultiset.hpp $(INC_DIR)/mio/mmap.hpp $(INC_DIR)/atomic_queue.h
	+. ./source_me.sh && cp $(MMMULTIMAP_DIR)/src/mmmultimap.hpp $(MMMULTIMAP_DIR)/src/mmmultiset.hpp $(CWD)/$(INC_DIR)/

$(INC_DIR)/ips4o.hpp: $(IPS4O_DIR)/ips4o.hpp $(IPS4O_DIR)/ips4o/*
	+. ./source_me.sh && cp -r $(IPS4O_DIR)/ips4o* $(CWD)/$(INC_DIR)/

# The xg repo has a cmake build system based all around external projects, and
# we need it to use our installed versions of everything instead.
# We also need to not build against GFAKluge
$(LIB_DIR)/libxg.a: $(XG_DIR)/src/*.hpp $(XG_DIR)/src/*.cpp $(INC_DIR)/mmmultimap.hpp $(INC_DIR)/ips4o.hpp $(LIB_DIR)/libhandlegraph.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(INC_DIR)/mio/mmap.hpp $(INC_DIR)/atomic_queue.h
	+rm -f $@
	+cp -r $(XG_DIR)/src/*.hpp $(CWD)/$(INC_DIR)
	+. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -DNO_GFAKLUGE -c -o $(XG_DIR)/xg.o $(XG_DIR)/src/xg.cpp $(FILTER)
	+ar rs $@ $(XG_DIR)/xg.o

# Auto-git-versioning

# We need to scope this variable here
GIT_VERSION_FILE_DEPS =
# Decide if .git exists and needs to be watched
ifeq ($(shell if [ -d .git ]; then echo present; else echo absent; fi),present)
    # If so, try and make a git version file
	GIT_VERSION_FILE_DEPS = .check-git
else
    # Just use the version file we have, if any
	GIT_VERSION_FILE_DEPS = .no-git
endif

# Build a real git version file.
# If it's not the same as the old one, replace the old one.
# If it is the same, do nothing and don't rebuild dependent targets.
.check-git:
	@echo "#define VG_GIT_VERSION \"$(shell git describe --always --tags 2>/dev/null || echo git-error)\"" > $(INC_DIR)/vg_git_version.hpp.tmp
	@diff $(INC_DIR)/vg_git_version.hpp.tmp $(INC_DIR)/vg_git_version.hpp >/dev/null 2>/dev/null || cp $(INC_DIR)/vg_git_version.hpp.tmp $(INC_DIR)/vg_git_version.hpp
	@rm -f $(INC_DIR)/vg_git_version.hpp.tmp

# Make sure the version file exists, if we weren't given one in our tarball
.no-git:
	@if [ ! -e $(INC_DIR)/vg_git_version.hpp ]; then \
		touch $(INC_DIR)/vg_git_version.hpp; \
	fi;

$(INC_DIR)/vg_git_version.hpp: $(GIT_VERSION_FILE_DEPS)
# Build an environment version file with this phony target.
# If it's not the same as the old one, replace the old one.
# If it is the same, do nothing and don't rebuild dependent targets.
.check-environment:
	@echo "#define VG_COMPILER_VERSION \"$(shell $(CXX) --version 2>/dev/null | head -n 1)\"" > $(INC_DIR)/vg_environment_version.hpp.tmp
	@echo "#define VG_OS \"$(shell uname)\"" >> $(INC_DIR)/vg_environment_version.hpp.tmp
	@echo "#define VG_BUILD_USER \"$(shell whoami)\"" >> $(INC_DIR)/vg_environment_version.hpp.tmp
	@echo "#define VG_BUILD_HOST \"$(shell hostname)\"" >> $(INC_DIR)/vg_environment_version.hpp.tmp
	@diff $(INC_DIR)/vg_environment_version.hpp.tmp $(INC_DIR)/vg_environment_version.hpp >/dev/null || cp $(INC_DIR)/vg_environment_version.hpp.tmp $(INC_DIR)/vg_environment_version.hpp
	@rm -f $(INC_DIR)/vg_environment_version.hpp.tmp

# The way to get the actual file is to maybe replace it.
$(INC_DIR)/vg_environment_version.hpp: .check-environment

###################################
## VG source code compilation begins here
####################################

$(OBJ_DIR)/version.o: $(SRC_DIR)/version.cpp $(SRC_DIR)/version.hpp $(INC_DIR)/vg_git_version.hpp $(INC_DIR)/vg_environment_version.hpp

########################
## Pattern Rules
########################

# Define a default rule for building objects from CPP files
# Depend on the .d file so we rebuild if dependency info is missing/deleted
# Make sure to touch the .o file after the compiler finishes so it is always newer than the .d file
# Use static pattern rules so the dependency files will not be ignored if the output exists
# See <https://stackoverflow.com/a/34983297>
$(OBJ) $(OBJ_DIR)/main.o: $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(SHARED_OBJ): $(SHARED_OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(SHARED_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -fPIC -c -o $@ $< $(FILTER)
	@touch $@
$(ALGORITHMS_OBJ): $(ALGORITHMS_OBJ_DIR)/%.o : $(ALGORITHMS_SRC_DIR)/%.cpp $(ALGORITHMS_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(ALGORITHMS_SHARED_OBJ): $(ALGORITHMS_SHARED_OBJ_DIR)/%.o : $(ALGORITHMS_SRC_DIR)/%.cpp $(ALGORITHMS_SHARED_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -fPIC -c -o $@ $< $(FILTER)
	@touch $@
$(IO_OBJ): $(IO_OBJ_DIR)/%.o : $(IO_SRC_DIR)/%.cpp $(IO_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(IO_SHARED_OBJ): $(IO_SHARED_OBJ_DIR)/%.o : $(IO_SRC_DIR)/%.cpp $(IO_SHARED_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -fPIC -c -o $@ $< $(FILTER)
	@touch $@
$(SUBCOMMAND_OBJ): $(SUBCOMMAND_OBJ_DIR)/%.o : $(SUBCOMMAND_SRC_DIR)/%.cpp $(SUBCOMMAND_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(UNITTEST_OBJ): $(UNITTEST_OBJ_DIR)/%.o : $(UNITTEST_SRC_DIR)/%.cpp $(UNITTEST_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(UNITTEST_SUPPORT_OBJ): $(UNITTEST_SUPPORT_OBJ_DIR)/%.o : $(UNITTEST_SUPPORT_SRC_DIR)/%.cpp $(UNITTEST_SUPPORT_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
	
# Config objects get individual rules
$(CONFIG_OBJ_DIR)/allocator_config_jemalloc.o: $(CONFIG_SRC_DIR)/allocator_config_jemalloc.cpp $(CONFIG_OBJ_DIR)/allocator_config_jemalloc.d $(DEPS) $(LIB_DIR)/libjemalloc.a
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(CONFIG_OBJ_DIR)/allocator_config_jemalloc_debug.o: $(CONFIG_SRC_DIR)/allocator_config_jemalloc_debug.cpp $(CONFIG_OBJ_DIR)/allocator_config_jemalloc_debug.d $(DEPS) $(LIB_DIR)/libjemalloc_debug.a
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@
$(CONFIG_OBJ_DIR)/allocator_config_system.o: $(CONFIG_SRC_DIR)/allocator_config_system.cpp $(CONFIG_OBJ_DIR)/allocator_config_system.d $(DEPS)
	. ./source_me.sh && $(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DEPGEN_FLAGS) -c -o $@ $< $(FILTER)
	@touch $@

# Use a fake rule to build .d files, so we don't complain if they don't exist.
$(OBJ_DIR)/%.d: ;
$(ALGORITHMS_OBJ_DIR)/%.d: ;
$(CONFIG_OBJ_DIR)/%.d: ;
$(IO_OBJ_DIR)/%.d: ;
$(SUBCOMMAND_OBJ_DIR)/%.d: ;
$(UNITTEST_OBJ_DIR)/%.d: ;

# Don't delete them.
.PRECIOUS: $(OBJ_DIR)/%.d $(ALGORITHMS_OBJ_DIR)/%.d $(CONFIG_OBJ_DIR)/%.d $(IO_OBJ_DIR)/%.d $(SUBCOMMAND_OBJ_DIR)/%.d $(UNITTEST_OBJ_DIR)/%.d

# Use no implicit rules
.SUFFIXES:

###################################
## VG source code compilation ends here
####################################


# Make directories before quitting target due to missing protoc.
# If we run the rest of the build without these, lib and include can become files.
# TODO: quitting if no protoc doesn't reliably stop the build.
.pre-build:
	@if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	@if [ ! -d $(UNITTEST_BIN_DIR) ]; then mkdir -p $(UNITTEST_BIN_DIR); fi
	@if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	@if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	@if [ ! -d $(SHARED_OBJ_DIR) ]; then mkdir -p $(SHARED_OBJ_DIR); fi
	@if [ ! -d $(ALGORITHMS_OBJ_DIR) ]; then mkdir -p $(ALGORITHMS_OBJ_DIR); fi
	@if [ ! -d $(ALGORITHMS_SHARED_OBJ_DIR) ]; then mkdir -p $(ALGORITHMS_SHARED_OBJ_DIR); fi
	@if [ ! -d $(CONFIG_OBJ_DIR) ]; then mkdir -p $(CONFIG_OBJ_DIR); fi
	@if [ ! -d $(IO_OBJ_DIR) ]; then mkdir -p $(IO_OBJ_DIR); fi
	@if [ ! -d $(IO_SHARED_OBJ_DIR) ]; then mkdir -p $(IO_SHARED_OBJ_DIR); fi
	@if [ ! -d $(SUBCOMMAND_OBJ_DIR) ]; then mkdir -p $(SUBCOMMAND_OBJ_DIR); fi
	@if [ ! -d $(UNITTEST_OBJ_DIR) ]; then mkdir -p $(UNITTEST_OBJ_DIR); fi
	@if [ ! -d $(UNITTEST_SUPPORT_OBJ_DIR) ]; then mkdir -p $(UNITTEST_SUPPORT_OBJ_DIR); fi
	@if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	@protoc --version >/dev/null 2>/dev/null || (echo "Error: protobuf compiler (protoc) not available!" ; exit 1)
	@if [ -e $(INC_DIR)/vg/vg.pb.h ] ; then \
		HEADER_VER=$$(cat $(INC_DIR)/vg/vg.pb.h | grep GOOGLE_PROTOBUF_VERSION | sed 's/[^0-9]*\([0-9]*\)[^0-9]*/\1/' | head -n1); \
		WORKDIR=$$(pwd); \
		TESTDIR=$$(mktemp -d); \
		echo 'syntax = "proto3";' > $${TESTDIR}/empty.proto; \
		protoc $${TESTDIR}/empty.proto --proto_path=$${TESTDIR} --cpp_out=$${TESTDIR}; \
		PROTOC_VER=$$(cat $${TESTDIR}/empty.pb.h | grep GOOGLE_PROTOBUF_VERSION | sed 's/[^0-9]*\([0-9]*\)[^0-9]*/\1/' | head -n1); \
		if [ "$${HEADER_VER}" != "$${PROTOC_VER}" ] ; then \
			echo "Protobuf version has changed!"; \
			echo "Headers are for $${HEADER_VER} but we make headers for $${PROTOC_VER}"; \
			echo "Need to rebuild libvgio"; \
			rm -f $(LIB_DIR)/libvgio.a; \
			rm -f $(INC_DIR)/vg/vg.pb.h; \
		fi; \
		rm $${TESTDIR}/empty.proto $${TESTDIR}/empty.pb.h $${TESTDIR}/empty.pb.cc; \
		rmdir $${TESTDIR}; \
	fi;

# A note about Protobuf:
# We have a lot of logic here to make sure that the protoc we have henerates headers with exactly the same
# version requirements as the headers we already have.
# If not, we regenerate them.
# Doesn't handle Protobuf 3.12.3 weirdness; just make clean if you change flavors of Protobuf 3.12.3.
	
	
	
# run .pre-build before we make anything at all.
-include .pre-build

# for rebuilding just vg
clean-vg:
	$(RM) -f $(BIN_DIR)/$(EXE)
	$(RM) -f $(UNITTEST_SUPPORT_OBJ_DIR)/*.o $(UNITTEST_SUPPORT_OBJ_DIR)/*.d
	$(RM) -f $(UNITTEST_OBJ_DIR)/*.o $(UNITTEST_OBJ_DIR)/*.d
	$(RM) -f $(SUBCOMMAND_OBJ_DIR)/*.o $(SUBCOMMAND_OBJ_DIR)/*.d
	$(RM) -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d
	$(RM) -f $(SHARED_OBJ_DIR)/*.o $(SHARED_OBJ_DIR)/*.d
	$(RM) -f $(ALGORITHMS_OBJ_DIR)/*.o $(ALGORITHMS_OBJ_DIR)/*.d
	$(RM) -f $(ALGORITHMS_SHARED_OBJ_DIR)/*.o $(ALGORITHMS_SHARED_OBJ_DIR)/*.d
	$(RM) -f $(IO_OBJ_DIR)/*.o $(IO_OBJ_DIR)/*.d
	$(RM) -f $(IO_SHARED_OBJ_DIR)/*.o $(IO_SHARED_OBJ_DIR)/*.d
	$(RM) -f $(INC_DIR)/vg_git_version.hpp $(INC_DIR)/vg_system_version.hpp

clean: clean-vcflib
	$(RM) -r $(UNITTEST_BIN_DIR)
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(UNITTEST_SUPPORT_OBJ_DIR)
	$(RM) -r $(UNITTEST_OBJ_DIR)
	$(RM) -r $(SUBCOMMAND_OBJ_DIR)
	$(RM) -r $(IO_SHARED_OBJ_DIR)
	$(RM) -r $(IO_OBJ_DIR)
	$(RM) -r $(ALGORITHMS_SHARED_OBJ_DIR)
	$(RM) -r $(ALGORITHMS_OBJ_DIR)
	$(RM) -r $(CONFIG_OBJ_DIR)
	$(RM) -r $(SHARED_OBJ_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r share/
	cd $(DEP_DIR) && cd htslib && $(MAKE) clean
	cd $(DEP_DIR) && cd tabixpp && rm -f tabix.o libtabixpp.a
	cd $(DEP_DIR) && cd sonLib && $(MAKE) clean
	cd $(DEP_DIR) && cd sparsehash && $(MAKE) clean || true
	cd $(DEP_DIR) && cd fastahack && $(MAKE) clean
	cd $(DEP_DIR) && cd gcsa2 && $(MAKE) clean
	cd $(DEP_DIR) && cd gbwt && $(MAKE) clean
	cd $(DEP_DIR) && cd gbwtgraph && $(MAKE) clean
	cd $(DEP_DIR) && cd kff-cpp-api && rm -Rf build
	cd $(DEP_DIR) && cd gssw && $(MAKE) clean
	cd $(DEP_DIR) && cd ssw && cd src && $(MAKE) clean
	cd $(DEP_DIR) && cd progress_bar && $(MAKE) clean
	cd $(DEP_DIR) && cd sdsl-lite && ./uninstall.sh || true
	cd $(DEP_DIR) && cd libVCFH && $(MAKE) clean
	cd $(DEP_DIR) && cd vcflib && $(MAKE) clean
	cd $(DEP_DIR) && cd sha1 && $(MAKE) clean
	cd $(DEP_DIR) && cd structures && $(MAKE) clean
	cd $(DEP_DIR) && cd jemalloc && $(MAKE) clean || true
	cd $(DEP_DIR) && cd sublinear-Li-Stephens && $(MAKE) clean
	cd $(DEP_DIR) && cd libhandlegraph && rm -Rf build CMakeCache.txt CMakeFiles
	cd $(DEP_DIR) && cd libvgio && rm -Rf build CMakeCache.txt CMakeFiles
	cd $(DEP_DIR) && cd raptor && cd build && find . -not \( -name '.gitignore' -or -name 'pkg.m4' \) -delete
	# lru_cache is never built because it is header-only
	# bash-tap is never built either

clean-vcflib:
	cd $(DEP_DIR) && cd vcflib && $(MAKE) clean
	rm -f $(LIB_DIR)/libvcfh.a
	cd $(INC_DIR) && rm -f BedReader.h convert.h join.h mt19937ar.h split.h Variant.h vec128int.h veclib_types.h

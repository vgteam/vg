DEP_DIR:=./deps
SRC_DIR:=src
ALGORITHMS_SRC_DIR:=$(SRC_DIR)/algorithms
UNITTEST_SRC_DIR:=$(SRC_DIR)/unittest
SUBCOMMAND_SRC_DIR:=$(SRC_DIR)/subcommand
BIN_DIR:=bin
OBJ_DIR:=obj
ALGORITHMS_OBJ_DIR:=$(OBJ_DIR)/algorithms
UNITTEST_OBJ_DIR:=$(OBJ_DIR)/unittest
SUBCOMMAND_OBJ_DIR:=$(OBJ_DIR)/subcommand
LIB_DIR:=lib
# INC_DIR must be a relative path
INC_DIR:=include
CPP_DIR:=cpp
CWD:=$(shell pwd)

EXE:=vg

all: $(BIN_DIR)/$(EXE)

# Magic dependencies (see <http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr>)
include $(wildcard $(OBJ_DIR)/*.d)
include $(wildcard $(ALGORITHMS_OBJ_DIR)/*.d)
include $(wildcard $(UNITTEST_OBJ_DIR)/*.d)
include $(wildcard $(SUBCOMMAND_OBJ_DIR)/*.d)

CXXFLAGS := -O3 -fopenmp -Werror=return-type -std=c++11 -ggdb -g -MMD -MP $(CXXFLAGS)

LD_INCLUDE_FLAGS:=-I$(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR) -I$(CWD)/$(UNITTEST_SRC_DIR) -I$(CWD)/$(SUBCOMMAND_SRC_DIR) -I$(CWD)/$(CPP_DIR) -I$(CWD)/$(INC_DIR)/dynamic -I$(CWD)/$(INC_DIR)/sonLib $(shell pkg-config --cflags cairo)

LD_LIB_FLAGS:= -L$(CWD)/$(LIB_DIR) -lvcflib -lgssw -lssw -lprotobuf -lsublinearLS -lhts -lpthread -ljansson -lncurses -lgcsa2 -lgbwt -ldivsufsort -ldivsufsort64 -lvcfh -lgfakluge -lraptor2 -lsdsl -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -llz4 -lstructures -lvw -lboost_program_options -lallreduce
# Use pkg-config to find Cairo and all the libs it uses
LD_LIB_FLAGS += $(shell pkg-config --libs --static cairo)


ifeq ($(shell uname -s),Darwin)
	# We may need libraries from Macports
    # TODO: where does Homebrew keep libraries?
    ifeq ($(shell if [ -d /opt/local/lib ];then echo 1;else echo 0;fi), 1)
        # Use /opt/local/lib if present
	    LD_LIB_FLAGS += -L/opt/local/lib
    endif

    ifeq ($(shell if [ -d /usr/local/lib ];then echo 1;else echo 0;fi), 1)
        # Use /usr/local/lib if present.
        LD_LIB_FLAGS += -L/usr/local/lib
    endif

else
    # We are not running on OS X
	# We can also have a normal Unix rpath
	LD_LIB_FLAGS += -Wl,-rpath,$(CWD)/$(LIB_DIR)
	# Make sure to allow backtrace access to all our symbols, even those which are not exported.
	# Absolutely no help in a static build.
	LD_LIB_FLAGS += -rdynamic

	# We want to link against the elfutils libraries
	LD_LIB_FLAGS += -ldwfl -ldw -ldwelf -lelf -lebl
endif

# These libs need to come after libdw if used, because libdw depends on them
LD_LIB_FLAGS += -ldl -llzma

# Sometimes we need to filter the assembler output. The assembler can run during
# ./configure scripts, compiler calls, or $(MAKE) calls (other than $(MAKE)
# install). So we just stick $(FILTER) at the end of all such commands.
ifeq ($(shell uname -s),Darwin)
    # We need to apply a filter to all our build command output. This discards
    # all the assembler warnings which can overwhelm Travis log storage.
    FILTER=2>&1 | python $(CWD)/scripts/filter-noisy-assembler-warnings.py
    # For the filter to work and not just swallow errors we also need to turn on
    # pipefail in the shell
    SHELL=/bin/bash -o pipefail
else
    # No filter
    FILTER=
endif

ROCKSDB_PORTABLE=PORTABLE=1 # needed to build rocksdb without weird assembler options
# TODO: configure RPATH-equivalent on OS X for finding libraries without environment variables at runtime

# RocksDB's dependecies depend on whether certain compression libraries
# happen to be installed on the build system. Define a lazy macro to
# detect these from its self-configuration. It has to be lazy because
# the configuration (make_config.mk) won't exist until after RocksDB
# is built by this Makefile.
LD_LIB_FLAGS += -lrocksdb
ROCKSDB_LDFLAGS = $(shell grep PLATFORM_LDFLAGS deps/rocksdb/make_config.mk | cut -d '=' -f2 | sed s/-ljemalloc// | sed s/-ltcmalloc// | sed s/-ltbb//)

# When building statically, we need to tell the linker not to bail if it sees multiple definitions.
# libc on e.g. our Jenkins host does not define malloc as weak, so tcmalloc can't override it in a static build.
# TODO: Why did this problem only begin to happen when libvw was added?
STATIC_FLAGS=-static -static-libstdc++ -static-libgcc -Wl,--allow-multiple-definition 

# These are put into libvg. Grab everything except main.
OBJ = $(filter-out $(OBJ_DIR)/main.o,$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp)))
# And all the algorithms
ALGORITHMS_OBJ = $(patsubst $(ALGORITHMS_SRC_DIR)/%.cpp,$(ALGORITHMS_OBJ_DIR)/%.o,$(wildcard $(ALGORITHMS_SRC_DIR)/*.cpp))

# These aren't put into libvg. But they do go into the main vg binary to power its self-test.
UNITTEST_OBJ = $(patsubst $(UNITTEST_SRC_DIR)/%.cpp,$(UNITTEST_OBJ_DIR)/%.o,$(wildcard $(UNITTEST_SRC_DIR)/*.cpp))

# These aren't put into libvg, but they provide subcommand implementations for the vg bianry
SUBCOMMAND_OBJ = $(patsubst $(SUBCOMMAND_SRC_DIR)/%.cpp,$(SUBCOMMAND_OBJ_DIR)/%.o,$(wildcard $(SUBCOMMAND_SRC_DIR)/*.cpp))

RAPTOR_DIR:=deps/raptor
PROTOBUF_DIR:=deps/protobuf
GPERF_DIR:=deps/gperftools
SDSL_DIR:=deps/sdsl-lite
SNAPPY_DIR:=deps/snappy
ROCKSDB_DIR:=deps/rocksdb
GCSA2_DIR:=deps/gcsa2
GBWT_DIR:=deps/gbwt
PROGRESS_BAR_DIR:=deps/progress_bar
FASTAHACK_DIR:=deps/fastahack
FERMI_DIR:=deps/fermi-lite
HTSLIB_DIR:=deps/htslib
VCFLIB_DIR:=deps/vcflib
GSSW_DIR:=deps/gssw
SPARSEHASH_DIR:=deps/sparsehash
SPARSEPP_DIR:=deps/sparsepp
SHA1_DIR:=deps/sha1
DYNAMIC_DIR:=deps/DYNAMIC
SSW_DIR:=deps/ssw/src
LINLS_DIR:=deps/sublinear-Li-Stephens
STRUCTURES_DIR:=deps/structures
BACKWARD_CPP_DIR:=deps/backward-cpp
ELFUTILS_DIR:=deps/elfutils
BOOST_DIR:=deps/boost-subset
VOWPALWABBIT_DIR:=deps/vowpal_wabbit

# Dependencies that go into libvg's archive
# These go in libvg but come from dependencies
DEP_OBJ =
DEP_OBJ += $(OBJ_DIR)/vg.pb.o 
DEP_OBJ += $(OBJ_DIR)/progress_bar.o
DEP_OBJ += $(OBJ_DIR)/sha1.o
DEP_OBJ += $(OBJ_DIR)/Fasta.o


# These are libraries that we need to build before we link vg.
# It would be nice to dump their contents into libvg to make it stand-alone.
# But that requires fancy ar scripting.
# If you just pass them to ar it puts the library *file* in libvg where nothing can read it.
LIB_DEPS =
LIB_DEPS += $(LIB_DIR)/libprotobuf.a
LIB_DEPS += $(LIB_DIR)/libsdsl.a
LIB_DEPS += $(LIB_DIR)/libssw.a
LIB_DEPS += $(LIB_DIR)/libsnappy.a
LIB_DEPS += $(LIB_DIR)/librocksdb.a
LIB_DEPS += $(LIB_DIR)/libgcsa2.a
LIB_DEPS += $(LIB_DIR)/libgbwt.a
LIB_DEPS += $(LIB_DIR)/libhts.a
LIB_DEPS += $(LIB_DIR)/libvcflib.a
LIB_DEPS += $(LIB_DIR)/libgssw.a
LIB_DEPS += $(LIB_DIR)/libvcfh.a
LIB_DEPS += $(LIB_DIR)/libgfakluge.a
LIB_DEPS += $(LIB_DIR)/libsonlib.a
LIB_DEPS += $(LIB_DIR)/libpinchesandcacti.a
LIB_DEPS += $(LIB_DIR)/libraptor2.a
LIB_DEPS += $(LIB_DIR)/libfml.a
LIB_DEPS += $(LIB_DIR)/libsublinearLS.a
LIB_DEPS += $(LIB_DIR)/libstructures.a
LIB_DEPS += $(LIB_DIR)/libvw.a
LIB_DEPS += $(LIB_DIR)/liballreduce.a
LIB_DEPS += $(LIB_DIR)/libboost_program_options.a
ifneq ($(shell uname -s),Darwin)
	# On non-Mac (i.e. Linux), where ELF binaries are used, pull in libdw which
	# backward-cpp will use.
	LIB_DEPS += $(LIB_DIR)/libdw.a
	LIB_DEPS += $(LIB_DIR)/libdwfl.a
	LIB_DEPS += $(LIB_DIR)/libdwelf.a
	LIB_DEPS += $(LIB_DIR)/libebl.a
	LIB_DEPS += $(LIB_DIR)/libelf.a
endif

# common dependencies to build before all vg src files
DEPS = $(LIB_DEPS)
DEPS += $(CPP_DIR)/vg.pb.h
DEPS += $(INC_DIR)/gcsa/gcsa.h
DEPS += $(INC_DIR)/gbwt/dynamic_gbwt.h
DEPS += $(INC_DIR)/lru_cache.h
DEPS += $(INC_DIR)/dynamic.hpp
DEPS += $(INC_DIR)/sparsehash/sparse_hash_map
DEPS += $(INC_DIR)/sparsepp/spp.h
DEPS += $(INC_DIR)/gfakluge.hpp
DEPS += $(INC_DIR)/sha1.hpp
DEPS += $(INC_DIR)/progress_bar.hpp
DEPS += $(INC_DIR)/backward.hpp

ifneq ($(shell uname -s),Darwin)
	DEPS += $(LIB_DIR)/libtcmalloc_minimal.a
	LD_LIB_FLAGS += -ltcmalloc_minimal
endif

.PHONY: clean get-deps deps test set-path static docs .pre-build .check-environment

$(BIN_DIR)/vg: $(OBJ_DIR)/main.o $(LIB_DIR)/libvg.a $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o $(BIN_DIR)/vg $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) -lvg $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

static: $(OBJ_DIR)/main.o $(LIB_DIR)/libvg.a $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/vg $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) -lvg $(STATIC_FLAGS) $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

$(LIB_DIR)/libvg.a: $(OBJ) $(ALGORITHMS_OBJ) $(DEP_OBJ) $(DEPS)
	rm -f $@
	ar rs $@ $(OBJ) $(ALGORITHMS_OBJ) $(DEP_OBJ)

# We have system-level deps to install
get-deps:
	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake libtool jq samtools curl unzip redland-utils librdf-dev cmake pkg-config wget bc gtk-doc-tools raptor2-utils rasqal-utils bison flex gawk libgoogle-perftools-dev liblz4-dev liblzma-dev libcairo2-dev libpixman-1-dev libffi-dev

# And we have submodule deps to build
deps: $(DEPS)

test: $(BIN_DIR)/vg $(LIB_DIR)/libvg.a test/build_graph $(BIN_DIR)/shuf
	. ./source_me.sh && cd test && prove -v t

docs: $(SRC_DIR)/*.cpp $(SRC_DIR)/*.hpp $(SUBCOMMAND_SRC_DIR)/*.cpp $(SUBCOMMAND_SRC_DIR)/*.hpp $(UNITTEST_SRC_DIR)/*.cpp $(UNITTEST_SRC_DIR)/*.hpp $(CPP_DIR)/vg.pb.cc
	doxygen
	cd doc && sphinx-build -b html . sphinx

# Hack to use gshuf or shuf as appropriate to the platform when testing
$(BIN_DIR)/shuf:
ifeq ($(shell uname -s),Darwin)
	ln -s `which gshuf` $(BIN_DIR)/shuf
else
	ln -s `which shuf` $(BIN_DIR)/shuf
endif

# Make sure we have protoc built, and the protobuf lib, both of which come from the same command using this fake intermediate
bin/protoc: .rebuild-protobuf
$(LIB_DIR)/libprotobuf.a: .rebuild-protobuf
	# intermediate targets don't trigger a rebuild just because they're missing.
.INTERMEDIATE: .rebuild-protobuf
	# Make sure to delete outdated libs and headers before rebuilding
	# Outdated headers can get picked up during the build
.rebuild-protobuf: deps/protobuf/src/google/protobuf/*.cc
	rm -rf $(LIB_DIR)/libprotobuf* $(LIB_DIR)/libprotoc*
	rm -Rf include/google/protobuf/
	+. ./source_me.sh && cd $(PROTOBUF_DIR) && ./autogen.sh && export DIST_LANG=cpp && ./configure --prefix="$(CWD)" $(FILTER) && $(MAKE) $(FILTER) && $(MAKE) install && export PATH=$(CWD)/bin:$$PATH

test/build_graph: test/build_graph.cpp $(LIB_DIR)/libvg.a $(CPP_DIR)/vg.pb.h $(SRC_DIR)/json2pb.h $(SRC_DIR)/vg.hpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o test/build_graph test/build_graph.cpp $(LD_INCLUDE_FLAGS) -lvg $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS) $(FILTER)

# remove annoying large alloc messages from tcmalloc
$(GPERF_DIR)/src/tcmalloc.cc.bak:
	cp $(GPERF_DIR)/src/tcmalloc.cc $(GPERF_DIR)/src/tcmalloc.cc.bak
	sed 's/printer.printf("tcmalloc: large alloc/return; printer.printf("tcmalloc: large alloc/' $(GPERF_DIR)/src/tcmalloc.cc.bak >$(GPERF_DIR)/src/tcmalloc.cc

$(LIB_DIR)/libtcmalloc_minimal.a: $(GPERF_DIR)/src/tcmalloc.cc.bak
	+. ./source_me.sh && cd $(GPERF_DIR) && ./autogen.sh && ./configure --prefix=`pwd` $(FILTER) && $(MAKE) $(FILTER) && $(MAKE) install && cp -r lib/* $(CWD)/$(LIB_DIR)/ && cp -r bin/* $(CWD)/$(BIN_DIR)/ && cp -r include/* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libsdsl.a: $(SDSL_DIR)/lib/*.cpp $(SDSL_DIR)/include/sdsl/*.hpp
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cd $(SDSL_DIR) && AS_INTEGRATED_ASSEMBLER=1 BUILD_PORTABLE=1 ./install.sh $(CWD) $(FILTER)
else
	+. ./source_me.sh && cd $(SDSL_DIR) && BUILD_PORTABLE=1 ./install.sh $(CWD) $(FILTER)
endif	


$(LIB_DIR)/libssw.a: $(SSW_DIR)/*.c $(SSW_DIR)/*.h
	+. ./source_me.sh && cd $(SSW_DIR) && $(MAKE) $(FILTER) && ar rs $(CWD)/$(LIB_DIR)/libssw.a ssw.o ssw_cpp.o && cp ssw_cpp.h ssw.h $(CWD)/$(LIB_DIR)

$(LIB_DIR)/libsnappy.a: $(SNAPPY_DIR)/*.cc $(SNAPPY_DIR)/*.h
	+. ./source_me.sh && cd $(SNAPPY_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) $(FILTER) && $(MAKE) $(FILTER) && $(MAKE) install

$(LIB_DIR)/librocksdb.a: $(LIB_DIR)/libtcmalloc_minimal.a $(LIB_DIR)/libsnappy.a $(ROCKSDB_DIR)/db/*.cc $(ROCKSDB_DIR)/db/*.h
	+. ./source_me.sh && cd $(ROCKSDB_DIR) && $(ROCKSDB_PORTABLE) DISABLE_JEMALLOC=1 $(MAKE) static_lib $(FILTER) && mv librocksdb.a $(CWD)/${LIB_DIR}/ && cp -r include/* $(CWD)/$(INC_DIR)/

$(INC_DIR)/gcsa/gcsa.h: $(LIB_DIR)/libgcsa2.a

$(LIB_DIR)/libgcsa2.a: $(LIB_DIR)/libsdsl.a $(wildcard $(GCSA2_DIR)/*.cpp) $(wildcard $(GCSA2_DIR)/include/gcsa/*.h)
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cd $(GCSA2_DIR) && AS_INTEGRATED_ASSEMBLER=1 $(MAKE) libgcsa2.a $(FILTER) && mv libgcsa2.a $(CWD)/$(LIB_DIR) && cp -r include/gcsa $(CWD)/$(INC_DIR)/
else
	+. ./source_me.sh && cd $(GCSA2_DIR) && $(MAKE) libgcsa2.a $(FILTER) && mv libgcsa2.a $(CWD)/$(LIB_DIR) && cp -r include/gcsa $(CWD)/$(INC_DIR)/
endif

$(INC_DIR)/gbwt/dynamic_gbwt.h: $(LIB_DIR)/libgbwt.a
$(LIB_DIR)/libgbwt.a: $(LIB_DIR)/libsdsl.a $(wildcard $(GBWT_DIR)/*.cpp) $(wildcard $(GBWT_DIR)/include/gbwt/*.h)
ifeq ($(shell uname -s),Darwin)
	+. ./source_me.sh && cd $(GBWT_DIR) && AS_INTEGRATED_ASSEMBLER=1 $(MAKE) libgbwt.a $(FILTER) && mv libgbwt.a $(CWD)/$(LIB_DIR) && cp -r include/gbwt $(CWD)/$(INC_DIR)/
else
	+. ./source_me.sh && cd $(GBWT_DIR) && $(MAKE) libgbwt.a $(FILTER) && mv libgbwt.a $(CWD)/$(LIB_DIR) && cp -r include/gbwt $(CWD)/$(INC_DIR)/
endif

$(INC_DIR)/progress_bar.hpp: $(PROGRESS_BAR_DIR)/progress_bar.hpp
	+cp $(PROGRESS_BAR_DIR)/progress_bar.hpp $(CWD)/$(INC_DIR)

$(OBJ_DIR)/progress_bar.o: $(PROGRESS_BAR_DIR)/*.hpp $(PROGRESS_BAR_DIR)/*.cpp
	+cd $(PROGRESS_BAR_DIR) && $(MAKE) $(FILTER) && cp progress_bar.o $(CWD)/$(OBJ_DIR)

$(OBJ_DIR)/Fasta.o: $(FASTAHACK_DIR)/*.h $(FASTAHACK_DIR)/*.cpp
	+cd $(FASTAHACK_DIR) && $(MAKE) $(FILTER) && mv Fasta.o $(CWD)/$(OBJ_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

$(LIB_DIR)/libhts.a: $(HTSLIB_DIR)/*.c $(HTSLIB_DIR)/*.h
	+cd $(HTSLIB_DIR) && $(MAKE) lib-static $(FILTER) && mv libhts.a $(CWD)/$(LIB_DIR) && cp *.h $(CWD)/$(INC_DIR) && cp -r htslib $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libvcflib.a: $(VCFLIB_DIR)/src/*.cpp $(VCFLIB_DIR)/src/*.hpp $(VCFLIB_DIR)/intervaltree/*.cpp $(VCFLIB_DIR)/intervaltree/*.h $(VCFLIB_DIR)/tabixpp/*.cpp $(VCFLIB_DIR)/tabixpp/*.hpp $(VCFLIB_DIR)/tabixpp/htslib/*.c $(VCFLIB_DIR)/tabixpp/htslib/htslib/*.h $(VCFLIB_DIR)/tabixpp/htslib/cram/*.c $(VCFLIB_DIR)/tabixpp/htslib/cram/*.h
	+. ./source_me.sh && cd $(VCFLIB_DIR) && $(MAKE) libvcflib.a $(FILTER) && cp lib/* $(CWD)/$(LIB_DIR)/ && cp include/* $(CWD)/$(INC_DIR)/ && cp intervaltree/*.h $(CWD)/$(INC_DIR)/ && cp src/*.h* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libgssw.a: $(GSSW_DIR)/src/gssw.c $(GSSW_DIR)/src/gssw.h
	+cd $(GSSW_DIR) && $(MAKE) $(FILTER) && cp lib/* $(CWD)/$(LIB_DIR)/ && cp obj/* $(CWD)/$(OBJ_DIR) && cp src/*.h $(CWD)/$(INC_DIR)

$(INC_DIR)/lru_cache.h: $(DEP_DIR)/lru_cache/*.h $(DEP_DIR)/lru_cache/*.cc
	+cd $(DEP_DIR)/lru_cache && $(MAKE) $(FILTER) && cp *.h* $(CWD)/$(INC_DIR)/

$(INC_DIR)/dynamic.hpp: $(DYNAMIC_DIR)/include/*.hpp $(DYNAMIC_DIR)/include/internal/*.hpp 
	+cat $(DYNAMIC_DIR)/include/dynamic.hpp | sed 's%<internal/%<dynamic/%' >$(INC_DIR)/dynamic.hpp && cp -r $(CWD)/$(DYNAMIC_DIR)/include/internal $(CWD)/$(INC_DIR)/dynamic

$(INC_DIR)/sparsehash/sparse_hash_map: $(wildcard $(SPARSEHASH_DIR)/**/*.cc) $(wildcard $(SPARSEHASH_DIR)/**/*.h) 
	+cd $(SPARSEHASH_DIR) && ./autogen.sh && LDFLAGS="-L/opt/local/lib" ./configure --prefix=$(CWD) $(FILTER) && $(MAKE) $(FILTER) && $(MAKE) install

$(INC_DIR)/sparsepp/spp.h: $(wildcard $(SPARSEHASH_DIR)/sparsepp/*.h)
	+cp -r $(SPARSEPP_DIR)/sparsepp $(INC_DIR)/

#$(INC_DIR)/Variant.h
$(LIB_DIR)/libvcfh.a: $(DEP_DIR)/libVCFH/*.cpp $(DEP_DIR)/libVCFH/*.hpp 
	+cd $(DEP_DIR)/libVCFH && $(MAKE) $(FILTER) && cp libvcfh.a $(CWD)/$(LIB_DIR)/ && cp vcfheader.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libgfakluge.a: $(INC_DIR)/gfakluge.hpp $(DEP_DIR)/gfakluge/src/*.hpp $(DEP_DIR)/gfakluge/src/*.cpp
	+cd $(DEP_DIR)/gfakluge && $(MAKE) libgfakluge.a $(FILTER) && cp libgfakluge.a $(CWD)/$(LIB_DIR)/

$(INC_DIR)/gfakluge.hpp: $(DEP_DIR)/gfakluge/src/gfakluge.hpp
	+cp $(DEP_DIR)/gfakluge/src/*.hpp $(CWD)/$(INC_DIR)/ && cp $(DEP_DIR)/gfakluge/src/tinyFA/*.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libsonlib.a: $(CWD)/$(DEP_DIR)/sonLib/C/inc/*.h $(CWD)/$(DEP_DIR)/sonLib/C/impl/*.c
	+cd $(DEP_DIR)/sonLib && kyotoTycoonLib="" $(MAKE) $(FILTER) && cp lib/sonLib.a $(CWD)/$(LIB_DIR)/libsonlib.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

$(LIB_DIR)/libpinchesandcacti.a: $(LIB_DIR)/libsonlib.a $(CWD)/$(DEP_DIR)/pinchesAndCacti/inc/*.h $(CWD)/$(DEP_DIR)/pinchesAndCacti/impl/*.c
	+cd $(DEP_DIR)/pinchesAndCacti && $(MAKE) $(FILTER) && cd $(CWD)/$(DEP_DIR)/sonLib && cp lib/stPinchesAndCacti.a $(CWD)/$(LIB_DIR)/libpinchesandcacti.a && cp lib/3EdgeConnected.a $(CWD)/$(LIB_DIR)/lib3edgeconnected.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

# When building raptor we need to make sure to pre-generate and fix up the lexer
$(LIB_DIR)/libraptor2.a: $(RAPTOR_DIR)/src/*.c $(RAPTOR_DIR)/src/*.h
	+cd $(RAPTOR_DIR)/build && cmake .. && rm -f src/turtle_parser.c && rm -f src/turtle_lexer.c && make turtle_lexer_tgt && make -f src/CMakeFiles/raptor2.dir/build.make src/turtle_lexer.c && sed -i.bak '/yycleanup/d' src/turtle_lexer.c && $(MAKE) $(FILTER) && cp src/libraptor2.a $(CWD)/$(LIB_DIR) && mkdir -p $(CWD)/$(INC_DIR)/raptor2 && cp src/*.h $(CWD)/$(INC_DIR)/raptor2

$(LIB_DIR)/libstructures.a: $(STRUCTURES_DIR)/src/include/structures/*.hpp $(STRUCTURES_DIR)/src/*.cpp 
	+. ./source_me.sh && cd $(STRUCTURES_DIR) && $(MAKE) lib/libstructures.a $(FILTER) && cp lib/libstructures.a $(CWD)/$(LIB_DIR)/ && cp -r src/include/structures $(CWD)/$(INC_DIR)/
    
# To build libvw we need to point it at our Boost, but then configure decides
# it needs to build vwdll, which depends on codecvt, which isn't actually
# shipped in the GCC 4.9 STL. So we hack vwdll AKA libvw_c_wrapper out of the
# build.
# Also, autogen.sh looks for Boost in the system, and who knows what it will do
# if it doesn't find it, so let it fail.
$(LIB_DIR)/libvw.a: $(LIB_DIR)/libboost_program_options.a $(VOWPALWABBIT_DIR)/* $(VOWPALWABBIT_DIR)/vowpalwabbit/*
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && sed -i -e 's/libvw_c_wrapper\.pc//g' Makefile.am
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && sed -i -e 's/libvw_c_wrapper\.la//g' vowpalwabbit/Makefile.am
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && sed -i -e '/libvw_c_wrapper\.pc/d' configure.ac
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && sed -i -e '/vwdll/d' Makefile.am
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && sed -i -e '/libvw_c_wrapper/d' vowpalwabbit/Makefile.am
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && (./autogen.sh || true)
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && ./configure --with-boost=$(CWD)
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && $(MAKE) $(FILTER)
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && cp vowpalwabbit/.libs/libvw.a vowpalwabbit/.libs/liballreduce.a $(CWD)/$(LIB_DIR)/
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && mkdir -p $(CWD)/$(INC_DIR)/vowpalwabbit
	+. ./source_me.sh && cd $(VOWPALWABBIT_DIR) && cp vowpalwabbit/*.h $(CWD)/$(INC_DIR)/vowpalwabbit/

$(LIB_DIR)/liballreduce.a: $(LIB_DIR)/libvw.a

$(LIB_DIR)/libboost_program_options.a: $(BOOST_DIR)/libs/program_options/src/* $(BOOST_DIR)/boost/program_options/*
	+. ./source_me.sh && cd $(BOOST_DIR) && ./bootstrap.sh --with-libraries=program_options --libdir=$(CWD)/$(LIB_DIR) --includedir=$(CWD)/$(INC_DIR) $(FILTER) && ./b2 --ignore-site-config --link=static install $(FILTER)
	+. ./source_me.sh && if [[ $(shell uname -s) == "Darwin" ]]; then install_name_tool -id $(CWD)/$(LIB_DIR)/libboost_program_options.dylib $(CWD)/$(LIB_DIR)/libboost_program_options.dylib; fi

$(INC_DIR)/sha1.hpp: $(SHA1_DIR)/sha1.hpp
	+cp $(SHA1_DIR)/*.h* $(CWD)/$(INC_DIR)/

$(INC_DIR)/backward.hpp: $(BACKWARD_CPP_DIR)/backward.hpp
	+cp $(BACKWARD_CPP_DIR)/backward.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libebl.a: $(LIB_DIR)/libelf.a

$(LIB_DIR)/libdw.a: $(LIB_DIR)/libelf.a

$(LIB_DIR)/libdwelf.a: $(LIB_DIR)/libelf.a

$(LIB_DIR)/libdwfl.a: $(LIB_DIR)/libelf.a

# We can't build elfutils from Git without "maintainer mode".
# There are some release-only headers or something that it complains it can't find otherwise.
# We also don't do a normal make and make install here because we don't want to build and install all the elfutils binaries and libasm.
$(LIB_DIR)/libelf.a: $(ELFUTILS_DIR)/libebl/* $(ELFUTILS_DIR)/libdw/* $(ELFUTILS_DIR)/libelf/* $(ELFUTILS_DIR)/src/*
	+cd $(ELFUTILS_DIR) && autoreconf -i -f && ./configure --enable-maintainer-mode --prefix=$(CWD) $(FILTER)
	+cd $(ELFUTILS_DIR)/libelf && $(MAKE) libelf.a $(FILTER)
	+cd $(ELFUTILS_DIR)/libebl && $(MAKE) libebl.a $(FILTER)
	+cd $(ELFUTILS_DIR)/libdw && $(MAKE) libdw.a known-dwarf.h $(FILTER)
	+cd $(ELFUTILS_DIR)/libdwfl && $(MAKE) libdwfl.a $(FILTER)
	+cd $(ELFUTILS_DIR)/libdwelf && $(MAKE) libdwelf.a $(FILTER)
	+cd $(ELFUTILS_DIR) && mkdir -p $(CWD)/$(INC_DIR)/elfutils && cp libdw/known-dwarf.h libdw/libdw.h libebl/libebl.h libelf/elf-knowledge.h version.h libdwfl/libdwfl.h libdwelf/libdwelf.h $(CWD)/$(INC_DIR)/elfutils && cp libelf/gelf.h libelf/libelf.h libdw/dwarf.h $(CWD)/$(INC_DIR) && cp libebl/libebl.a libdw/libdw.a libdwfl/libdwfl.a libdwelf/libdwelf.a libelf/libelf.a $(CWD)/$(LIB_DIR)/

$(OBJ_DIR)/sha1.o: $(SHA1_DIR)/sha1.cpp $(SHA1_DIR)/sha1.hpp
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(FILTER)

$(LIB_DIR)/libfml.a: $(FERMI_DIR)/*.h $(FERMI_DIR)/*.c
	cd $(FERMI_DIR) && $(MAKE) $(FILTER) && cp *.h $(CWD)/$(INC_DIR)/ && cp libfml.a $(CWD)/$(LIB_DIR)/

$(LIB_DIR)/libsublinearLS.a: $(LINLS_DIR)/src/*.cpp $(LINLS_DIR)/src/*.hpp $(LIB_DIR)/libhts.a
	cd $(LINLS_DIR) && INCLUDE_FLAGS="-I$(CWD)/$(INC_DIR)" $(MAKE) libs $(FILTER) && cp lib/libsublinearLS.a $(CWD)/$(LIB_DIR)/ && mkdir -p $(CWD)/$(INC_DIR)/sublinearLS && cp src/*.hpp $(CWD)/$(INC_DIR)/sublinearLS/

# Auto-git-versioning
$(INC_DIR)/vg_git_version.hpp: .git
	@echo "#define VG_GIT_VERSION \"$(shell git describe --always --tags || echo unknown)\"" > $@

# Build an environment version file with this phony target.
# If it's not the same as the old one, replace the old one.
# If it is the same, do nothing and don't rebuild dependent targets.
.check-environment:
	@echo "#define VG_COMPILER_VERSION \"$(shell $(CXX) --version | head -n 1)\"" > $(INC_DIR)/vg_environment_version.hpp.tmp
	@echo "#define VG_OS \"$(shell uname)\"" >> $(INC_DIR)/vg_environment_version.hpp.tmp
	@echo "#define VG_BUILD_USER \"$(shell whoami)\"" >> $(INC_DIR)/vg_environment_version.hpp.tmp
	@echo "#define VG_BUILD_HOST \"$(shell hostname)\"" >> $(INC_DIR)/vg_environment_version.hpp.tmp
	@diff $(INC_DIR)/vg_environment_version.hpp.tmp $(INC_DIR)/vg_environment_version.hpp >/dev/null || cp $(INC_DIR)/vg_environment_version.hpp.tmp $(INC_DIR)/vg_environment_version.hpp
	@rm -f $(INC_DIR)/vg_environment_version.hpp.tmp

# The way to get the actual file is to maybe replace it.
$(INC_DIR)/vg_environment_version.hpp: .check-environment
	

# Not important if .git isn't real
.git:

###################################
## VG source code compilation begins here
####################################

include/stream.hpp: src/stream.hpp
	cp src/stream.hpp include/stream.hpp

$(OBJ_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.o
	cp $(CPP_DIR)/vg.pb.o $(OBJ_DIR)/vg.pb.o

$(CPP_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.cc

$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h 

$(CPP_DIR)/vg.pb.h: $(LIB_DIR)/libprotobuf.a bin/protoc $(SRC_DIR)/vg.proto
	+. ./source_me.sh && ./bin/protoc $(SRC_DIR)/vg.proto --proto_path=$(SRC_DIR) --cpp_out=cpp
	+cp $@ $(INC_DIR)

$(OBJ_DIR)/version.o: $(SRC_DIR)/version.cpp $(SRC_DIR)/version.hpp $(INC_DIR)/vg_git_version.hpp $(INC_DIR)/vg_environment_version.hpp

########################
## Pattern Rules
########################

# Define a default rule for building objects from CPP files
# Depend on the .d file so we rebuild if dependency info is missing/deleted
# Use static pattern rules so the dependency files will not be ignored if the output exists
# See <https://stackoverflow.com/a/34983297>
$(OBJ) $(OBJ_DIR)/main.o: $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(FILTER)
$(ALGORITHMS_OBJ): $(ALGORITHMS_OBJ_DIR)/%.o : $(ALGORITHMS_SRC_DIR)/%.cpp $(ALGORITHMS_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(FILTER)
$(SUBCOMMAND_OBJ): $(SUBCOMMAND_OBJ_DIR)/%.o : $(SUBCOMMAND_SRC_DIR)/%.cpp $(SUBCOMMAND_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(FILTER)
$(UNITTEST_OBJ): $(UNITTEST_OBJ_DIR)/%.o : $(UNITTEST_SRC_DIR)/%.cpp $(UNITTEST_OBJ_DIR)/%.d $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(FILTER)
        
# Protobuf stuff builds into its same directory
$(CPP_DIR)/%.o : $(CPP_DIR)/%.cc $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(FILTER)

# Use a fake rule to build .d files, so we don't complain if they don't exist.
$(OBJ_DIR)/%.d: ;
$(ALGORITHMS_OBJ_DIR)/%.d: ;
$(SUBCOMMAND_OBJ_DIR)/%.d: ;
$(UNITTEST_OBJ_DIR)/%.d: ;

# Don't delete them.
.PRECIOUS: $(OBJ_DIR)/%.d $(ALGORITHMS_OBJ_DIR)/%.d $(SUBCOMMAND_OBJ_DIR)/%.d $(UNITTEST_OBJ_DIR)/%.d

# Use no implicit rules
.SUFFIXES:

###################################
## VG source code compilation ends here
####################################



.pre-build:
	@if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	@if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	@if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	@if [ ! -d $(ALGORITHMS_OBJ_DIR) ]; then mkdir -p $(ALGORITHMS_OBJ_DIR); fi
	@if [ ! -d $(UNITTEST_OBJ_DIR) ]; then mkdir -p $(UNITTEST_OBJ_DIR); fi
	@if [ ! -d $(SUBCOMMAND_OBJ_DIR) ]; then mkdir -p $(SUBCOMMAND_OBJ_DIR); fi
	@if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	@if [ ! -d $(CPP_DIR) ]; then mkdir -p $(CPP_DIR); fi

# run .pre-build before we make anything at all.
-include .pre-build

# for rebuilding just vg
clean-vg:
	$(RM) -r $(BIN_DIR)/vg
	$(RM) -r $(UNITTEST_OBJ_DIR)/*.o $(UNITTEST_OBJ_DIR)/*.d
	$(RM) -r $(SUBCOMMAND_OBJ_DIR)/*.o $(SUBCOMMAND_OBJ_DIR)/*.d
	$(RM) -r $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d
	$(RM) -r $(CPP_DIR)/*.o $(CPP_DIR)/*.d $(CPP_DIR)/*.cc $(CPP_DIR)/*.h
	$(RM) -f $(INC_DIR)/vg_git_version.hpp $(INC_DIR)/vg_system_version.hpp

clean: clean-rocksdb clean-protobuf clean-vcflib
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(UNITTEST_OBJ_DIR)
	$(RM) -r $(SUBCOMMAND_OBJ_DIR)
	$(RM) -r $(ALGORITHMS_OBJ_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r $(CPP_DIR)
	$(RM) -r share/
	cd $(DEP_DIR) && cd sparsehash && $(MAKE) clean
	cd $(DEP_DIR) && cd htslib && $(MAKE) clean
	cd $(DEP_DIR) && cd fastahack && $(MAKE) clean
	cd $(DEP_DIR) && cd gcsa2 && $(MAKE) clean
	cd $(DEP_DIR) && cd gbwt && $(MAKE) clean
	cd $(DEP_DIR) && cd gssw && $(MAKE) clean
	cd $(DEP_DIR) && cd progress_bar && $(MAKE) clean
	cd $(DEP_DIR) && cd sdsl-lite && ./uninstall.sh || true
	cd $(DEP_DIR) && cd libVCFH && $(MAKE) clean
	cd $(DEP_DIR) && cd gfakluge && $(MAKE) clean
	cd $(DEP_DIR) && cd sha1 && $(MAKE) clean
	cd $(DEP_DIR) && cd structures && $(MAKE) clean
	cd $(DEP_DIR) && cd gperftools && $(MAKE) clean
	cd $(DEP_DIR) && cd vowpal_wabbit && $(MAKE) clean
	rm -Rf $(RAPTOR_DIR)/build/*
	## TODO vg source code
	## TODO LRU_CACHE
	## TODO bash-tap

clean-rocksdb:
	cd $(DEP_DIR) && cd rocksdb && $(MAKE) clean
	rm -f $(LIB_DIR)/librocksdb.a 
	rm -rf $(INC_DIR)/rocksdb/

clean-protobuf:
	cd $(DEP_DIR) && cd protobuf && $(MAKE) clean
	rm -f $(LIB_DIR)/libprotobuf.a
	rm -rf $(INC_DIR)/google/protobuf/

clean-vcflib:
	cd $(DEP_DIR) && cd vcflib && $(MAKE) clean
	rm -f $(LIB_DIR)/libvcfh.a
	cd $(INC_DIR) && rm -f BedReader.h convert.h join.h mt19937ar.h split.h Variant.h vec128int.h veclib_types.h

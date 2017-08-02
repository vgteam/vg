DEP_DIR:=./deps
SRC_DIR:=src
ALGORITHMS_SRC_DIR:=$(SRC_DIR)/algorithms
UNITTEST_SRC_DIR:=$(SRC_DIR)/unittest
SUBCOMMAND_SRC_DIR:=$(SRC_DIR)/subcommand
BIN_DIR:=bin
OBJ_DIR:=obj
UNITTEST_OBJ_DIR:=$(OBJ_DIR)/unittest
SUBCOMMAND_OBJ_DIR:=$(OBJ_DIR)/subcommand
LIB_DIR:=lib
INC_DIR:=include
CPP_DIR:=cpp

EXE:=vg

CXXFLAGS:=-O3 -msse4.1 -fopenmp -std=c++11 -ggdb -g

CWD:=$(shell pwd)


LD_INCLUDE_FLAGS:=-I$(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR) -I$(CWD)/$(UNITTEST_SRC_DIR) -I$(CWD)/$(SUBCOMMAND_SRC_DIR) -I$(CWD)/$(CPP_DIR) -I$(CWD)/$(INC_DIR)/dynamic -I$(CWD)/$(INC_DIR)/sonLib -I$(CWD)/$(INC_DIR)/gcsa
LD_LIB_FLAGS:= -L$(CWD)/$(LIB_DIR) -lvcflib -lgssw -lssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lgcsa2 -lxg -ldivsufsort -ldivsufsort64 -lvcfh -lgfakluge -lraptor2 -lsupbub -lsdsl -lpinchesandcacti -l3edgeconnected -lsonlib -lfml

ifeq ($(shell uname -s),Darwin)
	# We may need libraries from Macports
	# TODO: where does Homebrew keep libraries?
	ifeq ($(shell if [ -d /opt/local/lib ];then echo 1;else echo 0;fi), 1)
	LD_LIB_FLAGS += -L/opt/local/lib
else
	CXXFLAGS += -march=native -mtune=native
endif
ifeq ($(shell if [ -d /usr/local/lib ];then echo 1;else echo 0;fi), 1)
	LD_LIB_FLAGS += -L/usr/local/lib
endif
else
	# We can also have a normal Unix rpath
	LD_LIB_FLAGS += -Wl,-rpath,$(CWD)/$(LIB_DIR)
endif

ROCKSDB_PORTABLE=PORTABLE=1 # needed to build rocksdb without weird assembler options
# TODO: configure RPATH-equivalent on OS X for finding libraries without environment variables at runtime

# RocksDB's dependecies depend on whether certain compression libraries
# happen to be installed on the build system. Define a lazy macro to
# detect these from its self-configuration. It has to be lazy because
# the configuration (make_config.mk) won't exist until after RocksDB
# is built by this Makefile.
LD_LIB_FLAGS += -lrocksdb
ROCKSDB_LDFLAGS = $(shell grep PLATFORM_LDFLAGS deps/rocksdb/make_config.mk | cut -d '=' -f2 | sed s/-ljemalloc// | sed s/-ltcmalloc//)

STATIC_FLAGS=-static -static-libstdc++ -static-libgcc

# These are put into libvg.
OBJ =
OBJ += $(OBJ_DIR)/gssw_aligner.o
OBJ += $(OBJ_DIR)/vg.o 
OBJ += $(OBJ_DIR)/index.o
OBJ += $(OBJ_DIR)/mem.o
OBJ += $(OBJ_DIR)/mapper.o
OBJ += $(OBJ_DIR)/region.o
OBJ += $(OBJ_DIR)/progress_bar.o
OBJ += $(OBJ_DIR)/vg_set.o
OBJ += $(OBJ_DIR)/utility.o
OBJ += $(OBJ_DIR)/path.o
OBJ += $(OBJ_DIR)/alignment.o
OBJ += $(OBJ_DIR)/edit.o
OBJ += $(OBJ_DIR)/sha1.o
OBJ += $(OBJ_DIR)/json2pb.o
OBJ += $(OBJ_DIR)/entropy.o
OBJ += $(OBJ_DIR)/pileup.o
OBJ += $(OBJ_DIR)/caller.o
OBJ += $(OBJ_DIR)/call2vcf.o
OBJ += $(OBJ_DIR)/genotyper.o
OBJ += $(OBJ_DIR)/genotypekit.o
OBJ += $(OBJ_DIR)/position.o
OBJ += $(OBJ_DIR)/deconstructor.o
OBJ += $(OBJ_DIR)/vectorizer.o
OBJ += $(OBJ_DIR)/sampler.o
OBJ += $(OBJ_DIR)/filter.o
OBJ += $(OBJ_DIR)/readfilter.o
OBJ += $(OBJ_DIR)/ssw_aligner.o
OBJ += $(OBJ_DIR)/bubbles.o
OBJ += $(OBJ_DIR)/translator.o
OBJ += $(OBJ_DIR)/version.o
OBJ += $(OBJ_DIR)/banded_global_aligner.o
OBJ += $(OBJ_DIR)/multipath_alignment.o
OBJ += $(OBJ_DIR)/phased_genome.o
OBJ += $(OBJ_DIR)/constructor.o
OBJ += $(OBJ_DIR)/progressive.o
OBJ += $(OBJ_DIR)/flow_sort.o
OBJ += $(OBJ_DIR)/homogenizer.o
OBJ += $(OBJ_DIR)/path_index.o
OBJ += $(OBJ_DIR)/phase_duplicator.o
OBJ += $(OBJ_DIR)/snarls.o
OBJ += $(OBJ_DIR)/feature_set.o
OBJ += $(OBJ_DIR)/simplifier.o
OBJ += $(OBJ_DIR)/chunker.o
OBJ += $(OBJ_DIR)/vcf_buffer.o
OBJ += $(OBJ_DIR)/variant_adder.o
OBJ += $(OBJ_DIR)/name_mapper.o
OBJ += $(OBJ_DIR)/graph_synchronizer.o
OBJ += $(OBJ_DIR)/vg_algorithms.o
OBJ += $(OBJ_DIR)/nested_traversal_finder.o
OBJ += $(OBJ_DIR)/option.o
OBJ += $(OBJ_DIR)/multipath_mapper.o
OBJ += $(OBJ_DIR)/haplotype_extracter.o
OBJ += $(OBJ_DIR)/gamsorter.o

# These aren't put into libvg. But they do go into the main vg binary to power its self-test.
UNITTEST_OBJ =
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/driver.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/distributions.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/genotypekit.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/readfilter.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/banded_global_aligner.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/pinned_alignment.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/vg.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/multipath_alignment.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/phased_genome.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/constructor.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/flow_sort_test.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/srpe_filter.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/phase_duplicator.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/snarls.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/feature_set.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/mapping.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/alignment.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/aligner.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/chunker.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/xg.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/vcf_buffer.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/path_index.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/vg_algorithms.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/union_find.o
UNITTEST_OBJ += $(UNITTEST_OBJ_DIR)/variant_adder.o

# These aren't put into libvg, but they provide subcommand implementations for the vg bianry
SUBCOMMAND_OBJ =
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/subcommand.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/construct_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/add_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/sift_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/srpe_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/simplify_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/index_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/mod_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/annotate_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/chunk_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/snarls_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/explode_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/call_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/genotype_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/trace_main.o 
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/gamsort_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/msga_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/map_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/align_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/surject_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/mpmap_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/find_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/sim_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/view_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/stats_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/vectorize_main.o
SUBCOMMAND_OBJ += $(SUBCOMMAND_OBJ_DIR)/filter_main.o

RAPTOR_DIR:=deps/raptor
PROTOBUF_DIR:=deps/protobuf
GPERF_DIR:=deps/gperftools
SDSL_DIR:=deps/sdsl-lite
SNAPPY_DIR:=deps/snappy
ROCKSDB_DIR:=deps/rocksdb
GCSA2_DIR:=deps/gcsa2
PROGRESS_BAR_DIR:=deps/progress_bar
FASTAHACK_DIR:=deps/fastahack
FERMI_DIR:=deps/fermi-lite
HTSLIB_DIR:=deps/htslib
VCFLIB_DIR:=deps/vcflib
XG_DIR:=deps/xg
GSSW_DIR:=deps/gssw
SPARSEHASH_DIR:=deps/sparsehash
SHA1_DIR:=deps/sha1
DYNAMIC_DIR:=deps/DYNAMIC
SSW_DIR:=deps/ssw/src
STATIC_FLAGS=-static -static-libstdc++ -static-libgcc

# Dependencies that go into libvg's archive
LIB_DEPS =
LIB_DEPS += $(LIB_DIR)/libprotobuf.a
LIB_DEPS += $(LIB_DIR)/libsdsl.a
LIB_DEPS += $(LIB_DIR)/libssw.a
LIB_DEPS += $(LIB_DIR)/libsnappy.a
LIB_DEPS += $(LIB_DIR)/librocksdb.a
LIB_DEPS += $(LIB_DIR)/libgcsa2.a
LIB_DEPS += $(OBJ_DIR)/Fasta.o
LIB_DEPS += $(LIB_DIR)/libhts.a
LIB_DEPS += $(LIB_DIR)/libxg.a
LIB_DEPS += $(LIB_DIR)/libvcflib.a
LIB_DEPS += $(LIB_DIR)/libgssw.a
LIB_DEPS += $(LIB_DIR)/libvcfh.a
LIB_DEPS += $(LIB_DIR)/libgfakluge.a
LIB_DEPS += $(LIB_DIR)/libsupbub.a
LIB_DEPS += $(LIB_DIR)/libsonlib.a
LIB_DEPS += $(LIB_DIR)/libpinchesandcacti.a
LIB_DEPS += $(LIB_DIR)/libraptor2.a
LIB_DEPS += $(LIB_DIR)/libfml.a

# common dependencies to build before all vg src files
DEPS = $(LIB_DEPS)
DEPS += $(CPP_DIR)/vg.pb.h
DEPS += $(INC_DIR)/gcsa/gcsa.h
DEPS += $(INC_DIR)/lru_cache.h
DEPS += $(INC_DIR)/dynamic.hpp
DEPS += $(INC_DIR)/sparsehash/sparse_hash_map
DEPS += $(INC_DIR)/gfakluge.hpp
DEPS += $(INC_DIR)/globalDefs.hpp
DEPS += $(INC_DIR)/sha1.hpp
DEPS += $(INC_DIR)/progress_bar.hpp

ifneq ($(shell uname -s),Darwin)
	DEPS += $(LIB_DIR)/libtcmalloc_minimal.a
	LD_LIB_FLAGS += -ltcmalloc_minimal
endif

.PHONY: clean get-deps test set-path static docs .pre-build

$(BIN_DIR)/vg: $(LIB_DIR)/libvg.a $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o $(BIN_DIR)/vg $(OBJ_DIR)/main.o $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) -lvg $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

# We also want to present the xg binary, which tests can use.
$(BIN_DIR)/xg: $(LIB_DIR)/libxg.a $(OBJ_DIR)/xg-main.o
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o $(BIN_DIR)/xg $(OBJ_DIR)/xg-main.o $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

static: $(OBJ_DIR)/main.o $(OBJ) $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/vg $(OBJ_DIR)/main.o $(OBJ) $(UNITTEST_OBJ) $(SUBCOMMAND_OBJ) $(STATIC_FLAGS) $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

$(LIB_DIR)/libvg.a: $(OBJ) $(DEPS)
	ar rs $@ $(OBJ) $(LIB_DEPS)

get-deps:
	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake libtool jq samtools curl unzip redland-utils librdf-dev cmake pkg-config wget bc gtk-doc-tools raptor2-utils rasqal-utils bison flex libgoogle-perftools-dev


test: $(BIN_DIR)/vg $(LIB_DIR)/libvg.a test/build_graph $(BIN_DIR)/shuf $(BIN_DIR)/xg
	. ./source_me.sh && cd test && $(MAKE)

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
	+. ./source_me.sh && cd $(PROTOBUF_DIR) && ./autogen.sh && ./configure --prefix="$(CWD)" && $(MAKE) && $(MAKE) install && export PATH=$(CWD)/bin:$$PATH

test/build_graph: test/build_graph.cpp $(LIB_DIR)/libvg.a $(CPP_DIR)/vg.pb.h $(SRC_DIR)/json2pb.h $(SRC_DIR)/vg.hpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o test/build_graph test/build_graph.cpp $(LD_INCLUDE_FLAGS) -lvg $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

# remove annoying large alloc messages from tcmalloc
$(GPERF_DIR)/src/tcmalloc.cc.bak:
	cp $(GPERF_DIR)/src/tcmalloc.cc $(GPERF_DIR)/src/tcmalloc.cc.bak
	sed 's/printer.printf("tcmalloc: large alloc/return; printer.printf("tcmalloc: large alloc/' $(GPERF_DIR)/src/tcmalloc.cc.bak >$(GPERF_DIR)/src/tcmalloc.cc

$(LIB_DIR)/libtcmalloc_minimal.a: $(GPERF_DIR)/src/tcmalloc.cc.bak
	+. ./source_me.sh && cd $(GPERF_DIR) && ./autogen.sh && ./configure --prefix=`pwd` && $(MAKE) && $(MAKE) install && cp -r lib/* $(CWD)/$(LIB_DIR)/ && cp -r bin/* $(CWD)/$(BIN_DIR)/ && cp -r include/* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libsdsl.a: $(SDSL_DIR)/lib/*.cpp $(SDSL_DIR)/include/sdsl/*.hpp
	+. ./source_me.sh && cd $(SDSL_DIR) && ./install.sh $(CWD)

$(LIB_DIR)/libssw.a:
	+. ./source_me.sh && cd $(SSW_DIR) && $(MAKE) && ar rs $(CWD)/$(LIB_DIR)/libssw.a ssw.o ssw_cpp.o && cp ssw_cpp.h ssw.h $(CWD)/$(LIB_DIR)

$(LIB_DIR)/libsnappy.a:
	+. ./source_me.sh && cd $(SNAPPY_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) && $(MAKE) && $(MAKE) install

$(LIB_DIR)/librocksdb.a: $(LIB_DIR)/libsnappy.a $(ROCKSDB_DIR)/db/*.cc $(ROCKSDB_DIR)/db/*.h
	+. ./source_me.sh && cd $(ROCKSDB_DIR) && $(ROCKSDB_PORTABLE) DISABLE_JEMALLOC=1 $(MAKE) static_lib && mv librocksdb.a $(CWD)/${LIB_DIR}/ && cp -r include/* $(CWD)/$(INC_DIR)/

$(INC_DIR)/gcsa/gcsa.h: $(LIB_DIR)/libgcsa2.a
$(LIB_DIR)/libgcsa2.a: $(LIB_DIR)/libsdsl.a $(wildcard $(GCSA2_DIR)/*.cpp) $(wildcard $(GCSA2_DIR)/*.hpp)
	+. ./source_me.sh && cd $(GCSA2_DIR) && cat Makefile | grep -v VERBOSE_STATUS_INFO >Makefile.quiet && $(MAKE) -f Makefile.quiet libgcsa2.a && mv libgcsa2.a $(CWD)/$(LIB_DIR) && cp -r include/gcsa $(CWD)/$(INC_DIR)/

$(INC_DIR)/progress_bar.hpp: $(PROGRESS_BAR_DIR)/progress_bar.hpp
	+cp $(PROGRESS_BAR_DIR)/progress_bar.hpp $(CWD)/$(INC_DIR)

$(OBJ_DIR)/progress_bar.o:
	+cd $(PROGRESS_BAR_DIR) && $(MAKE) && cp progress_bar.o $(CWD)/$(OBJ_DIR)

$(OBJ_DIR)/Fasta.o:
	+cd $(FASTAHACK_DIR) && $(MAKE) && mv Fasta.o $(CWD)/$(OBJ_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

$(LIB_DIR)/libhts.a:
	+cd $(HTSLIB_DIR) && $(MAKE) lib-static && mv libhts.a $(CWD)/$(LIB_DIR) && cp *.h $(CWD)/$(INC_DIR) && cp -r htslib $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libxg.a: $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libprotobuf.a $(CPP_DIR)/vg.pb.o $(INC_DIR)/dynamic.hpp $(XG_DIR)/src/xg.cpp $(XG_DIR)/src/xg.hpp $(INC_DIR)/sparsehash/sparse_hash_map 
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c $(XG_DIR)/src/xg.cpp -o $(CWD)/$(OBJ_DIR)/xg.o $(LD_INCLUDE_FLAGS) -I$(INC_DIR)/dynamic && ar rs $(CWD)/$(LIB_DIR)/libxg.a $(CWD)/$(OBJ_DIR)/xg.o $(CWD)/$(CPP_DIR)/vg.pb.o && cp $(XG_DIR)/src/*.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libvcflib.a: $(VCFLIB_DIR)/src/*.cpp $(VCFLIB_DIR)/src/*.hpp $(VCFLIB_DIR)/intervaltree/*.cpp $(VCFLIB_DIR)/intervaltree/*.h $(VCFLIB_DIR)/tabixpp/*.cpp $(VCFLIB_DIR)/tabixpp/*.hpp $(VCFLIB_DIR)/tabixpp/htslib/*.c $(VCFLIB_DIR)/tabixpp/htslib/htslib/*.h $(VCFLIB_DIR)/tabixpp/htslib/cram/*.c $(VCFLIB_DIR)/tabixpp/htslib/cram/*.h
	+. ./source_me.sh && cd $(VCFLIB_DIR) && $(MAKE) libvcflib.a && cp lib/* $(CWD)/$(LIB_DIR)/ && cp include/* $(CWD)/$(INC_DIR)/ && cp intervaltree/*.h $(CWD)/$(INC_DIR)/ && cp src/*.h* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libgssw.a: $(GSSW_DIR)/src/gssw.c $(GSSW_DIR)/src/gssw.h
	+cd $(GSSW_DIR) && $(MAKE) && cp lib/* $(CWD)/$(LIB_DIR)/ && cp obj/* $(CWD)/$(OBJ_DIR) && cp src/*.h $(CWD)/$(INC_DIR)

$(INC_DIR)/lru_cache.h:
	+cd $(DEP_DIR)/lru_cache && $(MAKE) && cp *.h* $(CWD)/$(INC_DIR)/

$(INC_DIR)/dynamic.hpp:
	+cat $(DYNAMIC_DIR)/include/dynamic.hpp | sed 's%<internal/%<dynamic/%' >$(INC_DIR)/dynamic.hpp && cp -r $(CWD)/$(DYNAMIC_DIR)/include/internal $(CWD)/$(INC_DIR)/dynamic

$(INC_DIR)/sparsehash/sparse_hash_map:
	+cd $(SPARSEHASH_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) && $(MAKE) && $(MAKE) install

$(LIB_DIR)/libvcfh.a:
	+cd $(DEP_DIR)/libVCFH && $(MAKE) && cp libvcfh.a $(CWD)/$(LIB_DIR)/ && cp vcfheader.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libgfakluge.a: $(INC_DIR)/gfakluge.hpp
	+cd $(DEP_DIR)/gfakluge && $(MAKE) && cp libgfakluge.a $(CWD)/$(LIB_DIR)/

$(INC_DIR)/gfakluge.hpp:
	cp $(DEP_DIR)/gfakluge/src/gfakluge.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libsupbub.a: $(LIB_DIR)/libsdsl.a $(INC_DIR)/globalDefs.hpp
	+. ./source_me.sh && cd $(DEP_DIR)/superbubbles && $(MAKE) && cp libsupbub.a $(CWD)/$(LIB_DIR)/

$(LIB_DIR)/libsonlib.a: $(CWD)/$(DEP_DIR)/sonLib/C/inc/*.h $(CWD)/$(DEP_DIR)/sonLib/C/impl/*.c
	+cd $(DEP_DIR)/sonLib && kyotoTycoonLib="" $(MAKE) && cp lib/sonLib.a $(CWD)/$(LIB_DIR)/libsonlib.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

$(LIB_DIR)/libpinchesandcacti.a: $(LIB_DIR)/libsonlib.a $(CWD)/$(DEP_DIR)/pinchesAndCacti/inc/*.h $(CWD)/$(DEP_DIR)/pinchesAndCacti/impl/*.c
	+cd $(DEP_DIR)/pinchesAndCacti && $(MAKE) && cd $(CWD)/$(DEP_DIR)/sonLib && cp lib/stPinchesAndCacti.a $(CWD)/$(LIB_DIR)/libpinchesandcacti.a && cp lib/3EdgeConnected.a $(CWD)/$(LIB_DIR)/lib3edgeconnected.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

$(INC_DIR)/globalDefs.hpp: $(LIB_DIR)/libsdsl.a
	cp $(DEP_DIR)/superbubbles/*.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libraptor2.a:
	+cd $(RAPTOR_DIR)/build && cmake .. && $(MAKE) && cp src/libraptor2.a $(CWD)/$(LIB_DIR) && mkdir -p $(CWD)/$(INC_DIR)/raptor2 && cp src/*.h $(CWD)/$(INC_DIR)/raptor2

$(INC_DIR)/sha1.hpp: $(SHA1_DIR)/sha1.hpp
	+cp $(SHA1_DIR)/*.h* $(CWD)/$(INC_DIR)/

$(OBJ_DIR)/sha1.o: $(SHA1_DIR)/sha1.cpp $(SHA1_DIR)/sha1.hpp
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS)

$(LIB_DIR)/libfml.a:
	cd $(FERMI_DIR) && $(MAKE) && cp *.h $(CWD)/$(INC_DIR)/ && cp libfml.a $(CWD)/$(LIB_DIR)/

# Auto-versioning
$(INC_DIR)/vg_git_version.hpp: .git
	echo "#define VG_GIT_VERSION \"$(shell git describe --always --tags || echo unknown)\"" > $@

# Not important if .git isn't real
.git:

###################################
## VG source code compilation begins here
####################################

include/stream.hpp: src/stream.hpp
	cp src/stream.hpp include/stream.hpp

$(CPP_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.cc

$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h 

$(CPP_DIR)/vg.pb.h: $(LIB_DIR)/libprotobuf.a bin/protoc $(SRC_DIR)/vg.proto
	+. ./source_me.sh && ./bin/protoc $(SRC_DIR)/vg.proto --proto_path=$(SRC_DIR) --cpp_out=cpp
	+cp $@ $(INC_DIR)

# For most of the .o files we just specify dependencies and use the pattern rule.

$(OBJ): $(DEPS)

$(OBJ_DIR)/vg.o: $(SRC_DIR)/vg.cpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/gssw_aligner.hpp $(DEPS)

$(OBJ_DIR)/banded_global_aligner.o: $(SRC_DIR)/banded_global_aligner.cpp $(SRC_DIR)/banded_global_aligner.hpp $(DEPS)

$(OBJ_DIR)/gssw_aligner.o: $(SRC_DIR)/gssw_aligner.cpp $(SRC_DIR)/gssw_aligner.hpp $(DEPS)

$(OBJ_DIR)/ssw_aligner.o: $(SRC_DIR)/ssw_aligner.cpp $(SRC_DIR)/ssw_aligner.hpp $(DEPS)

$(OBJ_DIR)/vg_set.o: $(SRC_DIR)/vg_set.cpp $(SRC_DIR)/vg_set.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/index.hpp $(DEPS)

$(OBJ_DIR)/mapper.o: $(SRC_DIR)/mapper.cpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/mem.hpp $(SRC_DIR)/vg.hpp $(DEPS)

$(OBJ_DIR)/mem.o: $(SRC_DIR)/mem.cpp $(SRC_DIR)/mem.hpp $(SRC_DIR)/vg.hpp $(DEPS)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(INC_DIR)/stream.hpp $(DEPS) $(SRC_DIR)/utility.hpp $(INC_DIR)/globalDefs.hpp $(SRC_DIR)/bubbles.hpp $(SRC_DIR)/genotyper.hpp $(SRC_DIR)/readfilter.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/index.hpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp

$(OBJ_DIR)/region.o: $(SRC_DIR)/region.cpp $(SRC_DIR)/region.hpp $(DEPS)

$(OBJ_DIR)/index.o: $(SRC_DIR)/index.cpp $(SRC_DIR)/index.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(INC_DIR)/stream.hpp $(DEPS)

$(OBJ_DIR)/utility.o: $(SRC_DIR)/utility.cpp $(SRC_DIR)/utility.hpp $(DEPS)

$(OBJ_DIR)/path.o: $(SRC_DIR)/path.cpp $(SRC_DIR)/path.hpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/utility.hpp $(DEPS)

$(OBJ_DIR)/edit.o: $(SRC_DIR)/edit.cpp $(SRC_DIR)/edit.hpp $(DEPS)

$(OBJ_DIR)/alignment.o: $(SRC_DIR)/alignment.cpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/alignment.hpp $(SRC_DIR)/edit.hpp $(SRC_DIR)/edit.cpp $(INC_DIR)/stream.hpp $(DEPS)

$(OBJ_DIR)/json2pb.o: $(SRC_DIR)/json2pb.cpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/bin2ascii.h $(DEPS)

$(OBJ_DIR)/entropy.o: $(SRC_DIR)/entropy.cpp $(SRC_DIR)/entropy.hpp $(DEPS)

$(OBJ_DIR)/pileup.o: $(SRC_DIR)/pileup.cpp $(SRC_DIR)/pileup.hpp $(INC_DIR)/stream.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/json2pb.h $(DEPS)

$(OBJ_DIR)/caller.o: $(SRC_DIR)/caller.cpp $(SRC_DIR)/caller.hpp $(SRC_DIR)/option.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(INC_DIR)/stream.hpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/pileup.hpp $(SRC_DIR)/snarls.hpp $(DEPS)

$(OBJ_DIR)/call2vcf.o: $(SRC_DIR)/call2vcf.cpp $(SRC_DIR)/caller.hpp $(SRC_DIR)/option.hpp $(SRC_DIR)/nodeside.hpp $(SRC_DIR)/nodetraversal.hpp $(SRC_DIR)/snarls.hpp $(SRC_DIR)/path_index.hpp $(SRC_DIR)/distributions.hpp $(SRC_DIR)/genotypekit.hpp $(SRC_DIR)/snarls.hpp $(DEPS)

$(OBJ_DIR)/genotyper.o: $(SRC_DIR)/genotyper.cpp $(SRC_DIR)/genotyper.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/path_index.hpp $(SRC_DIR)/json2pb.h $(DEPS) $(INC_DIR)/sparsehash/sparse_hash_map $(SRC_DIR)/bubbles.hpp $(SRC_DIR)/distributions.hpp $(SRC_DIR)/utility.hpp

$(OBJ_DIR)/genotypekit.o: $(SRC_DIR)/genotypekit.cpp $(SRC_DIR)/genotypekit.hpp $(DEPS) $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/distributions.hpp $(SRC_DIR)/snarls.hpp $(DEPS)

$(OBJ_DIR)/position.o: $(SRC_DIR)/position.cpp $(SRC_DIR)/position.hpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/json2pb.h $(DEPS)

$(OBJ_DIR)/version.o: $(SRC_DIR)/version.cpp $(SRC_DIR)/version.hpp $(INC_DIR)/vg_git_version.hpp

$(OBJ_DIR)/progressive.o: $(SRC_DIR)/progressive.cpp $(SRC_DIR)/progressive.hpp $(DEPS)

$(OBJ_DIR)/path_index.o: $(SRC_DIR)/path_index.cpp $(SRC_DIR)/path_index.hpp $(DEPS)

## TODO vcflib build loses variant.h
$(OBJ_DIR)/deconstructor.o: $(SRC_DIR)/deconstructor.cpp $(SRC_DIR)/deconstructor.hpp $(LIB_DIR)/libvcfh.a $(SRC_DIR)/bubbles.hpp $(DEPS)

$(OBJ_DIR)/vectorizer.o: $(SRC_DIR)/vectorizer.cpp $(SRC_DIR)/vectorizer.hpp $(DEPS)

$(OBJ_DIR)/sampler.o: $(SRC_DIR)/sampler.cpp $(SRC_DIR)/sampler.hpp $(SRC_DIR)/path.hpp $(DEPS)

$(OBJ_DIR)/filter.o: $(SRC_DIR)/filter.cpp $(SRC_DIR)/filter.hpp $(DEPS)

$(OBJ_DIR)/readfilter.o: $(SRC_DIR)/readfilter.cpp $(SRC_DIR)/readfilter.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(INC_DIR)/stream.hpp $(DEPS)

$(OBJ_DIR)/homogenizer.o: $(SRC_DIR)/homogenizer.cpp $(SRC_DIR)/homogenizer.hpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/bubbles.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/gssw_aligner.hpp $(SRC_DIR)/filter.hpp $(DEPS)

$(OBJ_DIR)/bubbles.o: $(SRC_DIR)/bubbles.cpp $(SRC_DIR)/bubbles.hpp $(DEPS)

$(OBJ_DIR)/translator.o: $(SRC_DIR)/translator.cpp $(SRC_DIR)/translator.hpp $(DEPS)

$(OBJ_DIR)/constructor.o: $(SRC_DIR)/constructor.cpp $(SRC_DIR)/constructor.hpp $(SRC_DIR)/vcf_buffer.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/name_mapper.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(OBJ_DIR)/chunker.o: $(SRC_DIR)/chunker.cpp $(SRC_DIR)/chunker.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(OBJ_DIR)/vcf_buffer.o: $(SRC_DIR)/vcf_buffer.cpp $(SRC_DIR)/vcf_buffer.hpp $(DEPS)

$(OBJ_DIR)/variant_adder.o: $(SRC_DIR)/variant_adder.cpp $(SRC_DIR)/variant_adder.hpp $(SRC_DIR)/name_mapper.hpp $(SRC_DIR)/vcf_buffer.hpp $(SRC_DIR)/graph_synchronizer.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/mapper.hpp $(DEPS)

$(OBJ_DIR)/name_mapper.o: $(SRC_DIR)/name_mapper.cpp $(DEPS)

$(OBJ_DIR)/graph_synchronizer.o: $(SRC_DIR)/graph_synchronizer.cpp $(SRC_DIR)/graph_synchronizer.hpp $(SRC_DIR)/algorithms/vg_algorithms.hpp $(SRC_DIR)/path_index.hpp $(SRC_DIR)/vg.hpp $(DEPS)

# We also build the main file from xg, so we can offer it as a command
# TODO: wrap it as a vg subcommand?
$(OBJ_DIR)/xg-main.o: $(XG_DIR)/src/main.cpp $(XG_DIR)/src/xg.hpp $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c $(XG_DIR)/src/main.cpp -o $(OBJ_DIR)/xg-main.o $(LD_INCLUDE_FLAGS) -I$(INC_DIR)/dynamic

$(OBJ_DIR)/phased_genome.o: $(SRC_DIR)/phased_genome.cpp $(SRC_DIR)/phased_genome.hpp $(DEPS)

$(OBJ_DIR)/multipath_alignment.o: $(SRC_DIR)/multipath_alignment.cpp $(SRC_DIR)/multipath_alignment.hpp $(DEPS)

$(OBJ_DIR)/snarls.o: $(SRC_DIR)/snarls.cpp $(SRC_DIR)/snarls.hpp $(DEPS)
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

$(OBJ_DIR)/flow_sort.o: $(SRC_DIR)/flow_sort.cpp $(SRC_DIR)/flow_sort.hpp $(SRC_DIR)/vg.hpp $(DEPS)

$(OBJ_DIR)/srpe.o: $(SRC_DIR)/srpe.cpp $(SRC_DIR)/srpe.hpp $(OBJ_DIR)/filter.o $(LIB_DIR)/libvcflib.a $(DEPS) $(LIB_DIR)/libfml.a

$(OBJ_DIR)/gamsorter.o: $(SRC_DIR)/gamsorter.cpp $(SRC_DIR)/gamsorter.hpp $(DEPS)

$(OBJ_DIR)/path_index.o: $(SRC_DIR)/path_index.cpp $(SRC_DIR)/path_index.hpp $(DEPS)

$(OBJ_DIR)/phase_duplicator.o: $(SRC_DIR)/phase_duplicator.cpp $(SRC_DIR)/phase_duplicator.hpp $(SRC_DIR)/types.hpp $(DEPS)

$(OBJ_DIR)/feature_set.o: $(SRC_DIR)/feature_set.cpp $(SRC_DIR)/feature_set.hpp $(SRC_DIR)/types.hpp $(DEPS)

$(OBJ_DIR)/simplifier.o: $(SRC_DIR)/simplifier.cpp $(SRC_DIR)/simplifier.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/feature_set.hpp $(SRC_DIR)/path.hpp $(SRC_DIR)/path_index.hpp $(SRC_DIR)/distributions.hpp $(DEPS)

$(OBJ_DIR)/nested_traversal_finder.o: $(SRC_DIR)/nested_traversal_finder.cpp $(SRC_DIR)/nested_traversal_finder.hpp $(SRC_DIR)/genotypekit.hpp $(SRC_DIR)/snarls.hpp $(DEPS)

$(OBJ_DIR)/option.o: $(SRC_DIR)/option.cpp $(SRC_DIR)/option.hpp $(DEPS)

$(OBJ_DIR)/multipath_mapper.o: $(SRC_DIR)/multipath_mapper.cpp $(SRC_DIR)/multipath_mapper.hpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/mem.hpp $(DEPS)

$(OBJ_DIR)/haplotype_extracter.o: $(SRC_DIR)/haplotype_extracter.cpp $(SRC_DIR)/haplotype_extracter.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/json2pb.h $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libxg.a $(CPP_DIR)/vg.pb.h
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

### Algorithms ###

$(OBJ_DIR)/vg_algorithms.o: $(ALGORITHMS_SRC_DIR)/vg_algorithms.cpp $(ALGORITHMS_SRC_DIR)/vg_algorithms.hpp $(DEPS)
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

###################################
## VG unit test compilation begins here
####################################

$(UNITTEST_OBJ_DIR)/driver.o: $(UNITTEST_SRC_DIR)/driver.cpp $(UNITTEST_SRC_DIR)/driver.hpp $(UNITTEST_SRC_DIR)/catch.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/distributions.o: $(UNITTEST_SRC_DIR)/distributions.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/distributions.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/srpe_filter.o: $(UNITTEST_SRC_DIR)/srpe_filter.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(DEPS)
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(UNITTEST_OBJ_DIR)/banded_global_aligner.o: $(UNITTEST_SRC_DIR)/banded_global_aligner.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/banded_global_aligner.hpp $(SRC_DIR)/gssw_aligner.hpp $(SRC_DIR)/gssw_aligner.cpp $(DEPS)

$(UNITTEST_OBJ_DIR)/pinned_alignment.o: $(UNITTEST_SRC_DIR)/pinned_alignment.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/gssw_aligner.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/genotypekit.o: $(UNITTEST_SRC_DIR)/genotypekit.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/genotypekit.hpp $(SRC_DIR)/snarls.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/readfilter.o: $(UNITTEST_SRC_DIR)/readfilter.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/readfilter.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/multipath_alignment.o: $(UNITTEST_SRC_DIR)/multipath_alignment.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/multipath_alignment.hpp $(SRC_DIR)/multipath_alignment.cpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/mem.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/phased_genome.o: $(UNITTEST_SRC_DIR)/phased_genome.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/phased_genome.hpp $(SRC_DIR)/phased_genome.cpp $(SRC_DIR)/snarls.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/vg.o: $(UNITTEST_SRC_DIR)/vg.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/constructor.o: $(UNITTEST_SRC_DIR)/constructor.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/constructor.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/name_mapper.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/flow_sort_test.o: $(UNITTEST_SRC_DIR)/flow_sort_test.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/phase_duplicator.o: $(UNITTEST_SRC_DIR)/phase_duplicator.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/phase_duplicator.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/feature_set.o: $(UNITTEST_SRC_DIR)/feature_set.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/feature_set.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/mapping.o: $(UNITTEST_SRC_DIR)/mapping.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/path.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/alignment.o: $(UNITTEST_SRC_DIR)/alignment.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/alignment.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/aligner.o: $(UNITTEST_SRC_DIR)/aligner.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/gssw_aligner.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/snarls.o: $(UNITTEST_SRC_DIR)/snarls.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/snarls.hpp $(DEPS)
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

$(UNITTEST_OBJ_DIR)/chunker.o: $(UNITTEST_SRC_DIR)/chunker.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/chunker.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/xg.o: $(UNITTEST_SRC_DIR)/xg.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/vcf_buffer.o: $(UNITTEST_SRC_DIR)/vcf_buffer.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/vcf_buffer.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/path_index.o: $(UNITTEST_SRC_DIR)/path_index.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/path_index.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/vg_algorithms.o: $(UNITTEST_SRC_DIR)/vg_algorithms.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(ALGORITHMS_SRC_DIR)/vg_algorithms.hpp $(ALGORITHMS_SRC_DIR)/vg_algorithms.cpp $(DEPS)

$(UNITTEST_OBJ_DIR)/union_find.o: $(UNITTEST_SRC_DIR)/union_find.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(UNITTEST_OBJ_DIR)/variant_adder.o: $(UNITTEST_SRC_DIR)/variant_adder.cpp $(UNITTEST_SRC_DIR)/catch.hpp $(SRC_DIR)/variant_adder.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/name_mapper.hpp $(DEPS)

###################################
## VG subcommand compilation begins here
####################################

$(SUBCOMMAND_OBJ_DIR)/subcommand.o: $(SUBCOMMAND_SRC_DIR)/subcommand.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/construct_main.o: $(SUBCOMMAND_SRC_DIR)/construct_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/constructor.hpp $(SRC_DIR)/name_mapper.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/simplify_main.o: $(SUBCOMMAND_SRC_DIR)/simplify_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/simplifier.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/feature_set.hpp $(SRC_DIR)/path.hpp $(SRC_DIR)/path_index.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/homogenize_main.o: $(SUBCOMMAND_SRC_DIR)/homogenize_main.cpp $(OBJ_DIR)/homogenizer.o $(OBJ_DIR)/filter.o $(OBJ_DIR)/mapper.hpp $(SRC_DIR)/mem.hpp $(OBJ_DIR)/bubbles.o $(OBJ_DIR)/vg.o $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/sift_main.o: $(SUBCOMMAND_SRC_DIR)/sift_main.cpp $(OBJ_DIR)/filter.o $(OBJ_DIR)/vg.o $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/srpe_main.o: $(SUBCOMMAND_SRC_DIR)/srpe_main.cpp $(OBJ_DIR)/srpe.o $(OBJ_DIR)/filter.o $(OBJ_DIR)/mapper.o $(OBJ_DIR)/vg.o $(DEPS)


$(SUBCOMMAND_OBJ_DIR)/gamsort_main.o: $(SUBCOMMAND_SRC_DIR)/gamsort_main.cpp $(OBJ_DIR)/gamsorter.o $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/index_main.o: $(SUBCOMMAND_SRC_DIR)/index_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/index.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/vg_set.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/path_index.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/mod_main.o: $(SUBCOMMAND_SRC_DIR)/mod_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/annotate_main.o: $(SUBCOMMAND_SRC_DIR)/annotate_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/chunk_main.o: $(SUBCOMMAND_SRC_DIR)/chunk_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/add_main.o: $(SUBCOMMAND_SRC_DIR)/add_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vcf_buffer.hpp $(SRC_DIR)/variant_adder.hpp $(SRC_DIR)/name_mapper.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/snarls_main.o: $(SUBCOMMAND_SRC_DIR)/snarls_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/distributions.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/explode_main.o: $(SUBCOMMAND_SRC_DIR)/explode_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/call_main.o: $(SUBCOMMAND_SRC_DIR)/call_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/option.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/utility.hpp $(SRC_DIR)/caller.hpp $(SRC_DIR)/nodeside.hpp $(SRC_DIR)/nodetraversal.hpp $(SRC_DIR)/snarls.hpp $(SRC_DIR)/path_index.hpp $(SRC_DIR)/distributions.hpp $(SRC_DIR)/progressive.hpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/pileup.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/genotype_main.o: $(SUBCOMMAND_SRC_DIR)/genotype_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/genotyper.hpp $(SRC_DIR)/genotyper.cpp $(SRC_DIR)/snarls.hpp $(SRC_DIR)/path_index.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/msga_main.o: $(SUBCOMMAND_SRC_DIR)/msga_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/mem.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/map_main.o: $(SUBCOMMAND_SRC_DIR)/map_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/mem.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/mpmap_main.o: $(SUBCOMMAND_SRC_DIR)/mpmap_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/multipath_mapper.hpp $(SRC_DIR)/mem.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/align_main.o: $(SUBCOMMAND_SRC_DIR)/align_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/gssw_aligner.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/surject_main.o: $(SUBCOMMAND_SRC_DIR)/surject_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/mapper.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/find_main.o: $(SUBCOMMAND_SRC_DIR)/find_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/sim_main.o: $(SUBCOMMAND_SRC_DIR)/sim_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/mapper.hpp $(SRC_DIR)/sampler.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/view_main.o: $(SUBCOMMAND_SRC_DIR)/view_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/stats_main.o: $(SUBCOMMAND_SRC_DIR)/stats_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/vectorize_main.o: $(SUBCOMMAND_SRC_DIR)/vectorize_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/vectorizer.hpp $(SRC_DIR)/mapper.hpp $(DEPS)

$(SUBCOMMAND_OBJ_DIR)/filter_main.o: $(SUBCOMMAND_SRC_DIR)/filter_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/stream.hpp $(SRC_DIR)/readfilter.hpp $(DEPS)


########################
## Pattern Rules
########################

# Define a default rule for building objects from CPP files
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS)
$(SUBCOMMAND_OBJ_DIR)/%.o : $(SUBCOMMAND_DIR)/%.cpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS)
$(UNITTEST_OBJ_DIR)/%.o : $(UNITTEST_SRC_DIR)/%.cpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS)
        
# Protobuf stuff builds into its same directory
$(CPP_DIR)/%.o : $(CPP_DIR)/%.cc
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS)

$(SUBCOMMAND_OBJ_DIR)/trace_main.o: $(SUBCOMMAND_SRC_DIR)/trace_main.cpp $(SUBCOMMAND_SRC_DIR)/subcommand.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/haplotype_extracter.hpp $(DEPS)
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS) $(ROCKSDB_LDFLAGS)

###################################
## VG source code compilation ends here
####################################



.pre-build:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(UNITTEST_OBJ_DIR) ]; then mkdir -p $(UNITTEST_OBJ_DIR); fi
	if [ ! -d $(SUBCOMMAND_OBJ_DIR) ]; then mkdir -p $(SUBCOMMAND_OBJ_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(CPP_DIR) ]; then mkdir -p $(CPP_DIR); fi

# run .pre-build before we make anything at all.
-include .pre-build

# for rebuilding just vg
clean-vg:
	$(RM) -r $(BIN_DIR)/vg
	$(RM) -r $(UNITTEST_OBJ_DIR)/*.o
	$(RM) -r $(SUBCOMMAND_OBJ_DIR)/*.o
	$(RM) -r $(OBJ_DIR)/*.o
	$(RM) -r $(CPP_DIR)/*.o $(CPP_DIR)/*.cc $(CPP_DIR)/*.h

clean:
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(UNITTEST_OBJ_DIR)
	$(RM) -r $(SUBCOMMAND_OBJ_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r $(CPP_DIR)
	$(RM) -r share/
	cd $(DEP_DIR) && cd protobuf && $(MAKE) clean
	cd $(DEP_DIR) && cd xg && $(MAKE) clean
	cd $(DEP_DIR) && cd vcflib && $(MAKE) clean
	cd $(DEP_DIR) && cd sparsehash && $(MAKE) clean
	cd $(DEP_DIR) && cd htslib && $(MAKE) clean
	cd $(DEP_DIR) && cd fastahack && $(MAKE) clean
	cd $(DEP_DIR) && cd gcsa2 && $(MAKE) clean
	cd $(DEP_DIR) && cd gssw && $(MAKE) clean
	cd $(DEP_DIR) && cd progress_bar && $(MAKE) clean
	cd $(DEP_DIR) && cd sdsl-lite && ./uninstall.sh || true
	cd $(DEP_DIR) && cd libVCFH && $(MAKE) clean
	cd $(DEP_DIR) && cd rocksdb && $(MAKE) clean
	cd $(DEP_DIR) && cd superbubbles && $(MAKE) clean
	cd $(DEP_DIR) && cd gfakluge && $(MAKE) clean
	rm -Rf $(RAPTOR_DIR)/build/*
	## TODO vg source code
	## TODO LRU_CACHE
	## TODO bash-tap
	## TODO sha1

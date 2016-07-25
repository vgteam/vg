DEP_DIR:=./deps
SRC_DIR:=src
BIN_DIR:=bin
OBJ_DIR:=obj
LIB_DIR:=lib
INC_DIR:=include
CPP_DIR:=cpp

EXE:=vg

CXX:=g++
CXXFLAGS:=-O3 -msse4.1 -fopenmp -std=c++11 -ggdb

CWD:=$(shell pwd)


LD_INCLUDE_FLAGS:=-I$(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR) -I$(CWD)/$(CPP_DIR) -I$(CWD)/$(INC_DIR)/dynamic -I$(CWD)/$(INC_DIR)/sonLib
LD_LIB_FLAGS:= -ggdb -L$(CWD)/$(LIB_DIR) -lvcflib -lgssw -lssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2 -lgcsa2 -lxg -ldivsufsort -ldivsufsort64 -lvcfh -lgfakluge -lraptor2 -lsupbub -lsdsl -lpinchesandcacti -l3edgeconnected  -lsonlib

ifeq ($(shell uname -s),Darwin)
    # We may need libraries from Macports
    # TODO: where does Homebrew keep libraries?
    ifeq ($(shell if [ -d /opt/local/lib ];then echo 1;else echo 0;fi), 1)
       LD_LIB_FLAGS += -L/opt/local/lib
    endif
    ifeq ($(shell if [ -d /usr/local/lib ];then echo 1;else echo 0;fi), 1)
       LD_LIB_FLAGS += -L/usr/local/lib
    endif
    ROCKSDB_PORTABLE=PORTABLE=1 # needed to build rocksdb without weird assembler options
else
    # Not on OS X, we can have librt
    LD_LIB_FLAGS += -lrt
endif

STATIC_FLAGS=-static -static-libstdc++ -static-libgcc

OBJ:=$(OBJ_DIR)/gssw_aligner.o $(OBJ_DIR)/vg.o cpp/vg.pb.o $(OBJ_DIR)/index.o $(OBJ_DIR)/mapper.o $(OBJ_DIR)/region.o $(OBJ_DIR)/progress_bar.o $(OBJ_DIR)/vg_set.o $(OBJ_DIR)/utility.o $(OBJ_DIR)/path.o $(OBJ_DIR)/alignment.o $(OBJ_DIR)/edit.o $(OBJ_DIR)/sha1.o $(OBJ_DIR)/json2pb.o $(OBJ_DIR)/entropy.o $(OBJ_DIR)/pileup.o $(OBJ_DIR)/caller.o $(OBJ_DIR)/genotyper.o $(OBJ_DIR)/position.o $(OBJ_DIR)/deconstructor.o $(OBJ_DIR)/vectorizer.o $(OBJ_DIR)/sampler.o $(OBJ_DIR)/filter.o $(OBJ_DIR)/ssw_aligner.o $(OBJ_DIR)/bubbles.o $(OBJ_DIR)/translator.o $(OBJ_DIR)/haplotypes.o

RAPTOR_DIR:=deps/raptor
PROTOBUF_DIR:=deps/protobuf
SDSL_DIR:=deps/sdsl-lite
SNAPPY_DIR:=deps/snappy
ROCKSDB_DIR:=deps/rocksdb
GCSA2_DIR:=deps/gcsa2
PROGRESS_BAR_DIR:=deps/progress_bar
FASTAHACK_DIR:=deps/fastahack
HTSLIB_DIR:=deps/htslib
VCFLIB_DIR:=deps/vcflib
XG_DIR:=deps/xg
GSSW_DIR:=deps/gssw
SPARSEHASH_DIR:=deps/sparsehash
SHA1_DIR:=deps/sha1
DYNAMIC_DIR:=deps/DYNAMIC
SSW_DIR:=deps/ssw/src
STATIC_FLAGS=-static -static-libstdc++ -static-libgcc

.PHONY: clean get-deps test set-path static .pre-build

$(BIN_DIR)/vg: $(LIB_DIR)/libvg.a $(OBJ_DIR)/main.o
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o $(BIN_DIR)/vg $(OBJ_DIR)/main.o -lvg $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

static: $(OBJ_DIR)/main.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/vg $(OBJ_DIR)/main.o $(OBJ) $(STATIC_FLAGS) $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(LIB_DIR)/libvg.a: $(OBJ)
	ar rs $@ $^

get-deps:
	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake libtool jq samtools curl unzip redland-utils librdf-dev cmake pkg-config wget bc gtk-doc-tools raptor2-utils rasqal-utils bison flex

test: $(BIN_DIR)/vg $(LIB_DIR)/libvg.a test/build_graph $(BIN_DIR)/shuf
	. ./source_me.sh && cd test && $(MAKE)

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
.rebuild-protobuf: deps/protobuf/src/google/protobuf/*.cc
	+. ./source_me.sh && cd $(PROTOBUF_DIR) && ./autogen.sh && ./configure --prefix="$(CWD)" && $(MAKE) && $(MAKE) install && export PATH=$(CWD)/bin:$$PATH

test/build_graph: test/build_graph.cpp $(LIB_DIR)/libvg.a $(CPP_DIR)/vg.pb.h $(SRC_DIR)/json2pb.h $(SRC_DIR)/vg.hpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o test/build_graph test/build_graph.cpp $(LD_INCLUDE_FLAGS) -lvg $(LD_LIB_FLAGS)

deps: $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libgssw.a $(LIB_DIR)/libgcsa2.a $(LIB_DIR)/libsnappy.a $(LIB_DIR)/libvcflib.a $(INC_DIR)/sparsehash/sparse_hash_map $(OBJ_DIR)/sha1.o $(LIB_DIR)/librocksdb.a $(LIB_DIR)/libhts.a $(LIB_DIR)/libxg.a $(LIB_DIR)/libsonlib.a $(LIB_DIR)/libpinchesandcacti.a 

$(LIB_DIR)/libsdsl.a:
	+. ./source_me.sh && cd $(SDSL_DIR) && ./install.sh $(CWD)

$(LIB_DIR)/libssw.a:
	+. ./source_me.sh && cd $(SSW_DIR) && $(MAKE) && ar rs $(CWD)/$(LIB_DIR)/libssw.a ssw.o ssw_cpp.o && cp ssw_cpp.h ssw.h $(CWD)/$(LIB_DIR)

$(LIB_DIR)/libsnappy.a:
	+. ./source_me.sh && cd $(SNAPPY_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) && $(MAKE) && $(MAKE) install

$(LIB_DIR)/librocksdb.a: $(LIB_DIR)/libsnappy.a
	+. ./source_me.sh && cd $(ROCKSDB_DIR) && $(ROCKSDB_PORTABLE) $(MAKE) static_lib && mv librocksdb.a $(CWD)/${LIB_DIR}/ && cp -r include/* $(CWD)/$(INC_DIR)/

$(INC_DIR)/gcsa.h: $(LIB_DIR)/libgcsa2.a
$(LIB_DIR)/libgcsa2.a: $(LIB_DIR)/libsdsl.a
	+. ./source_me.sh && cd $(GCSA2_DIR) && cat Makefile | grep -v VERBOSE_STATUS_INFO >Makefile.quiet && $(MAKE) -f Makefile.quiet libgcsa2.a && mv libgcsa2.a $(CWD)/$(LIB_DIR) && cp *.h* $(CWD)/$(INC_DIR)/

$(OBJ_DIR)/progress_bar.o:
	+cd $(PROGRESS_BAR_DIR) && $(MAKE) && cp progress_bar.o $(CWD)/$(OBJ_DIR) && cp *.h* $(CWD)/$(INC_DIR)

$(OBJ_DIR)/Fasta.o:
	+cd $(FASTAHACK_DIR) && $(MAKE) && mv Fasta.o $(CWD)/$(OBJ_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

$(LIB_DIR)/libhts.a:
	+cd $(HTSLIB_DIR) && $(MAKE) lib-static && mv libhts.a $(CWD)/$(LIB_DIR) && cp *.h $(CWD)/$(INC_DIR) && cp -r htslib $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libxg.a: $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libprotobuf.a $(CPP_DIR)/vg.pb.o $(INC_DIR)/dynamic.hpp $(XG_DIR)/src/xg.cpp
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c $(XG_DIR)/src/xg.cpp -o $(CWD)/$(OBJ_DIR)/xg.o $(LD_INCLUDE_FLAGS) -I$(INC_DIR)/dynamic && ar rs $(CWD)/$(LIB_DIR)/libxg.a $(CWD)/$(OBJ_DIR)/xg.o $(CWD)/$(CPP_DIR)/vg.pb.o && cp $(XG_DIR)/src/*.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libvcflib.a:
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
	cp $(DEP_DIR)/gfakluge/gfakluge.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libsupbub.a: $(LIB_DIR)/libsdsl.a $(INC_DIR)/globalDefs.hpp
	+. ./source_me.sh && cd $(DEP_DIR)/superbubbles && $(MAKE) && cp libsupbub.a $(CWD)/$(LIB_DIR)/

$(LIB_DIR)/libsonlib.a:
	+cd $(DEP_DIR)/sonLib && kyotoTycoonLib="" $(MAKE) && cp lib/sonLib.a $(CWD)/$(LIB_DIR)/libsonlib.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

$(LIB_DIR)/libpinchesandcacti.a: $(LIB_DIR)/libsonlib.a
	+cd $(DEP_DIR)/pinchesAndCacti && $(MAKE) && cd $(CWD)/$(DEP_DIR)/sonLib && cp lib/stPinchesAndCacti.a $(CWD)/$(LIB_DIR)/libpinchesandcacti.a && cp lib/3EdgeConnected.a $(CWD)/$(LIB_DIR)/lib3edgeconnected.a && mkdir -p $(CWD)/$(INC_DIR)/sonLib && cp lib/*.h $(CWD)/$(INC_DIR)/sonLib

$(INC_DIR)/globalDefs.hpp: $(LIB_DIR)/libsdsl.a
	cp $(DEP_DIR)/superbubbles/*.hpp $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libraptor2.a:
	+cd $(RAPTOR_DIR)/build && cmake .. && $(MAKE) && cp src/libraptor2.a $(CWD)/$(LIB_DIR) && mkdir -p $(CWD)/$(INC_DIR)/raptor2 && cp src/*.h $(CWD)/$(INC_DIR)/raptor2

$(INC_DIR)/sha1.hpp: $(OBJ_DIR)/sha1.o
$(OBJ_DIR)/sha1.o: $(SHA1_DIR)/sha1.cpp $(SHA1_DIR)/sha1.hpp
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) && cp $(SHA1_DIR)/*.h* $(CWD)/$(INC_DIR)/

# Auto-versioning
$(INC_DIR)/vg_git_version.hpp: .git
	echo "#define VG_GIT_VERSION \"$(shell git describe --always --tags || echo unknown)\"" > $@

###################################
## VG source code compilation begins here
####################################

include/stream.hpp:
	cp src/stream.hpp include/stream.hpp

$(CPP_DIR)/vg.pb.o: $(CPP_DIR)/vg.pb.cc
	+. ./source_me.sh && g++ -O3 -msse4.1 -fopenmp -std=c++11 -c -o cpp/vg.pb.o cpp/vg.pb.cc $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h 

$(CPP_DIR)/vg.pb.h: $(LIB_DIR)/libprotobuf.a bin/protoc src/vg.proto
	+. ./source_me.sh && ./bin/protoc $(SRC_DIR)/vg.proto --proto_path=$(SRC_DIR) --cpp_out=cpp
	+cp $@ $(INC_DIR)

$(OBJ_DIR)/vg.o: $(SRC_DIR)/vg.cpp $(SRC_DIR)/vg.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libvcflib.a $(FASTAHACK_DIR)/Fasta.o $(LIB_DIR)/libgssw.a $(INC_DIR)/sparsehash/sparse_hash_map $(INC_DIR)/lru_cache.h $(INC_DIR)/stream.hpp $(LIB_DIR)/libprotobuf.a $(OBJ_DIR)/progress_bar.o $(INC_DIR)/gcsa.h $(INC_DIR)/sha1.hpp $(LIB_DIR)/libgfakluge.a $(LIB_DIR)/libsupbub.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libraptor2.a $(LIB_DIR)/libpinchesandcacti.a $(LIB_DIR)/libsonlib.a
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/gssw_aligner.o: $(SRC_DIR)/gssw_aligner.cpp $(SRC_DIR)/gssw_aligner.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libgssw.a $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/libvcflib.a $(LIB_DIR)/libhts.a ${INC_DIR}/sha1.hpp
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/ssw_aligner.o: $(SRC_DIR)/ssw_aligner.cpp $(SRC_DIR)/ssw_aligner.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libssw.a $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/libvcflib.a
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/vg_set.o: $(SRC_DIR)/vg_set.cpp $(SRC_DIR)/vg_set.hpp $(SRC_DIR)/vg.hpp $(OBJ_DIR)/index.o $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libgssw.a $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libxg.a $(INC_DIR)/dynamic.hpp
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/mapper.o: $(SRC_DIR)/mapper.cpp $(SRC_DIR)/mapper.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libxg.a
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(INC_DIR)/vg_git_version.hpp $(LIB_DIR)/libvcflib.a $(OBJ_DIR)/Fasta.o $(LIB_DIR)/libgssw.a $(INC_DIR)/stream.hpp $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/librocksdb.a $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libxg.a $(INC_DIR)/gcsa.h $(LIB_DIR)/libhts.a $(INC_DIR)/sha1.hpp $(OBJ_DIR)/progress_bar.o $(INC_DIR)/lru_cache.h $(LIB_DIR)/libvcfh.a $(LIB_DIR)/libgfakluge.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libpinchesandcacti.a $(LIB_DIR)/libsonlib.a $(INC_DIR)/globalDefs.hpp $(SRC_DIR)/bubbles.hpp $(SRC_DIR)/genotyper.hpp
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/region.o: $(SRC_DIR)/region.cpp $(SRC_DIR)/region.hpp $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/index.o: $(SRC_DIR)/index.cpp $(SRC_DIR)/index.hpp $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/librocksdb.a $(LIB_DIR)/libxg.a $(LIB_DIR)/libsnappy.a $(CPP_DIR)/vg.pb.h
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/utility.o: $(SRC_DIR)/utility.cpp $(SRC_DIR)/utility.hpp $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(CPP_DIR)/vg.pb.h $(INC_DIR)/sha1.hpp $(CPP_DIR)/vg.pb.h
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/path.o: $(SRC_DIR)/path.cpp $(SRC_DIR)/path.hpp $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(CPP_DIR)/vg.pb.h $(OBJ_DIR)/utility.o $(LIB_DIR)/libgssw.a
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/edit.o: $(SRC_DIR)/edit.cpp $(SRC_DIR)/edit.hpp $(LIB_DIR)/libprotobuf.a $(CPP_DIR)/vg.pb.h
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/alignment.o: $(SRC_DIR)/alignment.cpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/alignment.hpp $(LIB_DIR)/libhts.a $(LIB_DIR)/libprotobuf.a  $(INC_DIR)/sparsehash/sparse_hash_map $(SRC_DIR)/edit.hpp $(SRC_DIR)/edit.cpp
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_INCLUDE_FLAGS)

$(OBJ_DIR)/json2pb.o: $(SRC_DIR)/json2pb.cpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/bin2ascii.h $(LIB_DIR)/libprotobuf.a
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/entropy.o: $(SRC_DIR)/entropy.cpp $(SRC_DIR)/entropy.hpp
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/pileup.o: $(SRC_DIR)/pileup.cpp $(SRC_DIR)/pileup.hpp $(CPP_DIR)/vg.pb.h $(INC_DIR)/stream.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/json2pb.h $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(LIB_DIR)/libgcsa2.a
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/caller.o: $(SRC_DIR)/caller.cpp $(SRC_DIR)/caller.hpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/vg.hpp $(INC_DIR)/stream.hpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/pileup.hpp $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/caller.o $(SRC_DIR)/caller.cpp $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/genotyper.o: $(SRC_DIR)/genotyper.cpp $(SRC_DIR)/genotyper.hpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/vg.hpp $(INC_DIR)/stream.hpp $(SRC_DIR)/json2pb.h $(LIB_DIR)/libprotobuf.a $(INC_DIR)/sparsehash/sparse_hash_map $(SRC_DIR)/bubbles.hpp $(SRC_DIR)/distributions.hpp $(SRC_DIR)/utility.hpp
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/genotyper.o $(SRC_DIR)/genotyper.cpp $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/position.o: $(SRC_DIR)/position.cpp $(SRC_DIR)/position.hpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/vg.hpp $(SRC_DIR)/json2pb.h $(LIB_DIR)/libprotobuf.a
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

## TODO vcflib build loses variant.h
$(OBJ_DIR)/deconstructor.o: $(SRC_DIR)/deconstructor.cpp $(SRC_DIR)/deconstructor.hpp $(LIB_DIR)/libvcflib.a $(LIB_DIR)/librocksdb.a $(INC_DIR)/gcsa.h $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libxg.a $(LIB_DIR)/libvcfh.a $(SRC_DIR)/bubbles.hpp
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/vectorizer.o: $(SRC_DIR)/vectorizer.cpp $(SRC_DIR)/vectorizer.hpp $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libxg.a $(LIB_DIR)/libprotobuf.a
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/sampler.o: $(SRC_DIR)/sampler.cpp $(SRC_DIR)/sampler.hpp $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libxg.a $(LIB_DIR)/libprotobuf.a
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/filter.o: $(SRC_DIR)/filter.cpp $(SRC_DIR)/filter.hpp $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libxg.a $(CPP_DIR)/vg.pb.h
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/bubbles.o: $(SRC_DIR)/bubbles.cpp $(SRC_DIR)/bubbles.hpp $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libxg.a $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libpinchesandcacti.a $(LIB_DIR)/libsonlib.a
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/translator.o: $(SRC_DIR)/translator.cpp $(SRC_DIR)/translator.hpp $(LIB_DIR)/libprotobuf.a $(CPP_DIR)/vg.pb.h
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)
	
$(OBJ_DIR)/haplotypes.o: $(SRC_DIR)/haplotypes.cpp $(SRC_DIR)/haplotypes.hpp $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libxg.a $(CPP_DIR)/vg.pb.h
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

.pre-build:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(CPP_DIR) ]; then mkdir -p $(CPP_DIR); fi

# run .pre-build before we make anything at all.
-include .pre-build

# for rebuilding just vg
clean-vg:
	$(RM) -r $(BIN_DIR)/vg
	$(RM) -r $(OBJ_DIR)/*
	$(RM) -r $(CPP_DIR)/*

clean:
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
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
	rm -R $(RAPTOR_DIR)/build/*
## TODO vg source code
## TODO LRU_CACHE
## TODO bash-tap
## TODO sha1

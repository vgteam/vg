.PHONY: all clean test get-deps

#TODO we need a way to manage paths
## perhaps a target fix_paths
## or make install_local / make install

## Compilers and compiler flags
#CC:=icc
#CFLAGS:="-O3 -xHost -xSSE4.1 -openmp"
#CXX:=icpc
#CXXFLAGS:="-O3 -xHost -openmp -xSSE4.1"
CXX:=g++
CXXFLAGS:=-O0 -msse4 -fopenmp -std=c++11
## Build directories for dependencies
VG_DIR:=vg
VCFLIB_DIR:=vcflib
PROTOBUF_DIR:=protobuf
ROCKS_DIR:=rocksdb
HTSLIB_DIR:=htslib
SPARSEHASH_DIR:=sparsehash
GCSA_DIR:=gcsa2
GSSW_DIR:=gssw
XG_DIR:=xg
SDSL_DIR:=sdsl-lite
SNAPPY_DIR:=snappy
FASTHACK_DIR:=fastahack
LRU_CACHE_DIR:=lru_cache
SHA1_DIR:=sha1
PROGRESS_BAR_DIR:=progress_bar

## Options
ROCKSDB_PORTABLE:=

## Copy all deps into this dir once built
DEP_DIR:=deps

BIN_DIR:=bin
LIB_DIR:=lib
INC_DIR:=include
SRC_DIR:=src

## Test dir
TEST_DIR:=test

## Place protobuf output here
PROTO_OUT_DIR:=cpp
## Vars for Object files
VG_SRC:=$(wildcard $(VG_DIR)/*.c)
VG_OBJ:=$(patsubst %.cpp, %.o, $(VG_SRC))

## Vars for lib files
LIB_VCFLIB:=libvcflib.a
LIB_GSSW:=libgssw.a
LIB_PROTOBUF:=libprotobuf.a
LIB_SNAPPY:=libsnappy.a
LIB_ROCKS:=librocksdb.a
LIB_SPARSEHASH:=sparsehash/build/include/sparsehash/sparse_hash_map
LIB_HTSLIB:=ibhts.a
LIB_GCSA:=libgcsa2.a
LIB_XG:=libxg.a
LIB_SDSLLITE:=build/include/sdsl
LIB_FASTAHACK:=

## Vars for Executable files
EXE:=$(BIN_DIR)/vg
## Linker flags
#INCLUDES:= -I./ -I$(PROTOBUF_OUT_DIR) -I$(VCFLIB_DIR) -I$(FASTAHACK_DIR) -I$(GSSW_DIR)/src -I$(PROTOBUF_DIR)/build/include -I$(ROCKS_DIR)/include -I$(LRU_CACHE_DIR) -I$(SHA1_DIR) -I$(XG_DIR) -I$(SDSL_DIR)/install/include -I$(GCSA_DIR) -I$(HTSLIB_DIR)
INCLUDES:= -I./include/
LIBS:= -L./lib/ -lvcflib -lgssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2 -lgcsa2 -lxg -lsdsl -ldivsufsort -ldivsufsort64 -lvg

## Switch for building on OS X
SYS=$(shell uname -s)
#ifeq ($(SYS),Darwin)
#
#else
#
#endif

all: $(EXE) $(LIB_DIR)/libvg.a


test: $(EXE) $(VG_DIR)/libvg.a build_graph
	cd $(TEST_DIR) && $(MAKE)

build_graph: $(TEST_DIR)/build_graph.cpp $(LIB_DIR)/libvg.a
	$(CXX) $(CXXFLAGS) $(TEST_DIR)/build_graph.cpp $(INCLUDES) $(LIBS) -o $(TEST_DIR)/build_graph

## Install fresh jq and link into path
## try ports for other dependencies
get-deps-mac:
	sudo port install bison bzip2 zlibc protobufcompiler

get-deps:
	sudo apt-get install -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev libtool curl unzip cmake pkg-config wget bc

deps: $(LIB_ROCKS) $(LIB_SNAPPY) $(LIB_GCSA) $(LIB_SDSLLITE) $(LIB_XG) $(LIB_HTSLIB)
	if [ ! -d $(DEP_DIR) ]; then mkdir $(DEP_DIR); fi
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	if [ ! -f $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	mv $(VCFLIB_DIR)/$(LIB_VCFLIB) $(DEP_DIR)
	mv $(SPARSEHASH_DIR)/$(LIB_SPARSEHASH) $(DEP_DIR)
	mv $(HTSLIB_DIR)/$(LIB_HTSLIB) $(DEP_DIR)
	mv $(PROTOBUF_DIR)/$(LIB_PROTOBUF) $(DEP_DIR)
	mv $(SNAPPY_DIR)/$(LIB_SNAPPY) $(DEP_DIR)
	mv $(HTSLIB_DIR)/$(LIB_HTSLIB) $(DEP_DIR)

$(LIB_HTSLIB):
	cd $(HTSLIB_DIR) && $(MAKE) lib-static

$(LIB_PROTOBUF): $(PROTOBUF_DIR)/src/google/protobuf/*cc $(PROTOBUF_DIR)/src/google/protobuf/*h
	cd $(PROTOBUF_DIR) && ./autogen.sh && ./configure --prefix=`pwd` && $(MAKE) && $(MAKE) install
	

$(LIB_ROCKS): $(ROCKS_DIR)/include/rocksdb/*.h $(ROCKS_DIR)/db/*.c $(ROCKS_DIR)/db/*.cc $(ROCKS_DIR)/db/*.h
	cd rocksdb && $(ROCKSDB_PORTABLE) $(MAKE) static_lib

$(LIB_XG): $(LIB_SDSLLITE) $(LIB_PROTOBUF) $(PROTO_OUT_DIR)/vg.pb.h
	cd $(XG_DIR) && $(CXX) $(CXXFLAGS) -I$(PROTO_OUT_DIR) -c -o xg.o xg.cpp 
	ar rs $(LIB_XG) xg.o

#$(XG_DIR)/xg.o: $(XG_DIR)/xg.cpp $(XG_DIR)/xg.hpp $(LIBSDSL) $(VG_DIR)/cpp/vg.pb.h
#	cd $(XG_DIR) && $(CXX) $(CXXFLAGS) -c -o xg.o xg.cpp $(INCLUDES)

$(LIB_SNAPPY): $(SNAPPY_DIR)/*cc $(SNAPPY_DIR)/*h
	cd $(SNAPPY_DIR) && ./autogen.sh && ./configure --prefix=`pwd` && $(MAKE) && $(MAKE) install

##$(GCSA_DIR)/*.cpp $(GCSA_DIR)/*.h
## $(SDSL_DIR)/lib/*.cpp $(SDSL_DIR)/include/*.hpp
$(LIB_DIR)/$(LIB_SDSLLITE): 
	$(DEP_DIR)/$(SDSL_DIR)/install.sh `pwd`
	touch $(LIB_DIR)/$(LIB_SDSLLITE)

$(LIB_GCSA): $(LIB_SDSLLITE)
	cd $(DEP_DIR)/$(GCSA_DIR) && $(MAKE) $(LIB_GCSA) && mv $(LIB_GCSA) ../../lib/
	touch $(LIB_DIR)/$(LIB_GCSA)

$(LIB_FASTAHACK): $(DEP_DIR)/$(FASTAHACK_DIR)/Fasta.h fastahack/Fasta.cpp
	cd $(FASTAHACK_DIR) && $(MAKE)

#progress_bar/progress_bar.cpp progress_bar/progress_bar.hpp
$(PROGRESS_BAR_DIR)/progress_bar.o: 
	cd $(PROGRESS_BAR_DIR) && $(MAKE)

$(PROTO_OUT_DIR)/vg.pb.cc: $(PROTO_OUT_DIR)/vg.pb.h
$(PROTO_OUT_DIR)/vg.pb.h: $(LIB_PROTOBUF)
	mkdir -p cpp
	$(PROTOBUF_DIR)/bin/protoc $(VG_DIR)/vg.proto --cpp_out=cpp






vg.o: $(PROTO_OUT_DIR)/vg.pb.h $(LIB_VCFLIB) $(FASTAHACK_DIR)/Fasta.o $(LIB_GSSW) $(LIB_SPARSEHASH) $(LRU_CACHE_DIR)/lru_cache.h stream.hpp $(LIB_PROTOBUF) $(LIB_SDSLLITE)
	$(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

gssw_aligner.o: gssw_aligner.cpp gssw_aligner.hpp cpp/vg.pb.h $(LIB_GSSW) $(LIB_PROTOBUF) $(LIB_SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o gssw_aligner.o gssw_aligner.cpp $(INCLUDES)

vg_set.o: vg_set.cpp vg_set.hpp vg.hpp index.hpp cpp/vg.pb.h $(LIBGSSW) $(LIBPROTOBUF) $(SPARSEHASH) $(SDSLLITE)
	$(CXX) $(CXXFLAGS) -c -o vg_set.o vg_set.cpp $(INCLUDES)

mapper.o: mapper.cpp mapper.hpp cpp/vg.pb.h $(LIBPROTOBUF) $(SPARSEHASH) $(SDSLLITE)
	$(CXX) $(CXXFLAGS) -c -o mapper.o mapper.cpp $(INCLUDES)

main.o: main.cpp $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) stream.hpp  $(LIBPROTOBUF) $(SPARSEHASH) $(SDSLLITE)
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

region.o: region.cpp region.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o region.o region.cpp $(INCLUDES)

index.o: index.cpp index.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o index.o index.cpp $(INCLUDES)

utility.o: utility.cpp utility.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o utility.o utility.cpp $(INCLUDES)

path.o: path.cpp path.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o path.o path.cpp $(INCLUDES)

edit.o: edit.cpp edit.hpp $(LIBPROTOBUF)
	$(CXX) $(CXXFLAGS) -c -o edit.o edit.cpp $(INCLUDES)

alignment.o: alignment.cpp alignment.hpp $(LIBHTS)  $(LIBPROTOBUF) $(SPARSEHASH) edit.hpp edit.cpp
	$(CXX) $(CXXFLAGS) -c -o alignment.o alignment.cpp $(INCLUDES)

sha1/sha1.o: sha1/sha1.cpp sha1/sha1.hpp
	$(CXX) $(CXXFLAGS) -c -o sha1/sha1.o sha1/sha1.cpp $(INCLUDES)

json2pb.o: json2pb.cpp json2pb.h bin2ascii.h $(LIBPROTOBUF)
	$(CXX) $(CXXFLAGS) -c -o json2pb.o json2pb.cpp $(INCLUDES)

entropy.o: entropy.cpp entropy.hpp
	$(CXX) $(CXXFLAGS) -c -o entropy.o entropy.cpp $(INCLUDES)

pileup.o: pileup.cpp pileup.hpp cpp/vg.pb.h vg.hpp stream.hpp json2pb.h $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o pileup.o pileup.cpp $(INCLUDES)

caller.o: caller.cpp caller.hpp cpp/vg.pb.h vg.hpp stream.hpp json2pb.h pileup.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o caller.o caller.cpp $(INCLUDES)





$(EXE): deps
	$(CXX) $(CXXFLAGS) -o $(EXE) $(INCLUDES) $(LIBS) $(STATIC_FLAGS) && mv $(VG_EXE) ../


clean:
	$(RM) $(DEP_DIR)
	cd $(VCFLIB_DIR) && $(MAKE) clean
	cd $(SNAPPY_DIR) && $(MAKE) clean
	cd $(PROTOBUF_DIR) && $(MAKE) clean
	cd $(GCSA_DIR) && $(MAKE) clean
	cd $(SPARSEHASH_DIR) && $(MAKE) clean
	cd $(ROCKS_DIR) && $(MAKE) clean
	cd $(SDSL_DIR) && ./uninstall.sh `pwd`

clobber-vg:
	$(RM) $(VG_OBJ)	
	$(RM) $(VG_EXE)
	$(RM) $(VG_LB)
	$(RM) $(PROTO_OUT_DIR)/*

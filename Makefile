DEP_DIR:=./deps
SRC_DIR:=src
BIN_DIR:=bin
OBJ_DIR:=obj
LIB_DIR:=lib
INC_DIR:=include

CXX:=g++
CXXFLAGS:=-O0 -msse4 -fopenmp -std=c++11

CWD:=$(shell pwd)

LD_INCLUDE_FLAGS:=-I$(INC_DIR)
LD_LIB_FLAGS:= -L$(LIB_DIR) -lvcflib -lgssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2 -lgcsa2 -lxg -lsdsl -ldivsufsort -ldivsufsort64 -lvg

OBJ:=gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o mapper.o region.o progress_bar/progress_bar.o vg_set.o utility.o path.o alignment.o edit.o sha1/sha1.o json2pb.o entropy.o pileup.o caller.o

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


$(shell ./source_me.sh)

.PHONY: clean get-deps

all: vg libvg.a

vg: deps

deps: protobuf sdsl-lite gssw gcsa2 snappy vcflib sparsehash sha1 rocksdb htslib

protobuf:
	cd $(PROTOBUF_DIR) && ./autogen.sh && ./configure --prefix="$(CWD)" && make -j 8 && make install

sdsl-lite:
	cd $(SDSL_DIR) && ./install.sh $(CWD)

snappy:
	cd $(SNAPPY_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) && $(MAKE) && $(MAKE) install

rocksdb:
	cd $(ROCKSDB_DIR) && $(MAKE) static_lib && mv librocksdb.a $(CWD)/lib

gcsa2: sdsl-lite
	cd $(GCSA2_DIR) && $(MAKE) && mv libgcsa2.a $(CWD)/lib

progress_bar:
	cd $(PROGRESS_BAR_DIR) && $(MAKE) && mv progress_bar.o $(OBJ_DIR) && cp progress_bar.hpp $(INC_DIR)

fastahack:
	cd $(FASTAHACK_DIR) && make && mv Fasta.o $(CWD)/$(OBJ_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

htslib:
	cd $(HTSLIB_DIR) && $(MAKE) lib-static && mv libhts.a $(CWD)/$(LIB_DIR) && cp *.h $(CWD)/$(INC_DIR)

xg: sdsl-lite protobuf
	cd $(XG_DIR) && $(MAKE) all && cp obj/xg.o $(CWD)/$(OBJ_DIR) && cp lib/libxg.a $(CWD)/$(LIB_DIR) && cp src/*.hpp $(CWD)/$(INC_DIR)

vcflib:
	cd $(VCFLIB_DIR) && $(MAKE) && cp lib/* $(CWD)/$(LIB_DIR) && cp src/*.h $(CWD)/$(INC_DIR)


.PHONY: all clean test get-deps

CXX=g++ -std=c++11 -fopenmp -g
CXXFLAGS=-O3
pb2json=pb2json/libpb2json.a
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
LIBGSSW=gssw/src/libgssw.a
LIBSNAPPY=snappy/libsnappy.a
LIBROCKSDB=rocksdb/librocksdb.a
SPARSEHASH=sparsehash/build/include/sparsehash/sparse_hash_map
INCLUDES=-I./ -Ipb2json -Icpp -I$(VCFLIB)/src -I$(VCFLIB) -Ifastahack -Igssw/src -Irocksdb/include -Iprogress_bar/ -Isparsehash/build/include
LDFLAGS=-L./ -Lpb2json -Lvcflib -Lgssw/src -Lsnappy -Lrocksdb -Lprogressbar -lpb2json -lvcflib -lgssw -lprotobuf -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2
LIBS=gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o mapper.o region.o progress_bar/progress_bar.o

all: vg libvg.a

get-deps:
	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev

test: vg libvg.a test/build_graph
	cd test && $(MAKE)

test/build_graph: test/build_graph.cpp
	$(CXX) test/build_graph.cpp $(INCLUDES) -lvg $(LDFLAGS) -o test/build_graph

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

$(LIBSNAPPY): snappy/*cc snappy/*h
	cd snappy && mkdir -p build && ./autogen.sh && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install
	cp snappy/build/lib/libsnappy.a snappy/

$(LIBROCKSDB): rocksdb/include/rocksdb/*.h rocksdb/db/*.c rocksdb/db/*.cc rocksdb/db/*.h
	cd rocksdb && $(MAKE) static_lib

progress_bar/progress_bar.o: progress_bar/progress_bar.cpp progress_bar/progress_bar.hpp
	cd progress_bar && make

$(pb2json):
	cd pb2json && $(MAKE) libpb2json.a

cpp/vg.pb.cc: cpp/vg.pb.h
cpp/vg.pb.h: vg.proto
	mkdir -p cpp
	protoc vg.proto --cpp_out=cpp

$(LIBVCFLIB): vcflib/src/Variant.h vcflib/src/Variant.cpp
	cd vcflib && $(MAKE) libvcflib.a

$(LIBGSSW): gssw/src/gssw.c gssw/src/gssw.h
	cd gssw/src && $(MAKE) libgssw.a

$(SPARSEHASH):
	cd sparsehash && mkdir -p build && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install

fastahack/Fasta.o: fastahack/Fasta.h fastahack/Fasta.cpp
	cd fastahack && $(MAKE)

cpp/vg.pb.o: cpp/vg.pb.h cpp/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

vg.o: vg.cpp vg.h cpp/vg.pb.h $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

gssw_aligner.o: gssw_aligner.cpp gssw_aligner.h cpp/vg.pb.h $(LIBGSSW)
	$(CXX) $(CXXFLAGS) -c -o gssw_aligner.o gssw_aligner.cpp $(INCLUDES)

mapper.o: mapper.cpp mapper.h cpp/vg.pb.h
	$(CXX) $(CXXFLAGS) -c -o mapper.o mapper.cpp $(INCLUDES)

main.o: main.cpp $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW)
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

region.o: region.cpp region.h
	$(CXX) $(CXXFLAGS) -c -o region.o region.cpp $(INCLUDES)

index.o: index.cpp index.h
	$(CXX) $(CXXFLAGS) -c -o index.o index.cpp $(INCLUDES)

vg: $(LIBS) $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW) $(LIBROCKSDB) $(LIBSNAPPY)
	$(CXX) $(CXXFLAGS) -o vg $(LIBS) $(INCLUDES) $(LDFLAGS)

libvg.a: vg
	ar rs libvg.a gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o mapper.o region.o progress_bar/progress_bar.o

clean-vg:
	rm -f vg
	rm -f cpp/*
	rm -f *.o
	cd progress_bar && make clean

clean: clean-vg
	cd test && $(MAKE) clean
	cd pb2json && $(MAKE) clean
	cd vcflib && $(MAKE) clean
	cd snappy && $(MAKE) clean
	rm -f snappy/libsnappy.a
	cd rocksdb && $(MAKE) clean

.PHONY: all clean test get-deps

CXX=g++ 
CXXFLAGS=-O3 -std=c++11 -fopenmp -g -march=native
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
LIBGSSW=gssw/src/libgssw.a
LIBPROTOBUF=protobuf/libprotobuf.a
LIBSNAPPY=snappy/libsnappy.a
LIBROCKSDB=rocksdb/librocksdb.a
SPARSEHASH=sparsehash/build/include/sparsehash/sparse_hash_map
LIBHTS=htslib/libhts.a
INCLUDES=-I./ -Icpp -I$(VCFLIB)/src -I$(VCFLIB) -Ifastahack -Igssw/src -Iprotobuf/build/include -Irocksdb/include -Iprogress_bar -Isparsehash/build/include -Ilru_cache -Ihtslib -Isha1
LDFLAGS=-L./ -Lvcflib -Lgssw/src -Lprotobuf -Lsnappy -Lrocksdb -Lprogressbar -Lhtslib -lvcflib -lgssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2
LIBS=gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o mapper.o region.o progress_bar/progress_bar.o vg_set.o utility.o path.o json.o alignment.o sha1/sha1.o pb2json.o entropy.o

all: vg libvg.a

get-deps:
	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake libtool jq

test: vg libvg.a test/build_graph
	cd test && $(MAKE)

test/build_graph: test/build_graph.cpp
	$(CXX) test/build_graph.cpp $(INCLUDES) -lvg $(LDFLAGS) -o test/build_graph

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

$(LIBPROTOBUF): protobuf/src/google/protobuf/*cc  protobuf/src/google/protobuf/*h
	cd protobuf && mkdir -p build && ./autogen.sh && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install
	cp protobuf/build/lib/libprotobuf.a protobuf/

$(LIBSNAPPY): snappy/*cc snappy/*h
	cd snappy && mkdir -p build && ./autogen.sh && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install
	cp snappy/build/lib/libsnappy.a snappy/

$(LIBROCKSDB): rocksdb/include/rocksdb/*.h rocksdb/db/*.c rocksdb/db/*.cc rocksdb/db/*.h
	cd rocksdb && $(MAKE) static_lib

progress_bar/progress_bar.o: progress_bar/progress_bar.cpp progress_bar/progress_bar.hpp
	cd progress_bar && make

cpp/vg.pb.cc: cpp/vg.pb.h
cpp/vg.pb.h: vg.proto $(LIBPROTOBUF)
	mkdir -p cpp
	protobuf/build/bin/protoc vg.proto --cpp_out=cpp

$(LIBVCFLIB): vcflib/src/Variant.h vcflib/src/Variant.cpp
	cd vcflib && $(MAKE) libvcflib.a

$(LIBGSSW): gssw/src/gssw.c gssw/src/gssw.h
	cd gssw/src && $(MAKE) libgssw.a

$(SPARSEHASH): sparsehash/build/include/sparsehash/dense_hash_map

sparsehash/build/include/sparsehash/dense_hash_map:
	cd sparsehash && mkdir -p build && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install

$(LIBHTS):
	cd htslib && $(MAKE) lib-static

fastahack/Fasta.o: fastahack/Fasta.h fastahack/Fasta.cpp
	cd fastahack && $(MAKE)

cpp/vg.pb.o: cpp/vg.pb.h cpp/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

vg.o: vg.cpp vg.hpp cpp/vg.pb.h $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) $(SPARSEHASH) lru_cache/lru_cache.h stream.hpp $(LIBPROTOBUF)
	$(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

gssw_aligner.o: gssw_aligner.cpp gssw_aligner.hpp cpp/vg.pb.h $(LIBGSSW) $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o gssw_aligner.o gssw_aligner.cpp $(INCLUDES)

vg_set.o: vg_set.cpp vg_set.hpp vg.hpp index.hpp cpp/vg.pb.h $(LIBGSSW) $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o vg_set.o vg_set.cpp $(INCLUDES)

mapper.o: mapper.cpp mapper.hpp cpp/vg.pb.h $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o mapper.o mapper.cpp $(INCLUDES)

main.o: main.cpp $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) stream.hpp  $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

region.o: region.cpp region.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o region.o region.cpp $(INCLUDES)

index.o: index.cpp index.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o index.o index.cpp $(INCLUDES)

utility.o: utility.cpp utility.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o utility.o utility.cpp $(INCLUDES)

path.o: path.cpp path.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o path.o path.cpp $(INCLUDES)

alignment.o: alignment.cpp alignment.hpp $(LIBHTS)  $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o alignment.o alignment.cpp $(INCLUDES)

json.o: json.cpp json.hpp $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o json.o json.cpp $(INCLUDES)

sha1/sha1.o: sha1/sha1.cpp sha1/sha1.hpp
	$(CXX) $(CXXFLAGS) -c -o sha1/sha1.o sha1/sha1.cpp $(INCLUDES)

pb2json.o: pb2json.cpp pb2json.h $(LIBPROTOBUF)
	$(CXX) $(CXXFLAGS) -c -o pb2json.o pb2json.cpp $(INCLUDES)

entropy.o: entropy.cpp entropy.hpp
	$(CXX) $(CXXFLAGS) -c -o entropy.o entropy.cpp $(INCLUDES)

vg: $(LIBS) $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) $(LIBROCKSDB) $(LIBSNAPPY) $(LIBHTS) $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -o vg $(LIBS) $(INCLUDES) $(LDFLAGS)

libvg.a: vg
	ar rs libvg.a gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o mapper.o region.o progress_bar/progress_bar.o utility.o path.o json.o alignment.o sha1/sha1.o pb2json.o

clean-vg:
	rm -f vg
	rm -f cpp/*
	rm -f *.o
	cd progress_bar && make clean

clean: clean-vg
	cd test && $(MAKE) clean
	cd vcflib && $(MAKE) clean
	cd snappy && $(MAKE) clean
	rm -f snappy/libsnappy.a
	cd protobuf && $(MAKE) clean && rm -rf build
	rm -f protobuf/libprotobuf.a
	cd rocksdb && $(MAKE) clean
	cd sparsehash && $(MAKE) clean && rm -rf build

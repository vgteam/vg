.PHONY: all clean test get-deps

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g # -ffast-math -funroll-loops
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
LIBGSSW=gssw/src/libgssw.a
LIBPROTOBUF=protobuf/libprotobuf.a
LIBSNAPPY=snappy/libsnappy.a
LIBROCKSDB=rocksdb/librocksdb.a
SPARSEHASH=sparsehash/build/include/sparsehash/sparse_hash_map
LIBHTS=htslib/libhts.a
LIBGCSA2=gcsa2/libgcsa2.a
LIBXG=xg/libxg.a
SDSLLITE=sdsl-lite/build/include/sdsl
INCLUDES=-I./ -Icpp -I$(VCFLIB)/src -I$(VCFLIB) -Ifastahack -Igssw/src -Iprotobuf/build/include -Irocksdb/include -Iprogress_bar -Isparsehash/build/include -Ilru_cache -Ihtslib -Isha1 -Isdsl-lite/install/include -Igcsa2 -Ixg
LDFLAGS=-L./ -Lvcflib -Lgssw/src -Lprotobuf -Lsnappy -Lrocksdb -Lprogressbar -Lhtslib -Lgcsa2 -Lsdsl-lite/install/lib -Lxg -lvcflib -lgssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2 -lgcsa2 -lxg -lsdsl -ldivsufsort -ldivsufsort64
LIBS=gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o mapper.o region.o progress_bar/progress_bar.o vg_set.o utility.o path.o alignment.o edit.o sha1/sha1.o json2pb.o entropy.o

#Some little adjustments to build on OSX
#(tested with gcc4.9 and jansson installed from MacPorts)
SYS=$(shell uname -s)
ifeq (${SYS},Darwin)
	CXXFLAGS:=$(CXXFLAGS) -msse2 #needed to link against gssw
	LDFLAGS:=$(LDFLAGS) -L/opt//local/lib/ # needed for macports jansson
	STATICFLAGS= # -static doesn't work on OSX unless libgcc compiled as static. 
	ROCKSDB_PORTABLE=PORTABLE=1 # needed to build rocksdb without weird assembler options
	CLEAN_SNAPPY_AG=sed -i -e "s/[[:<:]]libtoolize[[:>:]]/glibtoolize/g" autogen.sh
else
	LDFLAGS:=$(LDFLAGS) -lrt
	STATICFLAGS=-static -static-libstdc++ -static-libgcc
	ROCKSDB_PORTABLE=
	CLEAN_SNAPPY_AG=:
endif

all: vg libvg.a

get-deps:
	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake libtool jq samtools

test: vg libvg.a test/build_graph
	cd test && $(MAKE)

test/build_graph: test/build_graph.cpp libvg.a
	$(CXX) $(CXXFLAGS) test/build_graph.cpp $(INCLUDES) -lvg $(LDFLAGS) -o test/build_graph

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

$(LIBPROTOBUF): protobuf/src/google/protobuf/*cc  protobuf/src/google/protobuf/*h
	cd protobuf && mkdir -p build && ./autogen.sh && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install
	cp protobuf/build/lib/libprotobuf.a protobuf/

$(LIBSNAPPY): snappy/*cc snappy/*h
	cd snappy && mkdir -p build && $(CLEAN_SNAPPY_AG) && ./autogen.sh && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install
	cp snappy/build/lib/libsnappy.a snappy/

$(LIBROCKSDB): rocksdb/include/rocksdb/*.h rocksdb/db/*.c rocksdb/db/*.cc rocksdb/db/*.h
	cd rocksdb && $(ROCKSDB_PORTABLE) $(MAKE) static_lib

$(LIBGCSA2): gcsa2/*.cpp gcsa2/*.h $(SDSLLITE)
	cd gcsa2 && cat Makefile | grep -v VERBOSE_STATUS_INFO >Makefile.quiet && $(MAKE) -f Makefile.quiet libgcsa2.a
	touch $(LIBGCSA2)

$(LIBXG): xg/*.cpp xg/*.hpp $(SDSLLITE)
	cd xg && $(MAKE) libxg.a

$(SDSLLITE): sdsl-lite/lib/*.cpp sdsl-lite/include/sdsl/*.hpp
	cd sdsl-lite && mkdir -p install && ./install.sh `pwd`/install
	touch $(SDSLLITE)

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

vg.o: vg.cpp vg.hpp cpp/vg.pb.h $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) $(SPARSEHASH) lru_cache/lru_cache.h stream.hpp $(LIBPROTOBUF) $(SDSLLITE)
	$(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

gssw_aligner.o: gssw_aligner.cpp gssw_aligner.hpp cpp/vg.pb.h $(LIBGSSW) $(LIBPROTOBUF) $(SPARSEHASH)
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

alignment.o: alignment.cpp alignment.hpp $(LIBHTS)  $(LIBPROTOBUF) $(SPARSEHASH)
	$(CXX) $(CXXFLAGS) -c -o alignment.o alignment.cpp $(INCLUDES)

sha1/sha1.o: sha1/sha1.cpp sha1/sha1.hpp
	$(CXX) $(CXXFLAGS) -c -o sha1/sha1.o sha1/sha1.cpp $(INCLUDES)

json2pb.o: json2pb.cpp json2pb.h bin2ascii.h $(LIBPROTOBUF)
	$(CXX) $(CXXFLAGS) -c -o json2pb.o json2pb.cpp $(INCLUDES)

entropy.o: entropy.cpp entropy.hpp
	$(CXX) $(CXXFLAGS) -c -o entropy.o entropy.cpp $(INCLUDES)

vg: $(LIBS) $(LIBVCFLIB) $(fastahack/Fasta.o) $(LIBGSSW) $(LIBROCKSDB) $(LIBSNAPPY) $(LIBHTS) $(LIBPROTOBUF) $(LIBGCSA2) $(SPARSEHASH) $(SDSLLITE) $(LIBXG)
	$(CXX) $(CXXFLAGS) -o vg $(LIBS) $(INCLUDES) $(LDFLAGS) $(STATICFLAGS)

libvg.a: vg
	ar rs libvg.a $(LIBS)

clean-vg:
	rm -f vg
	rm -f cpp/*
	rm -f *.o
	rm -f libvg.a
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
	cd gcsa2 && $(MAKE) clean
	rm -f $(SDSLLITE) && cd sdsl-lite && ./uninstall.sh `pwd`/install

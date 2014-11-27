.PHONY: all clean

CXX=g++
CXXFLAGS=-O3
pb2json=pb2json/libpb2json.a
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
LIBGSSW=gssw/src/libgssw.a
LIBSNAPPY=snappy/libsnappy.a
#LIBHYPERLEVELDB=HyperLevelDB/libhyperleveldb.a
LIBLEVELDB=leveldb/libleveldb.a
INCLUDES=-Ipb2json -Icpp -I$(VCFLIB)/src -I$(VCFLIB) -Ifastahack -Igssw/src -Ileveldb/include
LDFLAGS=-Lpb2json -Lvcflib -Lgssw/src -Lsnappy -Lleveldb -lpb2json -lvcflib -lgssw -lprotobuf -lpthread -ljansson -lleveldb -lsnappy -lz
LIBS=gssw_aligner.o vg.o cpp/vg.pb.o main.o index.o

all: vg

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

$(LIBSNAPPY): snappy/*cc snappy/*h
	cd snappy && mkdir -p build && ./autogen.sh && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install
	cp snappy/build/lib/libsnappy.a snappy/

$(LIBLEVELDB): leveldb/include/leveldb/*.h leveldb/db/*.c leveldb/db/*.cc leveldb/db/*.h
	cd leveldb && $(MAKE) libleveldb.a

#$(LIBHYPERLEVELDB):
#	cd HyperLevelDB && autoreconf -i && mkdir -p build && ./configure --prefix=`pwd`/build/ && $(MAKE) && $(MAKE) install

#test: libsnappy.a libleveldb.a
#        g++ -pthread -Ileveldb/include test.cpp -o test -L. -lleveldb -lsnappy

$(pb2json):
	cd pb2json && $(MAKE) libpb2json.a

cpp/vg.pb.cc: cpp/vg.pb.h
cpp/vg.pb.h: vg.proto
	mkdir -p cpp
	protoc vg.proto --cpp_out=cpp

test: test.cpp cpp/vg.pb.o
	$(CXX) $(CXXFLAGS) test.cpp -o test cpp/vg.pb.o $(INCLUDES) $(LDFLAGS)

$(LIBVCFLIB): vcflib/src/Variant.h vcflib/src/Variant.cpp
	cd vcflib && $(MAKE) libvcflib.a

$(LIBGSSW): gssw/src/gssw.c gssw/src/gssw.h
	cd gssw/src && $(MAKE) libgssw.a

fastahack/Fasta.o: fastahack/Fasta.h fastahack/Fasta.cpp
	cd fastahack && $(MAKE)

cpp/vg.pb.o: cpp/vg.pb.h cpp/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

vg.o: vg.cpp vg.h cpp/vg.pb.h $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW)
	$(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

gssw_aligner.o: gssw_aligner.cpp gssw_aligner.h cpp/vg.pb.h $(LIBGSSW)
	$(CXX) $(CXXFLAGS) -c -o gssw_aligner.o gssw_aligner.cpp $(INCLUDES)

main.o: main.cpp $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW)
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

index.o: index.cpp index.h
	$(CXX) $(CXXFLAGS) -c -o index.o index.cpp $(INCLUDES)

vg: $(LIBS) $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW) $(LIBLEVELDB) $(LIBSNAPPY)
	$(CXX) $(CXXFLAGS) -o vg $(LIBS) $(INCLUDES) $(LDFLAGS)

clean:
	rm -f cpp/*
	rm -f test
	rm -f vg
	cd pb2json && $(MAKE) clean
	cd vcflib && $(MAKE) clean
	cd snappy && $(MAKE) clean
	rm snappy/libsnappy.a
	cd leveldb && $(MAKE) clean

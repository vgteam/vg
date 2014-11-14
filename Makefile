.PHONY: all clean

CXX=g++
CXXFLAGS=-O3
pb2json=pb2json/libpb2json.a
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
LIBGSSW=gssw/src/libgssw.a
INCLUDES=-Ipb2json -Icpp -I$(VCFLIB)/src -I$(VCFLIB)/ -Ifastahack/ -Igssw/src/
LDFLAGS=$(pb2json) -Lpb2json -lpb2json -Lvcflib -lvcflib -Lgssw/src -lgssw -lprotobuf -lpthread -ljansson -lz
LIBS=gssw_aligner.o vg.o cpp/vg.pb.o main.o

all: vg

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

vg: $(LIBS) $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json) $(LIBGSSW)
	$(CXX) $(CXXFLAGS) -o vg $(LIBS) $(INCLUDES) $(LDFLAGS)

clean:
	rm -f cpp/*
	rm -f test
	rm -f vg
	cd pb2json && $(MAKE) clean
	cd vcflib && $(MAKE) clean

.PHONY: all clean

CXX=g++
CXXFLAGS=-O3
pb2json=pb2json/libpb2json.a
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
INCLUDES=-Ipb2json -Icpp -I$(VCFLIB)/src -I$(VCFLIB)/ -Ifastahack/
LDFLAGS=$(pb2json) -Lpb2json -lpb2json -Lvcflib -lvcflib -lprotobuf -lpthread -ljansson -lz

all: vg

$(pb2json):
	cd pb2json && $(MAKE)

cpp/vg.pb.cc: cpp/vg.pb.h
cpp/vg.pb.h: vg.proto
	mkdir -p cpp
	protoc vg.proto --cpp_out=cpp

test: test.cpp cpp/vg.pb.o
	$(CXX) $(CXXFLAGS) test.cpp -o test cpp/vg.pb.o $(INCLUDES) $(LDFLAGS)

$(LIBVCFLIB): vcflib/src/Variant.h vcflib/src/Variant.cpp
	cd vcflib && $(MAKE) libvcflib.a

fastahack/Fasta.o: fastahack/Fasta.h fastahack/Fasta.cpp
	cd fastahack && $(MAKE)

cpp/vg.pb.o: cpp/vg.pb.h cpp/vg.pb.cc
	$(CXX) $(CXXFLAGS) -c -o cpp/vg.pb.o cpp/vg.pb.cc $(INCLUDES)

vg.o: vg.cpp vg.h cpp/vg.pb.h $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json)
	$(CXX) $(CXXFLAGS) -c -o vg.o vg.cpp $(INCLUDES)

main.o: main.cpp $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json)
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp $(INCLUDES)

vg: main.o vg.o cpp/vg.pb.o $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json)
	$(CXX) $(CXXFLAGS) -o vg vg.o cpp/vg.pb.o main.o $(INCLUDES) $(LDFLAGS)

clean:
	rm -f cpp/*
	cd pb2json && $(MAKE) clean
	cd vcflib && $(MAKE) clean

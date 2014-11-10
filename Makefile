.PHONY: all clean

CXX=g++
pb2json=pb2json/libpb2json.a
VCFLIB=vcflib
LIBVCFLIB=$(VCFLIB)/libvcflib.a
INCLUDES=-Ipb2json -Icpp -I$(VCFLIB)/src -I$(VCFLIB)/ -Ifastahack/
LDFLAGS=-Lpb2json -lpb2json -Lvcflib -lvcflib -lprotobuf -lpthread -ljansson -lz

all: $(pb2json) cpp/vg.bp.h test

$(pb2json):
	cd pb2json && $(MAKE)

cpp/vg.bp.h cpp/vg.bp.cpp: vg.proto
	mkdir -p cpp
	protoc vg.proto --cpp_out=cpp

test: test.cpp cpp/vg.bp.h
	$(CXX) test.cpp -o test cpp/vg.pb.cc $(pb2json) $(INCLUDES) $(LDFLAGS)

$(LIBVCFLIB): vcflib/src/Variant.h vcflib/src/Variant.cpp
	cd vcflib && $(MAKE) libvcflib.a

fastahack/Fasta.o: fastahack/Fasta.h fastahack/Fasta.cpp
	cd fastahack && $(MAKE)

vg: vg.cpp vg.h cpp/vg.bp.h $(LIBVCFLIB) $(fastahack/Fasta.o) $(pb2json)
	$(CXX) -o vg main.cpp vg.cpp cpp/vg.pb.cc $(INCLUDES) $(LDFLAGS)

clean:
	rm -f cpp/*
	cd pb2json && $(MAKE) clean
	cd vcflib && $(MAKE) clean

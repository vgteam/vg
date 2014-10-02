.PHONY: all clean

CXX=g++
pb2json=pb2json/libpb2json.a
INCLUDES=-Ipb2json -Icpp
LDFLAGS=-Lpb2json -lpb2json -lprotobuf -lpthread -ljansson

all: $(pb2json) cpp/vg.bp.h test

$(pb2json):
	cd pb2json && $(MAKE)

cpp/vg.bp.h cpp/vg.bp.cpp: vg.proto
	mkdir -p cpp
	protoc vg.proto --cpp_out=cpp

test: test.cpp cpp/vg.bp.h
	$(CXX) test.cpp -o test cpp/vg.pb.cc $(pb2json) $(INCLUDES) $(LDFLAGS)

clean:
	rm -f cpp/*
	cd pb2json && $(MAKE) clean

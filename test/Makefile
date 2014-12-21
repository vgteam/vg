.PHONY: all clean

CXX=g++ -std=c++11 -fopenmp
CXXFLAGS=-O3
VCFLIB=vcflib
INCLUDES=-I../ -I../pb2json -I../cpp -I../$(VCFLIB)/src -I../$(VCFLIB) -I../fastahack -I../gssw/src -I../rocksdb/include
LDFLAGS=-L../ -L../pb2json -L../vcflib -L../gssw/src -L../snappy -L../rocksdb -lvg -lpb2json -lvcflib -lgssw -lprotobuf -lpthread -ljansson -lrocksdb -lsnappy -lz -lbz2
vg=../vg

all: test clean

test: build_graph $(vg)
	prove -v t

$(vg):
	cd .. && make vg

build_graph: build_graph.cpp
	$(CXX) build_graph.cpp $(INCLUDES) $(LDFLAGS) -o build_graph

clean:
	rm -f build_graph

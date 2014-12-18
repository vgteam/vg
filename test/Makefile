.PHONY: all clean

CXX=g++
VCFLIB=vcflib
INCLUDES=-I../ -I../pb2json -I../cpp -I../$(VCFLIB)/src -I../$(VCFLIB) -I../fastahack -I../gssw/src -I../leveldb/include
LDFLAGS=-L../ -L../pb2json -L../vcflib -L../gssw/src -L../snappy -L../leveldb -lvg -lpb2json -lvcflib -lgssw -lprotobuf -lpthread -ljansson -lleveldb -lsnappy -lz
vg=../vg

all: test clean

test: make_graph $(vg)
	prove -v t

$(vg):
	cd .. && make vg

make_graph: make_graph.cpp
	$(CXX) make_graph.cpp $(INCLUDES) $(LDFLAGS) -o make_graph

clean:
	rm -f make_graph

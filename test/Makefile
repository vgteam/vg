.PHONY: all clean

CXX:=g++
CXXFLAGS:=-O3 -std=c++11 -fopenmp -g

vg:=../bin/vg

all: test clean

test: build_graph $(vg) fastahack
	prove -v t

fastahack:
	cd ../deps/vcflib/contrib/fastahack && $(MAKE)

$(vg):
	cd .. && $(MAKE) bin/vg

build_graph: build_graph.cpp
	cd .. && $(MAKE) test/build_graph

clean:
	rm -f build_graph

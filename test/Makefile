.PHONY: all clean

CXX=g++
CXXFLAGS=-O3 -std=c++11 -fopenmp -g

vg=../vg

all: test clean

test: build_graph $(vg)
	prove -v t

$(vg):
	cd .. && $(MAKE) vg

build_graph: build_graph.cpp
	cd .. && $(MAKE) test/build_graph

clean:
	rm -f build_graph

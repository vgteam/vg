.PHONY: all clean

CXX:=g++
CXXFLAGS:=-O3 -std=c++11 -fopenmp -g

vg:=../bin/vg

all: test clean

test: build_graph $(vg)
	prove -v t

$(vg):
	cd .. && $(MAKE) bin/vg

build_graph: build_graph.cpp
	cd .. && . ./source_me.sh && $(MAKE) test/build_graph

clean:
	rm -f build_graph

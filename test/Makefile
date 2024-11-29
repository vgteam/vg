.PHONY: all clean

CXX:=g++
CXXFLAGS:=-O3 -std=c++11 -fopenmp -g

vg:=../bin/vg

all: test clean

test: build_graph $(vg) vcf2tsv fastahack
	prove -v t

vcf2tsv:
	cd ../deps/vcflib && $(MAKE) bin/vcf2tsv

fastahack:
	cd ../deps/fastahack && $(MAKE)

$(vg):
	cd .. && $(MAKE) bin/vg

build_graph: build_graph.cpp
	cd .. && $(MAKE) test/build_graph

clean:
	rm -f build_graph

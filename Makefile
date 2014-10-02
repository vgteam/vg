.PHONY: all clean

all: cpp/vg.bp.h

cpp/vg.bp.h cpp/vg.bp.cpp: vg.proto
	mkdir -p cpp
	protoc vg.proto --cpp_out=cpp

clean:
	rm -f cpp/*

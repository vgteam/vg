# vg

[![Build Status](https://travis-ci.org/ekg/multinomial-ln.svg)](https://travis-ci.org/ekg/multinomial-ln)

## variant graph data structures, interchange formats, alignment, genotyping, and variant calling methods

### building

You'll need the protobuf and jansson development libraries installed on your server.

    sudo apt-get install libprotoc-dev libjansson-dev

Now, obtain the repo and its submodules:

    git clone --recursive https://github.com/ekg/vg.git

Then build with `make`, and run with `./vg`.

### Done

- data models for reference graph and alignments against it (vg.proto)
- local alignment against the graph (vg.cpp)
- index capable of storing large graphs on disk and efficiently retrieving subgraphs (index.cpp)
- protobuf, json, and dot format serialization (view)
- command-line interfaces: construct, view, index, find, align, paths (main.cpp)
- tap-compliant tests

### To do

- GFA input and output (efficient use requires bluntifying the graph, removing node-node overlaps)
- positional indexing (can be done on graph constructed from VCF+fasta reference)
- interface harmonization of in-memory (vg.cpp) and on-disk (index.cpp) graph representations
- per-node, per-sample quality and count information
- dynamic programming method to estimate path qualities given per-node qualities and counts
- genotype likelihood generation (given a source and sink, genotype paths)
- genotyping using freebayes genotyping model
- generalization to assembly graphs (although directional, nothing is intrinsically DAG-based except alignment)

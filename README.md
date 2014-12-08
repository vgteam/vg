# vg

[![Build Status](https://travis-ci.org/ekg/vg.svg)](https://travis-ci.org/ekg/vg)

## variant graph data structures, interchange formats, alignment, genotyping, and variant calling methods

If we know about variation in a given population, we should include that knowledge in our primary sequence analyses, or risk bias against things we've seen before. Reference bias is real. We can work around it by formulating our reference system as a graph: either an assembly, or a directed acyclic one similar to how we represent a multiple sequence alignment.

## Usage

### building

You'll need the protobuf and jansson development libraries installed on your server.

    sudo apt-get install libprotoc-dev libjansson-dev

Now, obtain the repo and its submodules:

    git clone --recursive https://github.com/ekg/vg.git

Then build with `make`, and run with `./vg`.

### What can I do?

Try building a graph and aligning to it:

    vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
    vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg

Note that you don't have to store the graph on disk at all, you can simply pipe it into the local aligner:

    vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -

A variety of commands are available:

- *construct*: graph construction
- *view*: conversion (protobuf/json/GFA)
- *index*: index features of the graph in a disk-backed key/value store
- *find*: use an index to find nodes, edges, kmers, or positions
- *paths*: traverse paths in the graph
- *align*: local alignment
- *stats*: metrics describing graph properties
- *join*: combine graphs

## Development

### Done

- data models for reference graph and alignments against it (vg.proto)
- local alignment against the graph (vg.cpp)
- index capable of storing large graphs on disk and efficiently retrieving subgraphs (index.cpp)
- protobuf, json, and dot format serialization (view)
- binary format for graph and alignments against it
- command-line interfaces: construct, view, index, find, align, paths (main.cpp)
- tap-compliant tests
- kmer-based indexing of the graph
- graph statistics
- subgraph decomposition
- k-path enumeration
- graph joining: combine subgraphs represented in a single or different .vg files

### To do

- global mapping against large graphs (depends on indexing of graph and kmers, which are done)
- alignment streams (via protobuf's ZeroCopyInputStream/ZeroCopyOutputStream interface)
- GFA input and output (efficient use requires bluntifying the graph, removing node-node overlaps)
- positional indexing (can be done on graph constructed from VCF+fasta reference)
- interface harmonization of in-memory (vg.cpp) and on-disk (index.cpp) graph representations
- per-node, per-sample quality and count information
- dynamic programming method to estimate path qualities given per-node qualities and counts
- genotype likelihood generation (given a source and sink, genotype paths)
- genotyping of paths using freebayes-like genotyping model
- generalization to assembly graphs (although directional, nothing is intrinsically DAG-based except alignment)

## License

MIT

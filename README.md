# vg

[![Build Status](https://travis-ci.org/ekg/vg.svg)](https://travis-ci.org/ekg/vg)

## variant graph data structures, interchange formats, alignment, genotyping, and variant calling methods

If we know about variation in a given population, we should include that knowledge in our primary sequence analyses, or risk bias against things we've seen before. Reference bias is real. We can work around it by formulating our reference system as a graph: either an assembly, or a directed acyclic one similar to how we represent a multiple sequence alignment.

![Variation graph](https://raw.githubusercontent.com/ekg/vg/master/figures/example.png)

## Usage

### building

You'll need the protobuf and jansson development libraries installed on your server.

    sudo apt-get install protobuf-compiler libprotoc-dev libjansson-dev

You can also run `make get-deps`.

Other libraries may be required. Please report any build difficulties.

Now, obtain the repo and its submodules:

    git clone --recursive https://github.com/ekg/vg.git

Then build with `make`, and run with `./vg`.

### What can I do?

Try building a graph and aligning to it:

```sh
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg
```

Note that you don't have to store the graph on disk at all, you can simply pipe it into the local aligner:

```sh
vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -
```

You can also index and then map reads against the index of the graph:

```sh
# construct the graph
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg

# store the graph in the index, and also index the kmers in the graph of size 11
vg index -s -k 11 x.vg

# align a read to the indexed version of the graph
# note that the graph file is not opened, but x.vg.index is assumed
vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -k 11 x.vg >alignment.json
```

A variety of commands are available:

- *construct*: graph construction
- *view*: conversion (dot/protobuf/json/GFA)
- *index*: index features of the graph in a disk-backed key/value store
- *find*: use an index to find nodes, edges, kmers, or positions
- *paths*: traverse paths in the graph
- *align*: local alignment
- *map*: global alignment (kmer-driven)
- *stats*: metrics describing graph properties
- *join*: combine graphs (parallel)
- *concat*: combine graphs (serial)
- *ids*: id manipulation

## Implementation notes

`vg` is based around a graph object (vg::VG) which has a native serialized representation that is almost identical on disk and in-memory, with the exception of adjacency indexes that are built when the object is parsed from a stream or file. These graph objects are the results of queries of larger indexes, or manipulation (for example joins or concatenations) of other graphs. vg is designed for interactive, stream-oriented use. You can, for instance, construct a graph, merge it with another one, and pipe the result into a local alignment process. The graph object can be stored in an index (vg::Index), aligned against directly (vg::GSSWAligner), or "mapped" against in a global sense (vg::Mapper), using an index of kmers.

Once constructed, a variant graph (.vg is the suggested file extension) is typically around the same size as the reference (FASTA) and uncompressed variant set (VCF) which were used to build it. The index, however, may be much larger, perhaps more than an order of magnitude. This is less of a concern as it is not loaded into memory, but could be a pain point as vg is scaled up to whole-genome mapping.

The serialization of very large graphs (>62MB) is enabled by the use of protocol buffer ZeroCopyStreams. Graphs are decomposed into sets of N (presently 10k) nodes, and these are written, with their edges, into graph objects that can be streamed into and out of vg. Graphs of unbounded size are possible using this approach.

## Development

- [x] data models for reference graph and alignments against it (vg.proto)
- [x] local alignment against the graph (vg.cpp)
- [x] index capable of storing large graphs on disk and efficiently retrieving subgraphs (index.cpp)
- [x] protobuf, json, and dot format serialization (view)
- [x] binary format for graph and alignments against it
- [x] command-line interfaces: construct, view, index, find, align, paths (main.cpp)
- [x] tap-compliant tests
- [x] kmer-based indexing of the graph
- [x] graph statistics
- [x] subgraph decomposition
- [x] k-path enumeration
- [x] limiting k-paths to only those crossing a certain number of nodes (this prevents kpath blowup in densely varying regions of the graph)
- [x] graph joining: combine subgraphs represented in a single or different .vg files
- [x] GFA output
- [x] global mapping against large graphs
- [x] non-recursive topological sort of graph
- [x] stable ID compaction
- [x] efficient construction for large DAGs
- [x] improve memory performance of kmer indexing for large graphs by storing incremental results of k-path generation
- [x] global alignment: retain and expand only the most-likely subgraphs
- [x] verify that snappy compression is enabled for index, and measure size for large graphs
- [x] move to rocksdb for better indexing performance on modern hardware (multiple cores, SSDs)
- [x] object streams (enable graphs > 60mb) and alignment streams (via protobuf's ZeroCopyInputStream/ZeroCopyOutputStream interface)
- [x] use dense_hash for improved memory and runtime efficiency with large graphs (or sparse_hash, if memory is at a premium--- but it's easy to switch and ideally we can design large-scale construction without loading entire whole-genome graphs into memory)
- [x] GFA input (efficient use requires bluntifying the graph, removing node-node overlaps), and probably default GFA output from vg view
- [x] index metadata (to quickly check if we have kmer index of size >=N)
- [x] use divide-and-conquer for graph fragment concatenation during construction
- [x] simplify mapping by setting a maximum node size in construction
- [x] kmer falloff in global alignment (if we can't find hits at a kmer size of K, try K-n; enabled by the sorted nature of the index's key-value backend)
- [x] positional indexing for improved global mapping (can be done on graph constructed from VCF+fasta reference)
- [x] index the kmers of large graphs in reasonable time (48 hours, 32 threads, 2500 samples in 1000 genomes phase 3)
- [x] compression of serialization format
- [ ] interface harmonization of in-memory (vg.cpp) and on-disk (index.cpp) graph representations
- [x] emded paths in serialized graph format (important especially in the case of reference paths)
- [ ] alignment serialization format
- [x] prune non-informative kmers from index, so as to save space
- [ ] per-node, per-sample quality and count information on graph
- [ ] should an alignment be a graph too? : express a sample's sequencing results as a labeled graph
- [x] multiple samples in one graph (colors) - solved by paths
- [ ] modify a graph using an alignment's path, adding new nodes as needed and updating the path/alignment to match
- [ ] dynamic programming method to estimate path qualities given per-node qualities and counts
- [ ] genotype likelihood generation (given a source and sink, genotype paths)
- [ ] genotyping of paths using freebayes-like genotyping model
- [ ] genotyping using dynamic programming genotyping model and compressed sequence results against graph
- [ ] generalization to assembly graphs (although directional, nothing is intrinsically DAG-based except alignment)
- [ ] [perceptual DNA hashing](http://arxiv-web3.library.cornell.edu/abs/1412.5517) to reduce index memory usage

## License

MIT

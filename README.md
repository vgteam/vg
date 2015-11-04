[![Stories in Ready](https://badge.waffle.io/ekg/vg.png?label=ready&title=Ready)](https://waffle.io/ekg/vg)
# vg

[![Join the chat at https://gitter.im/ekg/vg](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ekg/vg?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/ekg/vg.svg)](https://travis-ci.org/ekg/vg)

## variation graph data structures, interchange formats, alignment, genotyping, and variant calling methods

![Variation graph](https://raw.githubusercontent.com/ekg/vg/master/figures/example.png)

_Variation graphs_ provide a succinct encoding of the sequences of many genomes. A variation graph (in particular as implemented in vg) is composed of:

* _nodes_, which are labeled by sequences and ids
* _edges_, which connect two nodes via either of their respective ends
* _paths_, describe genomes, sequence alignments, and annotations (such as gene models and transcripts) as walks through nodes connected by edges

This model is similar to a number of sequence graphs that have been used in assembly and multiple sequence alignment. Paths provide coordinate systems relative to genomes encoded in the graph, allowing stable mappings to be produced even if the structure of the graph is changed. For visual documentation, please refer to a presentation on the topic: [Resequencing against a human whole genome variation graph](https://docs.google.com/presentation/d/1bbl2zY4qWQ0yYBHhoVuXb79HdgajRotIUa_VEn3kTpI/edit?usp=sharing) (April 14, 2015).

## Usage

### building

Before you begin, you'll need to install some basic tools if they are not already installed.

    sudo apt-get install git cmake pkg-config libncurses-dev libbz2-dev

You'll need the protobuf and jansson development libraries installed on your server.

    sudo apt-get install protobuf-compiler libprotoc-dev libjansson-dev automake libtool
    
Additionally, to run the tests, you will need jq.

    sudo apt-get install jq

You can also run `make get-deps`.

Other libraries may be required. Please report any build difficulties.

Now, obtain the repo and its submodules:

    git clone --recursive https://github.com/ekg/vg.git

Then build with `make`, and run with `./vg`.

#### building on Mac OS X

VG won't build with XCode's compiler (clang), but it should work with GCC 4.9.  One way to install the latter (and other dependencies) is to install [Mac Ports](https://www.macports.org/install.php), then run:

    sudo port install gcc49 libtool jansson jq

To make GCC 4.9 the default compiler, run (use `none` instead of `mp-gcc49` to revert back):

    sudo port select gcc mp-gcc49

VG can now be cloned and built as described above.

### Variation graph construction

The simplest thing to do with `vg` is to build a graph and align to it. At present, you'll want to use a reference and VCF file to do so. If you're working in the `test/` directory:

```sh
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
```

### Viewing, conversion

`vg view` provides a way to convert the graph into various formats:

```sh
# GFA output
vg view x.vg >x.gfa

# dot output suitable for graphviz
vg view -d x.vg >x.dot

# json version of binary alignments
vg view -a x.gam >x.json
```

### Alignment

As this is a small graph, you could align to it using a full-length partial order alignment:

```sh
vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg
```

Note that you don't have to store the graph on disk at all, you can simply pipe it into the local aligner:

```sh
vg construct -r small/x.fa -v small/x.vcf.gz | vg align -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -
```

Most commands allow the streaming of graphs into and out of `vg`.

### Mapping

If your graph is large, you want to use `vg index` to store the graph and `vg map` to align reads. `vg map` implements a kmer based seed and extend alignment model that is similar to that used in aligners like novoalign or MOSAIK. First an on-disk index is built with `vg index` which includes the graph itself and kmers of a particular size. When mapping, any kmer size shorter than that used in the index can be employed, and by default the mapper will decrease the kmer size to increase sensitivity when alignment at a particular _k_ fails.

```sh
# construct the graph
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg

# store the graph in the index, and also index the kmers in the graph of size 11
# you can provide a list of .vg files on the command line, which is useful if you
# have constructed a graph for each chromosome in a large reference
vg index -s -k 11 x.vg

# align a read to the indexed version of the graph
# note that the graph file is not opened, but x.vg.index is assumed
vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG x.vg >read.gam

# simulate a bunch of 150bp reads from the graph and map them
vg map -r <(vg sim -n 1000 -l 150 x.vg) x.vg >aln.gam

# surject the alignments back into the reference space of sequence "x", yielding a BAM file
vg surject -p x -b aln.gam >aln.bam
```

### Command line interface

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
- *kmers*: generate kmers from a graph
- *sim*: simulate reads by walking paths in the graph
- *mod*: various transformations of the graph
- *surject*: force graph alignments into a linear reference space

## Implementation notes

`vg` is based around a graph object (vg::VG) which has a native serialized representation that is almost identical on disk and in-memory, with the exception of adjacency indexes that are built when the object is parsed from a stream or file. These graph objects are the results of queries of larger indexes, or manipulation (for example joins or concatenations) of other graphs. vg is designed for interactive, stream-oriented use. You can, for instance, construct a graph, merge it with another one, and pipe the result into a local alignment process. The graph object can be stored in an index (vg::Index), aligned against directly (vg::GSSWAligner), or "mapped" against in a global sense (vg::Mapper), using an index of kmers.

Once constructed, a variation graph (.vg is the suggested file extension) is typically around the same size as the reference (FASTA) and uncompressed variant set (VCF) which were used to build it. The index, however, may be much larger, perhaps more than an order of magnitude. This is less of a concern as it is not loaded into memory, but could be a pain point as vg is scaled up to whole-genome mapping.

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
- [x] interface harmonization of in-memory (vg.cpp) and on-disk (index.cpp) graph representations
- [x] emded paths in serialized graph format (important especially in the case of reference paths)
- [x] alignment serialization format
- [x] prune non-informative kmers from index, so as to save space
- [x] per-node, per-sample quality and count information on graph
- [x] should an alignment be a graph too? : express a sample's sequencing results as a labeled graph
- [x] multiple samples in one graph (colors) - solved by paths
- [x] modify a graph using an alignment's path, adding new nodes as needed and updating the path/alignment to match
- [x] modify a graph using many alignments (update mod procedure to include many mappings)
- [x] use an alignment to add a named path to the graph
- [x] path range query from index (give me the subgraph corresponding to a particular genome location)
- [x] _surject_ alignments back into arbitrary path (such as GRCh37)
- [x] use htslib for BAM/CRAM/SAM i/o
- [ ] factor longer functions in main.cpp into library functions for easier reuse
- [x] sort alignments and output a sorted GAM stream
- [x] kmers input to GCSA: kmer, starting position = (node id, offset), previous characters, successive characters, successive positions
- [x] properly handle pairs in alignment and surjection
- [x] read paired-end data in FASTQ format
- [x] allow the efficient alignment of very long sequences by banding
- [ ] dynamic programming method to estimate path qualities given per-node qualities and counts
- [ ] genotype likelihood generation (given a source and sink, genotype paths)
- [ ] genotyping of paths using freebayes-like genotyping model
- [ ] genotyping using dynamic programming genotyping model and compressed sequence results against graph
- [x] generalization to assembly graphs (although directional, nothing is intrinsically DAG-based except alignment)
- [ ] [perceptual DNA hashing](http://arxiv-web3.library.cornell.edu/abs/1412.5517) to reduce index memory usage
- [ ] implement pBWT for quick haplotype lookup, which is useful to constrain kmer indexing space and also to support long haplotype driven genotyping
- [x] compressed bitvectors and intvectors from sdsl-lite for storage of genome paths and annotations

## License

MIT

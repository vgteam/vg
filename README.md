# vg

[![Join the chat at https://gitter.im/vgteam/vg](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/vgteam/vg?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![Build Status](https://travis-ci.org/vgteam/vg.svg?branch=master)](https://travis-ci.org/vgteam/vg) [![Performance Report](https://img.shields.io/badge/performance-report-brightgreen.svg)](http://cgl-pipeline-inputs.s3.amazonaws.com/vg_cgl/vg_ci/jenkins_reports/branch/master/index.html) [![Stories in Ready](https://badge.waffle.io/vgteam/vg.png?label=ready&title=Ready)](https://waffle.io/vgteam/vg)
[![Documentation Status](https://readthedocs.org/projects/vg/badge/?version=latest)](http://vg.readthedocs.io/en/latest/?badge=latest)

## variation graph data structures, interchange formats, alignment, genotyping, and variant calling methods

![Variation graph](https://raw.githubusercontent.com/vgteam/vg/master/doc/figures/vg_logo.png)

_Variation graphs_ provide a succinct encoding of the sequences of many genomes. A variation graph (in particular as implemented in vg) is composed of:

* _nodes_, which are labeled by sequences and ids
* _edges_, which connect two nodes via either of their respective ends
* _paths_, describe genomes, sequence alignments, and annotations (such as gene models and transcripts) as walks through nodes connected by edges

This model is similar to a number of sequence graphs that have been used in assembly and multiple sequence alignment. Paths provide coordinate systems relative to genomes encoded in the graph, allowing stable mappings to be produced even if the structure of the graph is changed.

![example variation graph](https://raw.githubusercontent.com/vgteam/vg/master/doc/figures/smallgraph.png)

## Usage

### building
e
Before you begin, you'll need to install some basic tools if they are not already installed. You'll need the protobuf and jansson development libraries installed on your server. Additionally, to run the tests, you will need `jq`, `bc` and `rs`.

    sudo apt-get install build-essential git cmake pkg-config libncurses-dev libbz2-dev  \
                         protobuf-compiler libprotoc-dev libjansson-dev automake libtool \
                         jq bc rs curl unzip redland-utils librdf-dev bison flex gawk \
                         lzma-dev liblzma-dev liblz4-dev  

You can also run `make get-deps`.

At present, you will need GCC version 4.9 or greater to compile vg. (Check your version with `gcc --version`.)

Other libraries may be required. Please report any build difficulties.

Note that a 64-bit OS is required. Ubuntu 16.04 should work.

Now, obtain the repo and its submodules:

    git clone --recursive https://github.com/vgteam/vg.git

Then build with `. ./source_me.sh && make static`, and run with `./bin/vg`.

#### building on Mac OS X

##### using Mac Ports

VG won't build with XCode's compiler (clang), but it should work with GCC >= 4.9.  One way to install the latter (and other dependencies) is to install [Mac Ports](https://www.macports.org/install.php), then run:

    sudo port install gcc7 libtool jansson jq cmake pkgconfig autoconf automake libtool coreutils samtools redland bison gperftools md5sha1sum rasqal gmake autogen cairo clang-3.8

To make GCC 7 the default compiler, run (use `none` instead of `mp-gcc7` to revert back):

    sudo port select gcc mp-gcc7

Some OSX users also need to have the MacPorts clang assembler for dependencies (use `none` instead of `mp-clang-3.8` to revert back):

    sudo port select clang mp-clang-3.8

VG can now be cloned and built:

    git clone --recursive https://github.com/vgteam/vg.git
    cd vg
    . ./source_me.sh && make
    
**Note that static binaries cannot yet be built for Mac.**

Our team has also successfully built vg on Mac with GCC versions 4.9, 5.3, 6, and 7.3.

##### using Homebrew

[Homebrew](http://brew.sh/) provides another package management solution for OSX, and may be preferable to some users over MacPorts.

```
brew install automake libtool jq jansson coreutils gcc49 samtools pkg-config cmake raptor bison lz4 xz
export PATH="/usr/local/opt/coreutils/libexec/gnubin:/usr/local/bin:$PATH"

# Force use of new version of bison
brew link bison --force

# Use glibtool/ize
export LIBTOOL=glibtool
export LIBTOOLIZE=glibtoolize
# Make symlinks to use gxx-4.9 instead of builtin gxx (CC and CXX not yet fully honored)
ln -s gcc-4.9 /usr/local/bin/gcc
ln -s g++-4.9 /usr/local/bin/g++

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH;
export LIBRARY_PATH=$LD_LIBRARY_PATH;

git clone --recursive https://github.com/vgteam/vg.git  
cd vg/  
. ./source_me.sh && make
```

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

# store the graph in the xg/gcsa index pair
vg index -x x.xg -g x.gcsa -k 16 x.vg

# align a read to the indexed version of the graph
# note that the graph file is not opened, but x.vg.index is assumed
vg map -s CTACTGACAGCAGAAGTTTGCTGTGAAGATTAAATTAGGTGATGCTTG -x x.xg -g x.gcsa >read.gam

# simulate a bunch of 150bp reads from the graph and map them
vg map -r <(vg sim -n 1000 -l 150 -x x.xg ) -x x.xg -g x.gcsa >aln.gam

# surject the alignments back into the reference space of sequence "x", yielding a BAM file
vg surject -x x.xg -b aln.gam >aln.bam

# or alternatively, surject them to BAM in the call to map
vg map -r <(vg sim -n 1000 -l 150 -x x.xg ) -x x.xg -g x.gcsa --surject-to bam >aln.bam
```
### Variant Calling

The following example shows how to construct a VCF file from a read alignment and graph.  Input must be split into chunks (see vg chunk) in order to run on whole genome.

```sh
# filter secondary and ambiguous read mappings out of the gam
vg filter alignment.gam -r 0.90 -fu -s 2 -o 0 -D 999 -x graph.xg > filtered.gam

# create an augmented graph by adding variation from the reads
vg augment graph.vg filtered.gam -q 10 -S aug_graph.support -Z aug_graph.trans -A aug_alignment.gam > aug_graph.vg

# Make calls by thresholding based on read support for graph path SEQ
vg call aug_graph.vg -b graph.vg -s aug_graph.support -z aug_graph.trans -r SEQ > calls.vcf

# Or Make calls using a Freebayes-like genotyping algorithm for graph path SEQ
vg genotype aug_graph.vg -G aug_alignment.gam -E -v -r SEQ > calls.vcf

# for comparison purposes, it's very useful to normalize the vcf output, especially for more complex graphs which can make large variant blocks that contain a lot of reference bases (Note: requires [vt](http://genome.sph.umich.edu/wiki/Vt)):
vt decompose_blocksub -a calls.vcf | vt normalize -r FASTA_FILE - > calls.clean.vcf

```

To produce a VCF file for a whole chromosome, the graph must be cut up along the reference genome and called in chunks.  `scripts/chunked_call` wraps this functionality to produce chromosome-sized VCFs in a single command line (from a GAM file and XG index)

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
- *msga*: construct a graph from an assembly of multiple sequences
- *validate*: determine if graph is valid
- *filter*: filter reads out of an alignment
- *augment*: adds variation from aligned reads into the graph
- *call/genotype*: call variants from an augmented graph

## Implementation notes

`vg` is a collection of tools based on a common data model (the variation graph) that is described by a protobuf schema (vg.proto). Data objects defined in vg.proto may be serialized via a stream pattern defined in stream.hpp. It is not necessary to write code in vg in order to interface with the algorithms defined here. Rather, it is sometimes simpler to write an external algorithm that reads and writes the same data formats.

## License

MIT

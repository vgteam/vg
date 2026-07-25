# file-info

This file contains extra text that will be added to the man pages generated with doc/vgmanmd.py
The `# description` section is added to the top of the page, and each `# subcommand` section will be added to given subcommand
When adding a new subcommand, add it to the appropriate section(s) in the description

# description

[vg](https://github.com/vgteam/vg) is a toolkit for variation graph data structures, interchange formats, alignment, genotyping, and variant calling methods.

For more in-depth explanations of tools and workflows, see the [general wiki page](https://github.com/vgteam/vg/wiki).

# synopsis

- **Graph construction and indexing**
    See the [wiki page](https://github.com/vgteam/vg/wiki/File-Types) for an overview of file types.
    - [`vg autoindex`](#autoindex): automatically construct a graph and indexes for a specific workflow (e.g. giraffe, rpvg). [wiki page](https://github.com/vgteam/vg/wiki/Automatic-indexing-for-read-mapping-and-downstream-inference)
    - [`vg construct`](#construct): manually construct a graph from a reference and variants. [wiki page](https://github.com/vgteam/vg/wiki/Construction)
    - [`vg gbwt`](#gbwt): manually build and manipulate GBWTs and indexes (GBWTgraph, GBZ, r-index). [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg haplotypes`](#haplotypes): haplotype sample a graph. Recommended for mapping with giraffe. [wiki page](https://github.com/vgteam/vg/wiki/Haplotype-Sampling)
    - [`vg index`](#index): manually build individual indexes (xg, distance, GCSA, etc). [wiki page](https://github.com/vgteam/vg/wiki/Index-Construction) 
    - [`vg minimizer`](#minimizer): manually build a minimizer index for mapping. 
    - [`vg rna`](#rna): construct splicing graphs and pantranscriptomes. [wiki page](https://github.com/vgteam/vg/wiki/Transcriptomic-analyses). Also see [rpvg](https://github.com/jonassibbesen/rpvg) 
- **Read mapping**
    - [`vg giraffe`](#giraffe): fast haplotype-aware short read alignment. [wiki page](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)
    - [`vg map`](#map): MEM-based read alignment.
    - [`vg mpmap`](#mpmap): splice-aware multipath alignment of short reads. [wiki page](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-vg-mpmap)
- **Downstream analyses**
    - [`vg call`](#call): call or genotype VCF variants. Uses [vg pack](#pack). [wiki page](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling)
    - [`vg deconstruct`](#deconstruct): create a VCF from variation in the graph. [wiki page](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)
    - [`vg depth`](#depth): calculate path coverage depth or aggregate alignment coverage Also see [vg pack](#pack)
    - [`vg genotype`](#genotype): alternative, less-tested, FreeBayes-like genotyper.
    - [`vg mcmc`](#mcmc): find haplotypes using MCMC methods on reads.
    - [`vg pack`](#pack): convert alignments to a compact coverage index. Used with [vg call](#call)
    - [`vg primers`](#primers): filter primers based on variation in a graph. [wiki page](https://github.com/vgteam/vg/wiki/Primer-Filter)
    - [`vg viz`](#viz): visualize a graph. [wiki page](https://github.com/vgteam/vg/wiki/Complex-graph-visualization#using-vg-viz)
- **Working with read alignments**
    - [`vg augment`](#augment): embed alignments into a graph.
    - [`vg filter`](#filter): filter alignments by properties.
    - [`vg gamcompare`](#gamcompare): compare GAM-format alignments to a truth.
    - [`vg gampcompare`](#gampcompare): compare [GAMP-format alignments](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-the-GAMP-format) to a truth.
    - [`vg gamsort`](#gamsort): sort a GAM/GAF file or index a sorted GAM file.
    - [`vg sim`](#sim): simulate reads from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Simulating-reads-with-vg-sim)
- **Graph and read statistics**
    - [`vg chains`](#chains): get top-level chains from graph. [wiki page](https://github.com/vgteam/vg/wiki/Snarls-and-chains)
    - [`vg dotplot`](#dotplot): get dotplot matrix for paths in graph.
    - [`vg filter`](#filter): get stats about alignments (use `--tsv-out`).
    - [`vg gbwt`](#gbwt): get stats about a GBWT. [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg paths`](#paths): get stats about the paths. [wiki page](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)
    - [`vg snarls`](#snarls): get a list of snarls within the graph. [wiki page](https://github.com/vgteam/vg/wiki/Snarls-and-chains)
    - [`vg stats`](#stats): get stats about the graph.
- **Manipulating a graph**
    - [`vg annotate`](#annotate): annotate a graph or alignments.
    - [`vg circularize`](#circularize): connect head and tail nodes to circularize paths.
    - [`vg clip`](#clip): remove variation from a graph.
    - [`vg combine`](#combine): merge graphs into a combined graph. Useful for per-chromosome graphs.
    - [`vg mask`](#mask): N-mask regions of a graph.
    - [`vg mod`](#mod): filter, transform, and edit the graph.
    - [`vg gbwt`](#gbwt): manipulate GBWTs and associated indexes. [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg ids`](#ids): manipulate graph node ids.
    - [`vg paths`](#paths): manipulate paths in a graph. [wiki page](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)
    - [`vg prune`](#prune): prune the graph for GCSA2 indexing.
    - [`vg simplify`](#simplify): simplify graph by removing small leaf sites or rare variants.
- **Conversion between formats**
    - [`vg convert`](#convert): convert between handle graph formats and GFA, and between alignment formats.
    - [`vg inject`](#inject): project alignments on a linear reference onto a graph (bam/sam/cram->gam/gaf).
    - [`vg paths`](#paths): extract a fasta from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Linear-references-and-vg)
    - [`vg surject`](#surject): project alignments on a graph onto a linear reference (gam/gaf->bam/sam/cram).
    - [`vg vectorize`](#vectorize): export GAM alignments to other vector formats.
    - [`vg view`](#view): convert between non-handle graph formats and alignment formats (dot, json, turtle...).
- **Subgraph extraction**
    - [`vg chunk`](#chunk): split a graph and/or alignment into smaller chunks.
    - [`vg find`](#find): use an index to find nodes, edges, kmers, paths, or positions.
    - [`vg paths`](#paths): subset a graph to keep or remove particular paths.
    - [`vg trace`](#trace): extract (sub-)haplotypes from a graph.
- **Extremely specific analyses** (i.e. don't use unless you know what you're doing)
    - [`vg align`](#align): run highly configurable alignment of a single sequence.
    - [`vg chain`](#chain): run an anchor chaining problem as specified in a JSON.
    - [`vg cluster`](#cluster): run only the seed-clustering step of an alignment.
    - [`vg zipcode`](#zipcode): test zipcode functionality on given GAM reads.
- **Developer tools**
    - [`vg benchmark`](#benchmark): run basic graph benchmarks.
    - [`vg bench-dist-query`](#bench-dist-query): run distance indexing benchmarks.
    - [`vg describe`](#describe): determine file format of input. [wiki page](https://github.com/vgteam/vg/wiki/File-Types)
    - [`vg help`](#help): display all subcommands.
    - [`vg test`](#test): run unit tests.
    - [`vg validate`](#validate): check if a graph is well-formed.
    - [`vg version`](#version): display current version number.

# bugs

Bugs can be reported at: https://github.com/vgteam/vg/issues

For technical support, please visit: https://www.biostars.org/tag/vg/


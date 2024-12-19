# file-info

This file contains extra text that will be added to the man pages generated with doc/vgmanmd.py
The `# description` section is added to the top of the page, and each `# subcommand` section will be added to given subcommand
When adding a new subcommand, add it to the appropriate section(s) in the description

# description

vg is a toolkit for variation graph data structures, interchange formats, alignment, genotyping, and variant calling methods.

For more in-depth explanations of tools and workflows, see the [vg wiki page](https://github.com/vgteam/vg/wiki)

# synopsis
This is an incomplete list of vg subcommands. For a complete list, run `vg help`.

- **Graph construction and indexing**
    See the [wiki page](https://github.com/vgteam/vg/wiki/Index-Types) for an overview of vg indexes.
    - [`vg autoindex`](#autoindex): automatically construct a graph and indexes for a specific workflow (e.g. giraffe, rpvg). [wiki page](https://github.com/vgteam/vg/wiki/Automatic-indexing-for-read-mapping-and-downstream-inference)
    - [`vg construct`](#construct): manually construct a graph from a reference and variants. [wiki page](https://github.com/vgteam/vg/wiki/Construction)
    - [`vg index`](#index): manually build individual indexes (xg, distance, GCSA, etc). [wiki page](https://github.com/vgteam/vg/wiki/Index-Construction) 
    - [`vg gbwt`](#gbwt): manually build and manipulate GBWTs and indexes (GBWTgraph, GBZ, r-index). [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg minimizer`](#minimizer): manually build a minimizer index for mapping. 
    - [`vg haplotypes`](#haplotypes): haplotype sample a graph. Recommended for mapping with giraffe. [wiki page](https://github.com/vgteam/vg/wiki/Haplotype-Sampling)
- **Read mapping**
    - [`vg giraffe`](#giraffe): fast haplotype-aware short read alignment. [wiki page](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)
    - [`vg mpmap`](#mpmap): splice-aware multipath alignment of short reads. [wiki page](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-vg-mpmap)
    - [`vg map`](#map): MEM-based read alignment. [wiki page](https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph)
- **Downstream analyses**
    - [`vg pack`](#pack): convert alignments to a compact coverage index. Used with [vg call](#call)
    - [`vg call`](#call): call or genotype VCF variants. Uses [vg pack](#pack). [wiki page](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling)
    - [`vg rna`](#rna): construct splicing graphs and pantranscriptomes. [wiki page](https://github.com/vgteam/vg/wiki/Transcriptomic-analyses). Also see [rpvg](https://github.com/jonassibbesen/rpvg) 
    - [`vg deconstruct`](#deconstruct): create a VCF from variation in the graph. [wiki page](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)
- **Working with read alignments**
    - [`vg gamsort`](#gamsort): sort a GAM/GAF file or index a sorted GAM file.
    - [`vg filter`](#filter): filter alignments by properties.
    - [`vg surject`](#surject): project alignments on a graph onto a linear reference (gam/gaf->bam/sam/cram).
    - [`vg inject`](#inject): project alignments on a linear reference onto a graph (bam/sam/cram->gam/gaf).
    - [`vg sim`](#sim): simulate reads from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Simulating-reads-with-vg-sim)
- **Graph and read statistics**
    - [`vg stats`](#stats): get stats about the graph.
    - [`vg paths`](#paths): get stats about the paths. [wiki page](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)
    - [`vg gbwt`](#gbwt): get stats about a GBWT.
    - [`vg filter`](#filter): get stats about alignments (use `--tsv-out`).
- **Manipulating a graph**
    - [`vg mod`](#mod): filter, transform, and edit the graph.
    - [`vg prune`](#prune): prune the graph for GCSA2 indexing.
    - [`vg ids`](#ids): manipulate graph node ids.
    - [`vg paths`](#paths): manipulate paths in a graph.
    - [`vg gbwt`](#gbwt): manipulate GBWTs and associated indexes. [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg annotate`](#annotate): annotate a graph or alignments.
- **Conversion between formats**
    - [`vg convert`](#convert): convert between handle graph formats and GFA, and between alignment formats.
    - [`vg view`](#view): convert between non-handle graph formats and alignment formats (dot, json, turtle...).
    - [`vg surject`](#surject): project alignments on a graph onto a linear reference (gam/gaf->bam/sam/cram).
    - [`vg inject`](#inject): project alignments on a linear reference onto a graph (bam/sam/cram->gam/gaf).
    - [`vg paths`](#paths): extract a fasta from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Extracting-a-FASTA-from-a-Graph)
- **Subgraph extraction**
    - [`vg chunk`](#chunk): split a graph and/or alignment into smaller chunks.
    - [`vg find`](#find): use an index to find nodes, edges, kmers, paths, or positions.

# annotate

Annotate alignments with graphs and graphs with alignments.

# autoindex

Mapping tool-oriented index construction from interchange formats.

# convert

Convert graphs between handle-graph compliant formats as well as GFA.

# find

Use an index to find nodes, edges, kmers, paths, or positions.

# ids

Manipulate node ids.

# pack

Convert alignments to a compact coverage index.

# paths

Traverse paths in the graph.

# view

format conversions for graphs and alignments

# filter

Filter alignments by properties.

# bugs

Bugs can be reported at: https://github.com/vgteam/vg/issues

For technical support, please visit: https://www.biostars.org/tag/vg/


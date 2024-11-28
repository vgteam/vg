# vg manpage

*Automatically made for vg version v1.61.0-36-g64d7e82e0 "Plodio".*



This is a redundant and incomplete list of subcommands of vg, organized by common uses. For a complete list of subcommands, run `vg help`.

For more in-depth explanations of tools and workflows, see the [vg wiki page](https://github.com/vgteam/vg/wiki)

- **Graph construction and indexing**
    See the [wiki page](https://github.com/vgteam/vg/wiki/Index-Types) for an overview of vg indexes.
    - [`vg autoindex`](#autoindex) automatically construct a graph and indexes for a specific workflow (e.g. giraffe, rpvg). [wiki page](https://github.com/vgteam/vg/wiki/Automatic-indexing-for-read-mapping-and-downstream-inference)
    - [`vg construct`](#construct) manually construct a graph from a reference and variants. [wiki page](https://github.com/vgteam/vg/wiki/Construction)
    - [`vg index`](#index) manually build individual indexes (xg, distance, GCSA, etc). [wiki page](https://github.com/vgteam/vg/wiki/Index-Construction) 
    - [`vg gbwt`](#gbwt) manually build and manipulate GBWTs and indexes (GBWTgraph, GBZ, r-index). [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg minimizer`](#minimizer) manually build a minimizer index for mapping. 
    - [`vg haplotypes`](#haplotypes) haplotype sample a graph. Recommended for mapping with giraffe. [wiki page](https://github.com/vgteam/vg/wiki/Haplotype-Sampling)
- **Read mapping**
    - [`vg giraffe`](#giraffe) fast haplotype-aware short read alignment. [wiki page](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)
    - [`vg mpmap`](#mpmap) splice-aware multipath alignment of short reads. [wiki page](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-vg-mpmap)
    - [`vg map`](#map) MEM-based read alignment. [wiki page](https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph)
- **Downstream analyses**
    - [`vg pack`](#pack) convert alignments to a compact coverage index. Used with [vg call](#call)
    - [`vg call`](#call) call or genotype VCF variants. Uses [vg pack](#pack). [wiki page](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling)
    - [`vg rna`](#rna) construct splicing graphs and pantranscriptomes. [wiki page](https://github.com/vgteam/vg/wiki/Transcriptomic-analyses). Also see [rpvg](https://github.com/jonassibbesen/rpvg) 
    - [`vg deconstruct`](#deconstruct) create a VCF from variation in the graph. [wiki page](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)
- **Working with read alignments**
    - [`vg gamsort`](#gamsort) sort a GAM/GAF file or index a sorted GAM file.
    - [`vg filter`](#filter) filter alignments by properties.
    - [`vg surject`](#surject) project alignments on a graph onto a linear reference (gam/gaf->bam/sam/cram).
    - [`vg inject`](#inject) project alignments on a linear reference onto a graph (bam/sam/cram->gam/gaf).
    - [`vg sim`](#sim) simulate reads from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Simulating-reads-with-vg-sim)
- **Graph and read statistics**
    - [`vg stats`](#stats) get stats about the graph.
    - [`vg paths`](#paths) get stats about the paths. [wiki page](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)
    - [`vg gbwt`](#gbwt) get stats about a GBWT.
    - [`vg filter`](#filter) get stats about alignments (use `--tsv-out`).
- **Manipulating a graph**
    - [`vg mod`](#mod) filter, transform, and edit the graph.
    - [`vg prune`](#prune) prune the graph for GCSA2 indexing.
    - [`vg ids`](#ids) manipulate graph node ids.
    - [`vg paths`](#paths) manipulate paths in a graph.
    - [`vg gbwt`](#gbwt) manipulate GBWTs and associated indexes. [wiki page](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand)
    - [`vg annotate`](#annotate) annotate a graph or alignments.
- **Conversion between formats**
    - [`vg convert`](#convert) convert between handle graph formats and GFA, and between alignment formats.
    - [`vg view`](#view) convert between non-handle graph formats and alignment formats (dot, json, turtle...).
    - [`vg surject`](#surject) project alignments on a graph onto a linear reference (gam/gaf->bam/sam/cram).
    - [`vg inject`](#inject) project alignments on a linear reference onto a graph (bam/sam/cram->gam/gaf).
    - [`vg paths`](#paths) extract a fasta from a graph. [wiki page](https://github.com/vgteam/vg/wiki/Extracting-a-FASTA-from-a-Graph)
- **Subgraph extraction**
    - [`vg chunk`](#chunk) split a graph and/or alignment into smaller chunks.
    - [`vg find`](#find) use an index to find nodes, edges, kmers, paths, or positions.





## annotate



Annotate alignments with graphs and graphs with alignments.





```
usage: vg annotate [options] >output.{gam,vg,tsv}
graph annotation options:
    -x, --xg-name FILE     xg index or graph to annotate (required)
    -b, --bed-name FILE    a BED file to convert to GAM. May repeat.
    -f, --gff-name FILE    a GFF3 file to convert to GAM. May repeat.
    -g, --ggff             output at GGFF subgraph annotation file instead of GAM (requires -s)
    -F, --gaf-output       output in GAF format rather than GAM
    -s, --snarls FILE      file containing snarls to expand GFF intervals into
alignment annotation options:
    -a, --gam FILE         file of Alignments to annotate (required)
    -x, --xg-name FILE     xg index of the graph against which the Alignments are aligned (required)
    -p, --positions        annotate alignments with reference positions
    -m, --multi-position   annotate alignments with multiple reference positions
    -l, --search-limit N   when annotating with positions, search this far for paths (default: read length)
    -b, --bed-name FILE    annotate alignments with overlapping region names from this BED. May repeat.
    -n, --novelty          output TSV table with header describing how much of each Alignment is novel
    -t, --threads          use the specified number of threads

```


## autoindex



Mapping tool-oriented index construction from interchange formats.





```
usage: vg autoindex [options]
options:
  output:
    -p, --prefix PREFIX    prefix to use for all output (default: index)
    -w, --workflow NAME    workflow to produce indexes for, can be provided multiple
                           times. options: map, mpmap, rpvg, giraffe (default: map)
  input data:
    -r, --ref-fasta FILE   FASTA file containing the reference sequence (may repeat)
    -v, --vcf FILE         VCF file with sequence names matching -r (may repeat)
    -i, --ins-fasta FILE   FASTA file with sequences of INS variants from -v
    -g, --gfa FILE         GFA file to make a graph from
    -x, --tx-gff FILE      GTF/GFF file with transcript annotations (may repeat)
    -H, --hap-tx-gff FILE  GTF/GFF file with transcript annotations of a named haplotype (may repeat)
  configuration:
    -f, --gff-feature STR  GTF/GFF feature type (col. 3) to add to graph (default: exon)
    -a, --gff-tx-tag STR   GTF/GFF tag (in col. 9) for transcript ID (default: transcript_id)
  logging and computation:
    -T, --tmp-dir DIR      temporary directory to use for intermediate files
    -M, --target-mem MEM   target max memory usage (not exact, formatted INT[kMG])
                           (default: 1/2 of available)
    -t, --threads NUM      number of threads (default: all available)
    -V, --verbosity NUM    log to stderr (0 = none, 1 = basic, 2 = debug; default 1)
    -h, --help             print this help message to stderr and exit

```


## call


```
usage: vg call [options] <graph> > output.vcf
Call variants or genotype known variants

support calling options:
    -k, --pack FILE          Supports created from vg pack for given input graph
    -m, --min-support M,N    Minimum allele support (M) and minimum site support (N) for call [default = 2,4]
    -e, --baseline-error X,Y Baseline error rates for Poisson model for small (X) and large (Y) variants [default= 0.005,0.01]
    -B, --bias-mode          Use old ratio-based genotyping algorithm as opposed to porbablistic model
    -b, --het-bias M,N       Homozygous alt/ref allele must have >= M/N times more support than the next best allele [default = 6,6]
GAF options:
    -G, --gaf               Output GAF genotypes instead of VCF
    -T, --traversals        Output all candidate traversals in GAF without doing any genotyping
    -M, --trav-padding N    Extend each flank of traversals (from -T) with reference path by N bases if possible
general options:
    -v, --vcf FILE          VCF file to genotype (must have been used to construct input graph with -a)
    -a, --genotype-snarls   Genotype every snarl, including reference calls (use to compare multiple samples)
    -A, --all-snarls        Genotype all snarls, including nested child snarls (like deconstruct -a)
    -c, --min-length N      Genotype only snarls with at least one traversal of length >= N
    -C, --max-length N      Genotype only snarls where all traversals have length <= N
    -f, --ref-fasta FILE    Reference fasta (required if VCF contains symbolic deletions or inversions)
    -i, --ins-fasta FILE    Insertions fasta (required if VCF contains symbolic insertions)
    -s, --sample NAME       Sample name [default=SAMPLE]
    -r, --snarls FILE       Snarls (from vg snarls) to avoid recomputing.
    -g, --gbwt FILE         Only call genotypes that are present in given GBWT index.
    -z, --gbz               Only call genotypes that are present in GBZ index (applies only if input graph is GBZ).
    -N, --translation FILE  Node ID translation (as created by vg gbwt --translation) to apply to snarl names in output
    -O, --gbz-translation   Use the ID translation from the input gbz to apply snarl names to snarl names and AT fields in output
    -p, --ref-path NAME     Reference path to call on (multipile allowed. defaults to all paths)
    -S, --ref-sample NAME   Call on all paths with given sample name (cannot be used with -p)
    -o, --ref-offset N      Offset in reference path (multiple allowed, 1 per path)
    -l, --ref-length N      Override length of reference in the contig field of output VCF
    -d, --ploidy N          Ploidy of sample.  Only 1 and 2 supported. (default: 2)
    -R, --ploidy-regex RULES    use the given comma-separated list of colon-delimited REGEX:PLOIDY rules to assign
                                ploidies to contigs not visited by the selected samples, or to all contigs simulated
                                from if no samples are used. Unmatched contigs get ploidy 2 (or that from -d).
    -n, --nested            Activate nested calling mode (experimental)
    -I, --chains            Call chains instead of snarls (experimental)
        --progress          Show progress
    -t, --threads N         number of threads to use

```


## chunk


```
usage: vg chunk [options] > [chunk.vg]
Splits a graph and/or alignment into smaller chunks

Graph chunks are saved to .vg files, read chunks are saved to .gam files, and haplotype annotations are 
saved to .annotate.txt files, of the form <BASENAME>-<i>-<region name or "ids">-<start>-<length>.<ext>.
The BASENAME is specified with -b and defaults to "./chunks".
For a single-range chunk (-p or -r), the graph data is sent to standard output instead of a file.

options:
    -x, --xg-name FILE       use this graph or xg index to chunk subgraphs
    -G, --gbwt-name FILE     use this GBWT haplotype index for haplotype extraction (for -T)
    -a, --gam-name FILE      chunk this gam file instead of the graph (multiple allowed)
    -g, --gam-and-graph      when used in combination with -a, both gam and graph will be chunked
    -F, --in-gaf             input alignment is a sorted bgzipped GAF
path chunking:
    -p, --path TARGET        write the chunk in the specified (0-based inclusive, multiple allowed)
                             path range TARGET=path[:pos1[-pos2]] to standard output
    -P, --path-list FILE     write chunks for all path regions in (line - separated file). format
                             for each as in -p (all paths chunked unless otherwise specified)
    -e, --input-bed FILE     write chunks for all (0-based end-exclusive) bed regions
    -S, --snarls FILE        write given path-range(s) and all snarls fully contained in them, as alternative to -c
id range chunking:
    -r, --node-range N:M     write the chunk for the specified node range to standard output
    -R, --node-ranges FILE   write the chunk for each node range in (newline or whitespace separated) file
    -n, --n-chunks N         generate this many id-range chunks, which are determined using the xg index
simple gam chunking:
    -m, --gam-split-size N   split gam (specified with -a, sort/index not required) up into chunks with at most N reads each
component chunking:
    -C, --components         create a chunk for each connected component.  if a targets given with (-p, -P, -r, -R), limit to components containing them
    -M, --path-components    create a chunk for each path in the graph's connected component
general:
    -s, --chunk-size N       create chunks spanning N bases (or nodes with -r/-R) for all input regions.
    -o, --overlap N          overlap between chunks when using -s [0]
    -E, --output-bed FILE    write all created chunks to a bed file
    -b, --prefix BASENAME    write output chunk files with the given base name. Files for chunk i will
                             be named: <BASENAME>-<i>-<name>-<start>-<length>.<ext> [./chunk]
    -c, --context-steps N    expand the context of the chunk this many node steps [1]
    -l, --context-length N   expand the context of the chunk by this many bp [0]
    -T, --trace              trace haplotype threads in chunks (and only expand forward from input coordinates).
                             Produces a .annotate.txt file with haplotype frequencies for each chunk.
    --no-embedded-haplotypes Don't load haplotypes from the graph. It is possible to -T without any haplotypes available.
    -f, --fully-contained    only return GAM alignments that are fully contained within chunk
    -O, --output-fmt         Specify output format (vg, pg, hg, gfa).  [pg (vg with -T)]
    -t, --threads N          for tasks that can be done in parallel, use this many threads [1]
    -h, --help

```


## construct


```
usage: vg construct [options] >new.vg
options:
construct from a reference and variant calls:
    -r, --reference FILE   input FASTA reference (may repeat)
    -v, --vcf FILE         input VCF (may repeat)
    -n, --rename V=F       match contig V in the VCFs to contig F in the FASTAs (may repeat)
    -a, --alt-paths        save paths for alts of variants by SHA1 hash
    -A, --alt-paths-plain  save paths for alts of variants by variant ID if possible, otherwise SHA1
                           (IDs must be unique across all input VCFs)
    -R, --region REGION    specify a VCF contig name or 1-based inclusive region (may repeat, if on different contigs)
    -C, --region-is-chrom  don't attempt to parse the regions (use when the reference
                           sequence name could be inadvertently parsed as a region)
    -z, --region-size N    variants per region to parallelize (default: 1024)
    -t, --threads N        use N threads to construct graph (defaults to numCPUs)
    -S, --handle-sv        include structural variants in construction of graph.
    -I, --insertions FILE  a FASTA file containing insertion sequences 
                           (referred to in VCF) to add to graph.
    -f, --flat-alts N      don't chop up alternate alleles from input VCF
    -l, --parse-max N      don't chop up alternate alleles from input VCF longer than N (default: 100)
    -i, --no-trim-indels   don't remove the 1bp reference base from alt alleles of indels.
    -N, --in-memory        construct the entire graph in memory before outputting it.
construct from a multiple sequence alignment:
    -M, --msa FILE         input multiple sequence alignment
    -F, --msa-format       format of the MSA file (options: fasta, clustal; default fasta)
    -d, --drop-msa-paths   don't add paths for the MSA sequences into the graph
shared construction options:
    -m, --node-max N       limit the maximum allowable node sequence size (default: 32)
                           nodes greater than this threshold will be divided
                           Note: nodes larger than ~1024 bp can't be GCSA2-indexed
    -p, --progress         show progress

```


## convert



Convert graphs between handle-graph compliant formats as well as GFA.





```
usage: vg convert [options] <input-graph>
input options:
    -g, --gfa-in           input in GFA format
    -r, --in-rgfa-rank N   import rgfa tags with rank <= N as paths [default=0]
    -b, --gbwt-in FILE     input graph is a GBWTGraph using the GBWT in FILE
        --ref-sample STR   change haplotypes for this sample to reference paths (may repeat)
gfa input options (use with -g):
    -T, --gfa-trans FILE   write gfa id conversions to FILE
output options:
    -v, --vg-out           output in VG's original Protobuf format [DEPRECATED: use -p instead].
    -a, --hash-out         output in HashGraph format
    -p, --packed-out       output in PackedGraph format [default]
    -x, --xg-out           output in XG format
    -f, --gfa-out          output in GFA format
    -H, --drop-haplotypes  do not include haplotype paths in the output
                           (useful with GBWTGraph / GBZ inputs)
gfa output options (use with -f):
    -P, --rgfa-path STR    write given path as rGFA tags instead of lines
                           (multiple allowed, only rank-0 supported)
    -Q, --rgfa-prefix STR  write paths with given prefix as rGFA tags instead of lines
                           (multiple allowed, only rank-0 supported)
    -B, --rgfa-pline       paths written as rGFA tags also written as lines
    -W, --no-wline         Write all paths as GFA P-lines instead of W-lines.
                           Allows handling multiple phase blocks and subranges used together.
    --gbwtgraph-algorithm  Always use the GBWTGraph library GFA algorithm.
                           Not compatible with other GFA output options or non-GBWT graphs.
    --vg-algorithm         Always use the VG GFA algorithm. Works with all options and graph types,
                           but can't preserve original GFA coordinates.
    --no-translation       When using the GBWTGraph algorithm, convert the graph directly to GFA.
                           Do not use the translation to preserve original coordinates.
alignment options:
    -G, --gam-to-gaf FILE  convert GAM FILE to GAF
    -F, --gaf-to-gam FILE  convert GAF FILE to GAM
general options:
    -t, --threads N        use N threads (defaults to numCPUs)

```


## deconstruct


```
usage: vg deconstruct [options] [-p|-P] <PATH> <GRAPH>
Outputs VCF records for Snarls present in a graph (relative to a chosen reference path).
options: 
    -p, --path NAME          A reference path to deconstruct against (multiple allowed).
    -P, --path-prefix NAME   All paths [excluding GBWT threads / non-reference GBZ paths] beginning with NAME used as reference (multiple allowed).
                             Other non-ref paths not considered as samples. 
    -r, --snarls FILE        Snarls file (from vg snarls) to avoid recomputing.
    -g, --gbwt FILE          consider alt traversals that correspond to GBWT haplotypes in FILE (not needed for GBZ graph input).
    -T, --translation FILE   Node ID translation (as created by vg gbwt --translation) to apply to snarl names and AT fields in output
    -O, --gbz-translation    Use the ID translation from the input gbz to apply snarl names to snarl names and AT fields in output
    -a, --all-snarls         Process all snarls, including nested snarls (by default only top-level snarls reported).
    -c, --context-jaccard N  Set context mapping size used to disambiguate alleles at sites with multiple reference traversals (default: 10000).
    -u, --untangle-travs     Use context mapping to determine the reference-relative positions of each step in allele traversals (AP INFO field).
    -K, --keep-conflicted    Retain conflicted genotypes in output.
    -S, --strict-conflicts   Drop genotypes when we have more than one haplotype for any given phase (set by default when using GBWT input).
    -C, --contig-only-ref    Only use the CONTIG name (and not SAMPLE#CONTIG#HAPLOTYPE etc) for the reference if possible (ie there is only one reference sample).
    -L, --cluster F          Cluster traversals whose (handle) Jaccard coefficient is >= F together (default: 1.0) [experimental]
    -n, --nested             Write a nested VCF, including special tags. [experimental]
    -R, --star-allele        Use *-alleles to denote alleles that span but do not cross the site. Only works with -n
    -t, --threads N          Use N threads
    -v, --verbose            Print some status messages


```


## filter



Filter alignments by properties.

`vg filter --tsv-out` can be used to produce a TSV file of user-specified fields from the GAM file. For example, 

`vg filter --tsv-out "name;mapping_quality" <alignment.gam>`

is the equivalent of 

`vg view -aj <alignment.gam> | jq -r '[.name,.mapping_quality] | @tsv'`

To find which fields are stored in a GAM file, use [`vg view`](#view) to view the GAM as a JSON file. 




```
vg: invalid option -- 'h'
usage: vg filter [options] <alignment.gam> > out.gam
Filter alignments by properties.

options:
    -M, --input-mp-alns        input is multipath alignments (GAMP) rather than GAM
    -n, --name-prefix NAME     keep only reads with this prefix in their names [default='']
    -N, --name-prefixes FILE   keep reads with names with one of many prefixes, one per nonempty line
    -e, --exact-name           match read names exactly instead of by prefix
    -a, --subsequence NAME     keep reads that contain this subsequence
    -A, --subsequences FILE    keep reads that contain one of these subsequences, one per nonempty line
    -p, --proper-pairs         keep reads that are annotated as being properly paired
    -P, --only-mapped          keep reads that are mapped
    -X, --exclude-contig REGEX drop reads with refpos annotations on contigs matching the given regex (may repeat)
    -F, --exclude-feature NAME drop reads with the given feature in the "features" annotation (may repeat)
    -s, --min-secondary N      minimum score to keep secondary alignment
    -r, --min-primary N        minimum score to keep primary alignment
    -O, --rescore              re-score reads using default parameters and only alignment information
    -f, --frac-score           normalize score based on length
    -u, --substitutions        use substitution count instead of score
    -o, --max-overhang N       filter reads whose alignments begin or end with an insert > N [default=99999]
    -m, --min-end-matches N    filter reads that don't begin with at least N matches on each end
    -S, --drop-split           remove split reads taking nonexistent edges
    -x, --xg-name FILE         use this xg index or graph (required for -S and -D)
    -v, --verbose              print out statistics on numbers of reads filtered by what.
    -V, --no-output            print out statistics (as above) but do not write out filtered GAM.
    -T, --tsv-out FIELD[;FIELD] do not write filtered gam but a tsv of the given fields
    -q, --min-mapq N           filter alignments with mapping quality < N
    -E, --repeat-ends N        filter reads with tandem repeat (motif size <= 2N, spanning >= N bases) at either end
    -D, --defray-ends N        clip back the ends of reads that are ambiguously aligned, up to N bases
    -C, --defray-count N       stop defraying after N nodes visited (used to keep runtime in check) [default=99999]
    -d, --downsample S.P       filter out all but the given portion 0.P of the reads. S may be an integer seed as in SAMtools
    -i, --interleaved          assume interleaved input. both ends will be filtered out if either fails filter
    -I, --interleaved-all      assume interleaved input. both ends will be filtered out if *both* fail filters
    -b, --min-base-quality Q:F filter reads with where fewer than fraction F bases have base quality >= PHRED score Q.
    -B, --annotation K[:V]     keep reads if the annotation is present. If a value is given, keep reads if the values are equal
                               similar to running jq 'select(.annotation.K==V)' on the json
    -c, --correctly-mapped     keep only reads that are marked as correctly-mapped
    -U, --complement           apply the complement of the filter implied by the other arguments.
    -t, --threads N            number of threads [1]

```


## find



Use an index to find nodes, edges, kmers, paths, or positions.





```
usage: vg find [options] >sub.vg
options:
graph features:
    -x, --xg-name FILE     use this xg index or graph (instead of rocksdb db)
    -n, --node ID          find node(s), return 1-hop context as graph
    -N, --node-list FILE   a white space or line delimited list of nodes to collect
        --mapping FILE     also include nodes that map to the selected node ids
    -e, --edges-end ID     return edges on end of node with ID
    -s, --edges-start ID   return edges on start of node with ID
    -c, --context STEPS    expand the context of the subgraph this many steps
    -L, --use-length       treat STEPS in -c or M in -r as a length in bases
    -P, --position-in PATH find the position of the node (specified by -n) in the given path
    -I, --list-paths       write out the path names in the index
    -r, --node-range N:M   get nodes from N to M
    -G, --gam GAM          accumulate the graph touched by the alignments in the GAM
    --connecting-start POS find the graph connecting from POS (node ID, + or -, node offset) to --connecting-end
    --connecting-end POS   find the graph connecting to POS (node ID, + or -, node offset) from --connecting-start
    --connecting-range INT traverse up to INT bases when going from --connecting-start to --connecting-end (default: 100)
subgraphs by path range:
    -p, --path TARGET      find the node(s) in the specified path range(s) TARGET=path[:pos1[-pos2]]
    -R, --path-bed FILE    read our targets from the given BED FILE
    -E, --path-dag         with -p or -R, gets any node in the partial order from pos1 to pos2, assumes id sorted DAG
    -W, --save-to PREFIX   instead of writing target subgraphs to stdout,
                           write one per given target to a separate file named PREFIX[path]:[start]-[end].vg
    -K, --subgraph-k K     instead of graphs, write kmers from the subgraphs
    -H, --gbwt FILE        when enumerating kmers from subgraphs, determine their frequencies in this GBWT haplotype index
alignments:
    -l, --sorted-gam FILE  use this sorted, indexed GAM file
    -F, --sorted-gaf FILE  use this sorted, indexed GAF file
    -o, --alns-on N:M      write alignments which align to any of the nodes between N and M (inclusive)
    -A, --to-graph VG      get alignments to the provided subgraph
sequences:
    -g, --gcsa FILE        use this GCSA2 index of the sequence space of the graph (required for sequence queries)
    -S, --sequence STR     search for sequence STR using
    -M, --mems STR         describe the super-maximal exact matches of the STR (gcsa2) in JSON
    -B, --reseed-length N  find non-super-maximal MEMs inside SMEMs of length at least N
    -f, --fast-reseed      use fast SMEM reseeding algorithm
    -Y, --max-mem N        the maximum length of the MEM (default: GCSA2 order)
    -Z, --min-mem N        the minimum length of the MEM (default: 1)
    -D, --distance         return distance on path between pair of nodes (-n). if -P not used, best path chosen heurstically
    -Q, --paths-named S    return all paths whose names are prefixed with S (multiple allowed)

```


## gamsort


```
gamsort: sort a GAM/GAF file, or index a sorted GAM file
Usage: gamsort [Options] gamfile
Options:
  -i / --index FILE       produce an index of the sorted GAM file
  -d / --dumb-sort        use naive sorting algorithm (no tmp files, faster for small GAMs)
  -p / --progress         Show progress.
  -G / --gaf-input        Input is a GAF file.
  -c / --chunk-size       Number of reads per chunk when sorting GAFs.
  -t / --threads          Use the specified number of threads.


```


## gbwt


```
usage: vg gbwt [options] [args]

Manipulate GBWTs. Input GBWTs are loaded from input args or built in earlier steps.
The input graph is provided with one of -x, -G, or -Z

General options:
    -x, --xg-name FILE      read the graph from FILE
    -o, --output FILE       write output GBWT to FILE
    -d, --temp-dir DIR      use directory DIR for temporary files
    -p, --progress          show progress and statistics

GBWT construction parameters (for steps 1 and 4):
        --buffer-size N     GBWT construction buffer size in millions of nodes (default 100)
        --id-interval N     store path ids at one out of N positions (default 1024)

Multithreading:
        --num-jobs N        use at most N parallel build jobs (for -v, -G, -A, -l, -P; default 4)
        --num-threads N     use N parallel search threads (for -b and -r; default 8)

Step 1: GBWT construction (requires -o and one of { -v, -G, -Z, -E, A }):
    -v, --vcf-input         index the haplotypes in the VCF files specified in input args in parallel
                            (inputs must be over different contigs; requires -x, implies -f)
                            (does not store graph contigs in the GBWT)
        --preset X          use preset X (available: 1000gp)
        --inputs-as-jobs    create one build job for each input instead of using first-fit heuristic
        --parse-only        store the VCF parses without building GBWTs
                            (use -o for the file name prefix; skips subsequent steps)
        --ignore-missing    do not warn when variants are missing from the graph
        --actual-phasing    do not interpret unphased homozygous genotypes as phased
        --force-phasing     replace unphased genotypes with randomly phased ones
        --discard-overlaps  skip overlapping alternate alleles if the overlap cannot be resolved
                            instead of creating a phase break
        --batch-size N      index the haplotypes in batches of N samples (default 200)
        --sample-range X-Y  index samples X to Y (inclusive, 0-based)
        --rename V=P        VCF contig V matches path P in the graph (may repeat)
        --vcf-variants      variants in the graph use VCF contig names instead of path names
        --vcf-region C:X-Y  restrict VCF contig C to coordinates X to Y (inclusive, 1-based; may repeat)
        --exclude-sample X  do not index the sample with name X (faster than -R; may repeat)
    -G, --gfa-input         index the walks or paths in the GFA file (one input arg)
        --max-node N        chop long segments into nodes of at most N bp (default 1024, use 0 to disable)
        --path-regex X      parse metadata as haplotypes from path names using regex X instead of vg-parser-compatible rules
        --path-fields X     parse metadata as haplotypes, mapping regex submatches to these fields instead of using vg-parser-compatible rules
        --translation FILE  write the segment to node translation table to FILE
    -Z, --gbz-input         extract GBWT and GBWTGraph from GBZ input (one input arg)
        --translation FILE  write the segment to node translation table to FILE
    -I, --gg-in FILE        load GBWTGraph from FILE and GBWT from input (one input arg) 
    -E, --index-paths       index the embedded non-alt paths in the graph (requires -x, no input args)
    -A, --alignment-input   index the alignments in the GAF files specified in input args (requires -x)
        --gam-format        the input files are in GAM format instead of GAF format

Step 2: Merge multiple input GBWTs (requires -o):
    -m, --merge             use the insertion algorithm
    -f, --fast              fast merging algorithm (node ids must not overlap)
    -b, --parallel          use the parallel algorithm
        --chunk-size N      search in chunks of N sequences (default 1)
        --pos-buffer N      use N MiB position buffers for each search thread (default 64)
        --thread-buffer N   use N MiB thread buffers for each search thread (default 256)
        --merge-buffers N   merge 2^N thread buffers into one file per merge job (default 6)
        --merge-jobs N      run N parallel merge jobs (default 4)

Step 3: Alter GBWT (requires -o and one input GBWT):
    -R, --remove-sample X   remove the sample with name X from the index (may repeat)
        --set-tag K=V       set a GBWT tag (may repeat)
        --set-reference X   set sample X as the reference (may repeat)

Step 4: Path cover GBWT construction (requires an input graph, -o, and one of { -a, -l, -P }):
    -a, --augment-gbwt      add a path cover of missing components (one input GBWT)
    -l, --local-haplotypes  sample local haplotypes (one input GBWT)
    -P, --path-cover        build a greedy path cover (no input GBWTs)
    -n, --num-paths N       find N paths per component (default 64 for -l, 16 otherwise)
    -k, --context-length N  use N-node contexts (default 4)
        --pass-paths        include named graph paths in local haplotype or greedy path cover GBWT

Step 5: GBWTGraph construction (requires an input graph and one input GBWT):
    -g, --graph-name FILE   build GBWTGraph and store it in FILE
        --gbz-format        serialize both GBWT and GBWTGraph in GBZ format (makes -o unnecessary)

Step 6: R-index construction (one input GBWT):
    -r, --r-index FILE      build an r-index and store it in FILE

Step 7: Metadata (one input GBWT):
    -M, --metadata          print basic metadata
    -C, --contigs           print the number of contigs
    -H, --haplotypes        print the number of haplotypes
    -S, --samples           print the number of samples
    -L, --list-names        list contig/sample names (use with -C or -S)
    -T, --path-names        list path names
        --tags              list GBWT tags

Step 8: Paths (one input GBWT):
    -c, --count-paths       print the number of paths
    -e, --extract FILE      extract paths in SDSL format to FILE


```


## giraffe


```
usage:
  vg giraffe -Z graph.gbz [-d graph.dist -m graph.min] <input options> [other options] > output.gam
  vg giraffe -Z graph.gbz --haplotype-name graph.hapl --kff-name sample.kff <input options> [other options] > output.gam

Fast haplotype-aware short read mapper.

basic options:
  -Z, --gbz-name FILE           map to this GBZ graph
  -d, --dist-name FILE          cluster using this distance index
  -m, --minimizer-name FILE     use this minimizer index
  -p, --progress                show progress
  -t, --threads INT             number of mapping threads to use
  -b, --parameter-preset NAME   set computational parameters (fast / default) [default]
  -h, --help                    print full help with all available options
input options:
  -G, --gam-in FILE             read and realign GAM-format reads from FILE
  -f, --fastq-in FILE           read and align FASTQ-format reads from FILE (two are allowed, one for each mate)
  -i, --interleaved             GAM/FASTQ input is interleaved pairs, for paired-end alignment
haplotype sampling:
  --haplotype-name FILE         sample from haplotype information in FILE
  --kff-name FILE               sample according to kmer counts in FILE
  --index-basename STR          name prefix for generated graph/index files (default: from graph name)
alternate graphs:
  -x, --xg-name FILE            map to this graph (if no -Z / -g), or use this graph for HTSLib output
  -g, --graph-name FILE         map to this GBWTGraph (if no -Z)
  -H, --gbwt-name FILE          use this GBWT index (when mapping to -x / -g)
output options:
  -N, --sample NAME             add this sample name
  -R, --read-group NAME         add this read group
  -o, --output-format NAME      output the alignments in NAME format (gam / gaf / json / tsv / SAM / BAM / CRAM) [gam]
  --ref-paths FILE              ordered list of paths in the graph, one per line or HTSlib .dict, for HTSLib @SQ headers
  --named-coordinates           produce GAM/GAF outputs in named-segment (GFA) space
  -P, --prune-low-cplx          prune short and low complexity anchors during linear format realignment
  -n, --discard                 discard all output alignments (for profiling)
  --output-basename NAME        write output to a GAM file beginning with the given prefix for each setting combination
  --report-name NAME            write a TSV of output file and mapping speed to the given file
  --show-work                   log how the mapper comes to its conclusions about mapping locations
Giraffe parameters:
  -A, --rescue-algorithm NAME   use algorithm NAME for rescue (none / dozeu / gssw) [dozeu]
  --fragment-mean FLOAT         force the fragment length distribution to have this mean (requires --fragment-stdev)
  --fragment-stdev FLOAT        force the fragment length distribution to have this standard deviation (requires --fragment-mean)
  --track-provenance            track how internal intermediate alignment candidates were arrived at
  --track-correctness           track if internal intermediate alignment candidates are correct (implies --track-provenance)
  -B, --batch-size INT          number of reads or pairs per batch to distribute to threads [512]
program options:
  --watchdog-timeout INT                           complain after INT seconds working on a read or read pair [10]
scoring options:
  --match INT                                      use this match score [1]
  --mismatch INT                                   use this mismatch penalty [4]
  --gap-open INT                                   use this gap open penalty [6]
  --gap-extend INT                                 use this gap extension penalty [1]
  --full-l-bonus INT                               the full-length alignment bonus [5]
result options:
  -M, --max-multimaps INT                          produce up to INT alignments for each read [1]
computational parameters:
  -c, --hit-cap INT                                use all minimizers with at most INT hits [10]
  -C, --hard-hit-cap INT                           ignore all minimizers with more than INT hits [500]
  -F, --score-fraction FLOAT                       select minimizers between hit caps until score is FLOAT of total [0.9]
  -U, --max-min INT                                use at most INT minimizers [500]
  --num-bp-per-min INT                             use maximum of number minimizers calculated by READ_LENGTH / INT and --max-min [1000]
  -D, --distance-limit INT                         cluster using this distance limit [200]
  -e, --max-extensions INT                         extend up to INT clusters [800]
  -a, --max-alignments INT                         align up to INT extensions [8]
  -s, --cluster-score FLOAT                        only extend clusters if they are within INT of the best score [50]
  -S, --pad-cluster-score FLOAT                    also extend clusters within INT of above threshold to get a second-best cluster [20]
  -u, --cluster-coverage FLOAT                     only extend clusters if they are within FLOAT of the best read coverage [0.3]
  -v, --extension-score INT                        only align extensions if their score is within INT of the best score [1]
  -w, --extension-set FLOAT                        only align extension sets if their score is within INT of the best score [20]
  -O, --no-dp                                      disable all gapped alignment
  -r, --rescue-attempts INT                        attempt up to INT rescues per read in a pair [15]
  -L, --max-fragment-length INT                    assume that fragment lengths should be smaller than INT when estimating the fragment length distribution [2000]
  --exclude-overlapping-min                        exclude overlapping minimizers
  --paired-distance-limit FLOAT                    cluster pairs of read using a distance limit FLOAT standard deviations greater than the mean [2]
  --rescue-subgraph-size FLOAT                     search for rescued alignments FLOAT standard deviations greater than the mean [4]
  --rescue-seed-limit INT                          attempt rescue with at most INT seeds [100]
long-read/chaining parameters:
  --align-from-chains                              chain up extensions to create alignments, instead of doing each separately
  --chaining-cluster-distance INT                  maximum distance to cluster over before chaining [100]
  --precluster-connection-coverage-threshold FLOAT threshold of precluster pair coverage below the base, after which to stop reseeding between preclusters [0.3]
  --min-precluster-connections INT                 minimum number of precluster connections to reseed over [10]
  --max-precluster-connections INT                 maximum number of precluster connections to reseed over [50]
  --max-lookback-bases INT                         maximum distance to look back when chaining [100]
  --min-lookback-items INT                         minimum items to consider coming from when chaining [1]
  --lookback-item-hard-cap INT                     maximum items to consider coming from when chaining [15]
  --chain-score-threshold FLOAT                    only align chains if their score is within this many points of the best score [100]
  --min-chains INT                                 ignore score threshold to get this many chains aligned [1]
  --chain-min-score INT                            do not align chains with less than this score [100]
  --max-chain-connection INT                       maximum distance across which to connect seeds when aligning a chain [100]
  --max-tail-length INT                            maximum length of a tail to align before forcing softclipping when aligning a chain [100]
  --max-dp-cells INT                               maximum number of alignment cells to allow in a tail with GSSW [16777216]

```


## haplotypes


```
Usage:
    vg haplotypes [options] -k kmers.kff -g output.gbz graph.gbz
    vg haplotypes [options] -H output.hapl graph.gbz
    vg haplotypes [options] -i graph.hapl -k kmers.kff -g output.gbz graph.gbz
    vg haplotypes [options] -i graph.hapl --vcf-input variants.vcf graph.gbz > output.tsv
    vg haplotypes [options] -i graph.hapl -k kmers.kff --extract M:N graph.gbz > output.fa

Haplotype sampling based on kmer counts.

Output files:
    -g, --gbz-output X        write the output GBZ to X
    -H, --haplotype-output X  write haplotype information to X

Input files:
    -d, --distance-index X    use this distance index (default: <basename>.dist)
    -r, --r-index X           use this r-index (default: <basename>.ri)
    -i, --haplotype-input X   use this haplotype information (default: generate)
    -k, --kmer-input X        use kmer counts from this KFF file (required for --gbz-output)

Options for generating haplotype information:
        --kmer-length N       kmer length for building the minimizer index (default: 29)
        --window-length N     window length for building the minimizer index (default: 11)
        --subchain-length N   target length (in bp) for subchains (default: 10000)
        --linear-structure    extend subchains to avoid haplotypes visiting them multiple times

Options for sampling haplotypes:
        --preset X            use preset X (default, haploid, diploid)
        --coverage N          kmer coverage in the KFF file (default: estimate)
        --num-haplotypes N    generate N haplotypes (default: 4)
                              sample from N candidates (with --diploid-sampling; default: 32)
        --present-discount F  discount scores for present kmers by factor F (default: 0.9)
        --het-adjustment F    adjust scores for heterozygous kmers by F (default: 0.05)
        --absent-score F      score absent kmers -F/+F (default: 0.8)
        --haploid-scoring     use a scoring model without heterozygous kmers
        --diploid-sampling    choose the best pair from the sampled haplotypes
        --include-reference   include named and reference paths in the output

Other options:
    -v, --verbosity N         verbosity level (0 = silent, 1 = basic, 2 = detailed, 3 = debug; default: 0)
    -t, --threads N           approximate number of threads (default: 8 on this system)

Developer options:
        --validate            validate the generated information (may be slow)
        --vcf-input X         map the variants in VCF file X to subchains
        --contig-prefix X     a prefix for transforming VCF contig names into GBWT contig names
        --extract M:N         extract haplotypes in chain M, subchain N in FASTA format
        --score-output X      write haplotype scores to X
        --classify X          classify kmers and write output to X


```


## ids



Manipulate node ids.





```
usage: vg ids [options] <graph1.vg> [graph2.vg ...] >new.vg
options:
    -c, --compact        minimize the space of integers used by the ids
    -i, --increment N    increase ids by N
    -d, --decrement N    decrease ids by N
    -j, --join           make a joint id space for all the graphs that are supplied
                         by iterating through the supplied graphs and incrementing
                         their ids to be non-conflicting (modifies original files)
    -m, --mapping FILE   create an empty node mapping for vg prune
    -s, --sort           assign new node IDs in (generalized) topological sort order

```


## index


```
usage: vg index [options] <graph1.vg> [graph2.vg ...]
Creates an index on the specified graph or graphs. All graphs indexed must 
already be in a joint ID space.
general options:
    -b, --temp-dir DIR        use DIR for temporary files
    -t, --threads N           number of threads to use
    -p, --progress            show progress
xg options:
    -x, --xg-name FILE        use this file to store a succinct, queryable version of the graph(s), or read for GCSA or distance indexing
    -L, --xg-alts             include alt paths in xg
gcsa options:
    -g, --gcsa-out FILE       output a GCSA2 index to the given file
    -f, --mapping FILE        use this node mapping in GCSA2 construction
    -k, --kmer-size N         index kmers of size N in the graph (default 16)
    -X, --doubling-steps N    use this number of doubling steps for GCSA2 construction (default 4)
    -Z, --size-limit N        limit temporary disk space usage to N gigabytes (default 2048)
    -V, --verify-index        validate the GCSA2 index using the input kmers (important for testing)
gam indexing options:
    -l, --index-sorted-gam    input is sorted .gam format alignments, store a GAI index of the sorted GAM in INPUT.gam.gai
vg in-place indexing options:
    --index-sorted-vg         input is ID-sorted .vg format graph chunks, store a VGI index of the sorted vg in INPUT.vg.vgi
snarl distance index options
    -j  --dist-name FILE      use this file to store a snarl-based distance index
        --snarl-limit N       don't store snarl distances for snarls with more than N nodes (default 10000)
                              if N is 0 then don't store distances, only the snarl tree
        --no-nested-distance  only store distances along the top-level chain

```


## map


```
vg: invalid option -- 'h'
usage: vg map [options] -d idxbase -f in1.fq [-f in2.fq] >aln.gam
Align reads to a graph.

graph/index:
    -d, --base-name BASE          use BASE.xg and BASE.gcsa as the input index pair
    -x, --xg-name FILE            use this xg index or graph (defaults to <graph>.vg.xg)
    -g, --gcsa-name FILE          use this GCSA2 index (defaults to <graph>.gcsa)
    -1, --gbwt-name FILE          use this GBWT haplotype index (defaults to <graph>.gbwt)
algorithm:
    -t, --threads N               number of compute threads to use
    -k, --min-mem INT             minimum MEM length (if 0 estimate via -e) [0]
    -e, --mem-chance FLOAT        set {-k} such that this fraction of {-k} length hits will by chance [5e-4]
    -c, --hit-max N               ignore MEMs who have >N hits in our index (0 for no limit) [2048]
    -Y, --max-mem INT             ignore mems longer than this length (unset if 0) [0]
    -r, --reseed-x FLOAT          look for internal seeds inside a seed longer than FLOAT*--min-seed [1.5]
    -u, --try-up-to INT           attempt to align up to the INT best candidate chains of seeds (1/2 for paired) [128]
    -l, --try-at-least INT        attempt to align at least the INT best candidate chains of seeds [1]
    -E, --approx-mq-cap INT       weight MQ by suffix tree based estimate when estimate less than FLOAT [0]
    --id-mq-weight N              scale mapping quality by the alignment score identity to this power [2]
    -W, --min-chain INT           discard a chain if seeded bases shorter than INT [0]
    -C, --drop-chain FLOAT        drop chains shorter than FLOAT fraction of the longest overlapping chain [0.45]
    -n, --mq-overlap FLOAT        scale MQ by count of alignments with this overlap in the query with the primary [0]
    -P, --min-ident FLOAT         accept alignment only if the alignment identity is >= FLOAT [0]
    -H, --max-target-x N          skip cluster subgraphs with length > N*read_length [100]
    -w, --band-width INT          band width for long read alignment [256]
    -O, --band-overlap INT        band overlap for long read alignment [{-w}/8]
    -J, --band-jump INT           the maximum number of bands of insertion we consider in the alignment chain model [128]
    -B, --band-multi INT          consider this many alignments of each band in banded alignment [16]
    -Z, --band-min-mq INT         treat bands with less than this MQ as unaligned [0]
    -I, --fragment STR            fragment length distribution specification STR=m:μ:σ:o:d [5000:0:0:0:1]
                                  max, mean, stdev, orientation (1=same, 0=flip), direction (1=forward, 0=backward)
    -U, --fixed-frag-model        don't learn the pair fragment model online, use {-I} without update
    -p, --print-frag-model        suppress alignment output and print the fragment model on stdout as per {-I} format
    --frag-calc INT               update the fragment model every INT perfect pairs [10]
    --fragment-x FLOAT            calculate max fragment size as frag_mean+frag_sd*FLOAT [10]
    --mate-rescues INT            attempt up to INT mate rescues per pair [64]
    -S, --unpaired-cost INT       penalty for an unpaired read pair [17]
    --no-patch-aln                do not patch banded alignments by locally aligning unaligned regions
    --xdrop-alignment             use X-drop heuristic (much faster for long-read alignment)
    --max-gap-length              maximum gap length allowed in each contiguous alignment (for X-drop alignment) [40]
scoring:
    -q, --match INT               use this match score [1]
    -z, --mismatch INT            use this mismatch penalty [4]
    --score-matrix FILE           read a 4x4 integer substitution scoring matrix from a file
    -o, --gap-open INT            use this gap open penalty [6]
    -y, --gap-extend INT          use this gap extension penalty [1]
    -L, --full-l-bonus INT        the full-length alignment bonus [5]
    --drop-full-l-bonus           remove the full length bonus from the score before sorting and MQ calculation
    -a, --hap-exp FLOAT           the exponent for haplotype consistency likelihood in alignment score [1]
    --recombination-penalty FLOAT use this log recombination penalty for GBWT haplotype scoring [20.7]
    -A, --qual-adjust             perform base quality adjusted alignments (requires base quality input)
preset:
    -m, --alignment-model STR     use a preset alignment scoring model, either "short" (default) or "long" (for ONT/PacBio)
                                  "long" is equivalent to `-u 2 -L 63 -q 1 -z 2 -o 2 -y 1 -w 128 -O 32`
input:
    -s, --sequence STR            align a string to the graph in graph.vg using partial order alignment
    -V, --seq-name STR            name the sequence using this value (for graph modification with new named paths)
    -T, --reads FILE              take reads (one per line) from FILE, write alignments to stdout
    -b, --hts-input FILE          align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout
    -G, --gam-input FILE          realign GAM input
    -f, --fastq FILE              input fastq or (2-line format) fasta, possibly compressed, two are allowed, one for each mate
    -F, --fasta FILE              align the sequences in a FASTA file that may have multiple lines per reference sequence
    -i, --interleaved             fastq or GAM is interleaved paired-ended
    -N, --sample NAME             for --reads input, add this sample
    -R, --read-group NAME         for --reads input, add this read group
output:
    -j, --output-json             output JSON rather than an alignment stream (helpful for debugging)
    -%, --gaf                     output alignments in GAF format
    --surject-to TYPE             surject the output into the graph's paths, writing TYPE := bam |sam | cram
    --ref-paths FILE              ordered list of paths in the graph, one per line or HTSlib .dict, for HTSLib @SQ headers
    --buffer-size INT             buffer this many alignments together before outputting in GAM [512]
    -X, --compare                 realign GAM input (-G), writing alignment with "correct" field set to overlap with input
    -v, --refpos-table            for efficient testing output a table of name, chr, pos, mq, score
    -K, --keep-secondary          produce alignments for secondary input alignments in addition to primary ones
    -M, --max-multimaps INT       produce up to INT alignments for each read [1]
    -Q, --mq-max INT              cap the mapping quality at INT [60]
    --exclude-unaligned           exclude reads with no alignment
    -D, --debug                   print debugging information about alignment to stderr
    --log-time                    print runtime to stderr

```


## minimizer


```
usage: vg minimizer [options] -d graph.dist -o graph.min graph

Builds a (w, k)-minimizer index or a (k, s)-syncmer index of the threads in the GBWT
index. The graph can be any HandleGraph, which will be transformed into a GBWTGraph.
The transformation can be avoided by providing a GBWTGraph or a GBZ graph.

Required options:
    -d, --distance-index X  annotate the hits with positions in this distance index
    -o, --output-name X     store the index to file X

Minimizer options:
    -k, --kmer-length N     length of the kmers in the index (default 29, max 31)
    -w, --window-length N   choose the minimizer from a window of N kmers (default 11)
    -c, --closed-syncmers   index closed syncmers instead of minimizers
    -s, --smer-length N     use smers of length N in closed syncmers (default 18)

Weighted minimizers:
    -W, --weighted          use weighted minimizers
        --threshold N       downweight kmers with more than N hits (default 500)
        --iterations N      downweight frequent kmers by N iterations (default 3)
        --fast-counting     use the fast kmer counting algorithm (default)
        --save-memory       use the space-efficient kmer counting algorithm
        --hash-table N      use 2^N-cell hash tables for kmer counting (default: guess)

Other options:
    -l, --load-index X      load the index from file X and insert the new kmers into it
                            (overrides minimizer / weighted minimizer options)
    -g, --gbwt-name X       use the GBWT index in file X (required with a non-GBZ graph)
    -p, --progress          show progress information
    -t, --threads N         use N threads for index construction (default 8)
                            (using more than 16 threads rarely helps)
        --no-dist           build the index without distance index annotations (not recommended)


```


## mod


```
usage: vg mod [options] <graph.vg> >[mod.vg]
Modifies graph, outputs modified on stdout.

options:
    -P, --label-paths       don't edit with -i alignments, just use them for labeling the graph
    -c, --compact-ids       should we sort and compact the id space? (default false)
    -b, --break-cycles      use an approximate topological sort to break cycles in the graph
    -n, --normalize         normalize the graph so that edges are always non-redundant
                            (nodes have unique starting and ending bases relative to neighbors,
                            and edges that do not introduce new paths are removed and neighboring
                            nodes are merged)
    -U, --until-normal N    iterate normalization until convergence, or at most N times
    -z, --nomerge-pre STR   do not let normalize (-n, -U) zip up any pair of nodes that both belong to path with prefix STR
    -E, --unreverse-edges   flip doubly-reversing edges so that they are represented on the
                            forward strand of the graph
    -s, --simplify          remove redundancy from the graph that will not change its path space
    -d, --dagify-step N     copy strongly connected components of the graph N times, forwarding
                            edges from old to new copies to convert the graph into a DAG
    -w, --dagify-to N       copy strongly connected components of the graph forwarding
                            edges from old to new copies to convert the graph into a DAG
                            until the shortest path through each SCC is N bases long
    -L, --dagify-len-max N  stop a dagification step if the unrolling component has this much sequence
    -f, --unfold N          represent inversions accessible up to N from the forward
                            component of the graph
    -O, --orient-forward    orient the nodes in the graph forward
    -N, --remove-non-path   keep only nodes and edges which are part of paths
    -A, --remove-path       keep only nodes and edges which are not part of any path
    -k, --keep-path NAME    keep only nodes and edges in the path
    -R, --remove-null       removes nodes that have no sequence, forwarding their edges
    -g, --subgraph ID       gets the subgraph rooted at node ID, multiple allowed
    -x, --context N         steps the subgraph out by N steps (default: 1)
    -p, --prune-complex     remove nodes that are reached by paths of --length which
                            cross more than --edge-max edges
    -S, --prune-subgraphs   remove subgraphs which are shorter than --length
    -l, --length N          for pruning complex regions and short subgraphs
    -X, --chop N            chop nodes in the graph so they are not more than N bp long
    -u, --unchop            where two nodes are only connected to each other and by one edge
                            replace the pair with a single node that is the concatenation of their labels
    -e, --edge-max N        only consider paths which make edge choices at <= this many points
    -M, --max-degree N      unlink nodes that have edge degree greater than N
    -m, --markers           join all head and tails nodes to marker nodes
                            ('###' starts and '$$$' ends) of --length, for debugging
    -y, --destroy-node ID   remove node with given id
    -a, --cactus            convert to cactus graph representation
    -v, --sample-vcf FILE   for a graph with allele paths, compute the sample graph from the given VCF
    -G, --sample-graph FILE subset an augmented graph to a sample graph using a Locus file
    -t, --threads N         for tasks that can be done in parallel, use this many threads

```


## mpmap


```
usage: vg mpmap [options] -x graph.xg -g index.gcsa [-f reads1.fq [-f reads2.fq] | -G reads.gam] > aln.gamp
Multipath align reads to a graph.

basic options:
graph/index:
  -x, --graph-name FILE     graph (required; XG format recommended but other formats are valid, see `vg convert`) 
  -g, --gcsa-name FILE      use this GCSA2/LCP index pair for MEMs (required; both FILE and FILE.lcp, see `vg index`)
  -d, --dist-name FILE      use this snarl distance index for clustering (recommended, see `vg index`)
  -s, --snarls FILE         align to alternate paths in these snarls (unnecessary if providing -d, see `vg snarls`)
input:
  -f, --fastq FILE          input FASTQ (possibly gzipped), can be given twice for paired ends (for stdin use -)
  -i, --interleaved         input contains interleaved paired ends
algorithm presets:
  -n, --nt-type TYPE        sequence type preset: 'DNA' for genomic data, 'RNA' for transcriptomic data [RNA]
  -l, --read-length TYPE    read length preset: 'very-short', 'short', or 'long' (approx. <50bp, 50-500bp, and >500bp) [short]
  -e, --error-rate TYPE     error rate preset: 'low' or 'high' (approx. PHRED >20 and <20) [low]
output:
  -F, --output-fmt TYPE     format to output alignments in: 'GAMP for' multipath alignments, 'GAM' or 'GAF' for single-path
                            alignments, 'SAM', 'BAM', or 'CRAM' for linear reference alignments (may also require -S) [GAMP]
  -S, --ref-paths FILE      paths in the graph either 1) one per line in a text file, or 2) in an HTSlib .dict, to treat as
                            reference sequences for HTSlib formats (see -F) [all paths]
  -N, --sample NAME         add this sample name to output
  -R, --read-group NAME     add this read group to output
  -p, --suppress-progress   do not report progress to stderr
computational parameters:
  -t, --threads INT         number of compute threads to use [all available]

advanced options:
algorithm:
  -X, --not-spliced         do not form spliced alignments, even if aligning with --nt-type 'rna'
  -M, --max-multimaps INT   report (up to) this many mappings per read [10 rna / 1 dna]
  -a, --agglomerate-alns    combine separate multipath alignments into one (possibly disconnected) alignment
  -r, --intron-distr FILE   intron length distribution (from scripts/intron_length_distribution.py)
  -Q, --mq-max INT          cap mapping quality estimates at this much [60]
  -b, --frag-sample INT     look for this many unambiguous mappings to estimate the fragment length distribution [1000]
  -I, --frag-mean FLOAT     mean for a pre-determined fragment length distribution (also requires -D)
  -D, --frag-stddev FLOAT   standard deviation for a pre-determined fragment length distribution (also requires -I)
  -G, --gam-input FILE      input GAM (for stdin, use -)
  -u, --map-attempts INT    perform (up to) this many mappings per read (0 for no limit) [24 paired / 64 unpaired]
  -c, --hit-max INT         use at most this many hits for any match seeds (0 for no limit) [1024 DNA / 100 RNA]
scoring:
  -A, --no-qual-adjust      do not perform base quality adjusted alignments even when base qualities are available
  -q, --match INT           use this match score [1]
  -z, --mismatch INT        use this mismatch penalty [4 low error, 1 high error]
  -o, --gap-open INT        use this gap open penalty [6 low error, 1 high error]
  -y, --gap-extend INT      use this gap extension penalty [1]
  -L, --full-l-bonus INT    add this score to alignments that align each end of the read [mismatch+1 short, 0 long]
  -w, --score-matrix FILE   read a 4x4 integer substitution scoring matrix from a file (in the order ACGT)
  -m, --remove-bonuses      remove full length alignment bonuses in reported scores

```


## pack



Convert alignments to a compact coverage index.





```
usage: vg pack [options]
options:
    -x, --xg FILE          use this basis graph (any format accepted, does not have to be xg)
    -o, --packs-out FILE   write compressed coverage packs to this output file
    -i, --packs-in FILE    begin by summing coverage packs from each provided FILE
    -g, --gam FILE         read alignments from this GAM file (could be '-' for stdin)
    -a, --gaf FILE         read alignments from this GAF file (could be '-' for stdin)
    -d, --as-table         write table on stdout representing packs
    -D, --as-edge-table    write table on stdout representing edge coverage
    -u, --as-qual-table    write table on stdout representing average node mapqs
    -e, --with-edits       record and write edits rather than only recording graph-matching coverage
    -b, --bin-size N       number of sequence bases per CSA bin [default: inf]
    -n, --node ID          write table for only specified node(s)
    -N, --node-list FILE   a white space or line delimited list of nodes to collect
    -Q, --min-mapq N       ignore reads with MAPQ < N and positions with base quality < N [default: 0]
    -c, --expected-cov N   expected coverage.  used only for memory tuning [default : 128]
    -s, --trim-ends N      ignore the first and last N bases of each read
    -t, --threads N        use N threads (defaults to numCPUs)

```


## paths



Traverse paths in the graph.





```
usage: vg paths [options]
options:
  input:
    -x, --xg FILE            use the paths and haplotypes in this graph FILE. Supports GBZ haplotypes.
                             (Also accepts -v, --vg)
    -g, --gbwt FILE          use the threads in the GBWT index in FILE
                             (graph also required for most output options; -g takes priority over -x)
  output graph (.vg format)
    -V, --extract-vg         output a path-only graph covering the selected paths
    -d, --drop-paths         output a graph with the selected paths removed
    -r, --retain-paths       output a graph with only the selected paths retained
    -n, --normalize-paths    output a graph where all equivalent paths in a site a merged (using selected paths to snap to if possible)
  output path data:
    -X, --extract-gam        print (as GAM alignments) the stored paths in the graph
    -A, --extract-gaf        print (as GAF alignments) the stored paths in the graph
    -L, --list               print (as a list of names, one per line) the path (or thread) names
    -E, --lengths            print a list of path names (as with -L) but paired with their lengths
    -M, --metadata           print a table of path names and their metadata
    -C, --cyclicity          print a list of path names (as with -L) but paired with flag denoting the cyclicity
    -F, --extract-fasta      print the paths in FASTA format
    -c, --coverage           print the coverage stats for selected paths (not including cylces)
  path selection:
    -p, --paths-file FILE    select the paths named in a file (one per line)
    -Q, --paths-by STR       select the paths with the given name prefix
    -S, --sample STR         select the haplotypes or reference paths for this sample
    -a, --variant-paths      select the variant paths added by 'vg construct -a'
    -G, --generic-paths      select the generic, non-reference, non-haplotype paths
    -R, --reference-paths    select the reference paths
    -H, --haplotype-paths    select the haplotype paths paths
  configuration:
    -o, --overlay            apply a ReferencePathOverlayHelper to the graph
    -t, --threads N          number of threads to use [all available]. applies only to snarl finding within -n

```


## prune


```
usage: vg prune [options] <graph.vg> >[output.vg]

Prunes the complex regions of the graph for GCSA2 indexing. Pruning the graph
removes embedded paths.

Pruning parameters:
    -k, --kmer-length N    kmer length used for pruning
                           defaults: 24 with -P; 24 with -r; 24 with -u
    -e, --edge-max N       remove the edges on kmers making > N edge choices
                           defaults: 3 with -P; 3 with -r; 3 with -u
    -s, --subgraph-min N   remove subgraphs of < N bases
                           defaults: 33 with -P; 33 with -r; 33 with -u
    -M, --max-degree N     if N > 0, remove nodes with degree > N before pruning
                           defaults: 0 with -P; 0 with -r; 0 with -u

Pruning modes (-P, -r, and -u are mutually exclusive):
    -P, --prune            simply prune the graph (default)
    -r, --restore-paths    restore the edges on non-alt paths
    -u, --unfold-paths     unfold non-alt paths and GBWT threads
    -v, --verify-paths     verify that the paths exist after pruning
                           (potentially very slow)

Unfolding options:
    -g, --gbwt-name FILE   unfold the threads from this GBWT index
    -m, --mapping FILE     store the node mapping for duplicates in this file (required with -u)
    -a, --append-mapping   append to the existing node mapping

Other options:
    -p, --progress         show progress
    -t, --threads N        use N threads (default: 8)
    -d, --dry-run          determine the validity of the combination of options


```


## rna


```

usage: vg rna [options] graph.[vg|pg|hg|gbz] > splicing_graph.[vg|pg|hg]

General options:
    -t, --threads INT          number of compute threads to use [1]
    -p, --progress             show progress
    -h, --help                 print help message

Input options:
    -n, --transcripts FILE     transcript file(s) in gtf/gff format; may repeat
    -m, --introns FILE         intron file(s) in bed format; may repeat
    -y, --feature-type NAME    parse only this feature type in the gtf/gff (parses all if empty) [exon]
    -s, --transcript-tag NAME  use this attribute tag in the gtf/gff file(s) as id [transcript_id]
    -l, --haplotypes FILE      project transcripts onto haplotypes in GBWT index file
    -z, --gbz-format           input graph is in GBZ format (contains both a graph and haplotypes (GBWT index))

Construction options:
    -j, --use-hap-ref          use haplotype paths in GBWT index as reference sequences (disables projection)
    -e, --proj-embed-paths     project transcripts onto embedded haplotype paths
    -c, --path-collapse TYPE   collapse identical transcript paths across no|haplotype|all paths [haplotype]
    -k, --max-node-length INT  chop nodes longer than maximum node length (0 disables chopping) [0]
    -d, --remove-non-gene      remove intergenic and intronic regions (deletes all paths in the graph)
    -o, --do-not-sort          do not topological sort and compact the graph
    -r, --add-ref-paths        add reference transcripts as embedded paths in the graph
    -a, --add-hap-paths        add projected transcripts as embedded paths in the graph

Output options:
    -b, --write-gbwt FILE      write pantranscriptome transcript paths as GBWT index file
    -v, --write-hap-gbwt FILE  write input haplotypes as a GBWT with node IDs matching the output graph
    -f, --write-fasta FILE     write pantranscriptome transcript sequences as fasta file
    -i, --write-info FILE      write pantranscriptome transcript info table as tsv file
    -q, --out-exclude-ref      exclude reference transcripts from pantranscriptome output
    -g, --gbwt-bidirectional   use bidirectional paths in GBWT index construction


```


## sim


```
usage: vg sim [options]
Samples sequences from the xg-indexed graph.

basic options:
    -x, --xg-name FILE          use the graph in FILE (required)
    -n, --num-reads N           simulate N reads or read pairs
    -l, --read-length N         simulate reads of length N
    -r, --progress              show progress information
output options:
    -a, --align-out             write alignments in GAM-format
    -J, --json-out              write alignments in json
    --multi-position            annotate alignments with multiple reference positions
simulation parameters:
    -F, --fastq FILE            match the error profile of NGS reads in FILE, repeat for paired reads (ignores -l,-f)
    -I, --interleaved           reads in FASTQ (-F) are interleaved read pairs
    -s, --random-seed N         use this specific seed for the PRNG
    -e, --sub-rate FLOAT        base substitution rate (default 0.0)
    -i, --indel-rate FLOAT      indel rate (default 0.0)
    -d, --indel-err-prop FLOAT  proportion of trained errors from -F that are indels (default 0.01)
    -S, --scale-err FLOAT       scale trained error probabilities from -F by this much (default 1.0)
    -f, --forward-only          don't simulate from the reverse strand
    -p, --frag-len N            make paired end reads with given fragment length N
    -v, --frag-std-dev FLOAT    use this standard deviation for fragment length estimation
    -N, --allow-Ns              allow reads to be sampled from the graph with Ns in them
    --max-tries N               attempt sampling operations up to N times before giving up [100]
    -t, --threads               number of compute threads (only when using FASTQ with -F) [1]
simulate from paths:
    -P, --path PATH             simulate from this path (may repeat; cannot also give -T)
    -A, --any-path              simulate from any path (overrides -P)
    -m, --sample-name NAME      simulate from this sample (may repeat; requires -g)
    -R, --ploidy-regex RULES    use the given comma-separated list of colon-delimited REGEX:PLOIDY rules to assign
                                ploidies to contigs not visited by the selected samples, or to all contigs simulated
                                from if no samples are used. Unmatched contigs get ploidy 2.
    -g, --gbwt-name FILE        use samples from this GBWT index
    -T, --tx-expr-file FILE     simulate from an expression profile formatted as RSEM output (cannot also give -P)
    -H, --haplo-tx-file FILE    transcript origin info table from vg rna -i (required for -T on haplotype transcripts)
    -u, --unsheared             sample from unsheared fragments
    -E, --path-pos-file FILE    output a TSV with sampled position on path of each read (requires -F)

```


## stats


```
usage: vg stats [options] [<graph file>]
options:
    -z, --size             size of graph
    -N, --node-count       number of nodes in graph
    -E, --edge-count       number of edges in graph
    -l, --length           length of sequences in graph
    -L, --self-loops       number of self-loops
    -s, --subgraphs        describe subgraphs of graph
    -H, --heads            list the head nodes of the graph
    -T, --tails            list the tail nodes of the graph
    -e, --nondeterm        list the nondeterministic edge sets
    -c, --components       print the strongly connected components of the graph
    -A, --is-acyclic       print if the graph is acyclic or not
    -n, --node ID          consider node with the given id
    -d, --to-head          show distance to head for each provided node
    -t, --to-tail          show distance to head for each provided node
    -a, --alignments FILE  compute stats for reads aligned to the graph
    -r, --node-id-range    X:Y where X and Y are the smallest and largest node id in the graph, respectively
    -o, --overlap PATH    for each overlapping path mapping in the graph write a table:
                              PATH, other_path, rank1, rank2
                          multiple allowed; limit comparison to those provided
    -O, --overlap-all     print overlap table for the cartesian product of paths
    -R, --snarls          print statistics for each snarl
        --snarl-contents  print out a table of <snarl, depth, parent, contained node ids>
    -C, --chains          print statistics for each chain
    -F, --format          graph format from {VG-Protobuf, PackedGraph, HashGraph, XG}. Can't detect Protobuf if graph read from stdin
    -D, --degree-dist     print degree distribution of the graph.
    -b, --dist-snarls FILE print the sizes and depths of the snarls in a given distance index.
    -p, --threads N       number of threads to use [all available]
    -v, --verbose         output longer reports

```


## surject


```
usage: vg surject [options] <aln.gam> >[proj.cram]
Transforms alignments to be relative to particular paths.

options:
  -x, --xg-name FILE       use this graph or xg index (required)
  -t, --threads N          number of threads to use
  -p, --into-path NAME     surject into this path or its subpaths (many allowed, default: reference, then non-alt generic)
  -F, --into-paths FILE    surject into path names listed in HTSlib sequence dictionary or path list FILE
  -i, --interleaved        GAM is interleaved paired-ended, so when outputting HTS formats, pair reads
  -M, --multimap           include secondary alignments to all overlapping paths instead of just primary
  -G, --gaf-input          input file is GAF instead of GAM
  -m, --gamp-input         input file is GAMP instead of GAM
  -c, --cram-output        write CRAM to stdout
  -b, --bam-output         write BAM to stdout
  -s, --sam-output         write SAM to stdout
  -l, --subpath-local      let the multipath mapping surjection produce local (rather than global) alignments
  -T, --max-tail-len N     only align up to N bases of read tails (default: 10000)
  -P, --prune-low-cplx     prune short and low complexity anchors during realignment
  -a, --max-anchors N      use no more than N anchors per target path (default: unlimited)
  -S, --spliced            interpret long deletions against paths as spliced alignments
  -A, --qual-adj           adjust scoring for base qualities, if they are available
  -N, --sample NAME        set this sample name for all reads
  -R, --read-group NAME    set this read group for all reads
  -f, --max-frag-len N     reads with fragment lengths greater than N will not be marked properly paired in SAM/BAM/CRAM
  -L, --list-all-paths     annotate SAM records with a list of all attempted re-alignments to paths in SS tag
  -C, --compression N      level for compression [0-9]
  -V, --no-validate        skip checking whether alignments plausibly are against the provided graph
  -w, --watchdog-timeout N warn when reads take more than the given number of seconds to surject

```


## view



format conversions for graphs and alignments





```
usage: vg view [options] [ <graph.vg> | <graph.json> | <aln.gam> | <read1.fq> [<read2.fq>] ]
options:
    -g, --gfa                  output GFA format (default)
    -F, --gfa-in               input GFA format, reducing overlaps if they occur
    -v, --vg                   output VG format [DEPRECATED, use vg convert instead]
    -V, --vg-in                input VG format only
    -j, --json                 output JSON format
    -J, --json-in              input JSON format
    -c, --json-stream          streaming conversion of a VG format graph in line delimited JSON format
                               (this cannot be loaded directly via -J)
    -G, --gam                  output GAM format (vg alignment format: Graph Alignment/Map)
    -Z, --translation-in       input is a graph translation description
    -t, --turtle               output RDF/turtle format (can not be loaded by VG)
    -T, --turtle-in            input turtle format.
    -r, --rdf_base_uri         set base uri for the RDF output
    -a, --align-in             input GAM format
    -A, --aln-graph GAM        add alignments from GAM to the graph
    -q, --locus-in             input stream is Locus format
    -z, --locus-out            output stream Locus format
    -Q, --loci FILE            input is Locus format for use by dot output
    -d, --dot                  output dot format
    -S, --simple-dot           simplify the dot output; remove node labels, simplify alignments
    -u, --noseq-dot            shows size information instead of sequence in the dot output
    -e, --ascii-labels         use labels for paths or superbubbles with char/colors rather than emoji
    -Y, --ultra-label          label nodes with emoji/colors that correspond to ultrabubbles
    -m, --skip-missing         skip mappings to nodes not in the graph when drawing alignments
    -C, --color                color nodes that are not in the reference path (DOT OUTPUT ONLY)
    -p, --show-paths           show paths in dot output
    -w, --walk-paths           add labeled edges to represent paths in dot output
    -n, --annotate-paths       add labels to normal edges to represent paths in dot output
    -M, --show-mappings        with -p print the mappings in each path in JSON
    -I, --invert-ports         invert the edge ports in dot so that ne->nw is reversed
    -s, --random-seed N        use this seed when assigning path symbols in dot output
    -b, --bam                  input BAM or other htslib-parseable alignments
    -f, --fastq-in             input fastq (output defaults to GAM). Takes two 
                               positional file arguments if paired
    -X, --fastq-out            output fastq (input defaults to GAM)
    -i, --interleaved          fastq is interleaved paired-ended
    -L, --pileup               output VG Pileup format
    -l, --pileup-in            input VG Pileup format
    -B, --distance-in          input distance index
    -R, --snarl-in             input VG Snarl format
    -E, --snarl-traversal-in   input VG SnarlTraversal format
    -K, --multipath-in         input VG MultipathAlignment format (GAMP)
    -k, --multipath            output VG MultipathAlignment format (GAMP)
    -D, --expect-duplicates    don't warn if encountering the same node or edge multiple times
    -x, --extract-tag TAG      extract and concatenate messages with the given tag
    --verbose                  explain the file being read with --extract-tag
    --threads N                for parallel operations use this many threads [1]

```



# Haplotype sampling example

This directory contains a small test case for haplotype sampling using the `vg haplotypes` subcommand.
It is primarily intended for automated testing.
This readme contains instructions for replicating the tests manually.

Graph `micb-kir3dl1.gfa` consists of two subgraphs of the HPRC Minigraph-Cactus v1.1 graph:

* GRCh38#chr6:31498145-31511124 (micb)
* GRCh38#chr19:54816468-54830778 (kir3dl1)

File `HG003.fq.gz` contains single-end Illumina NovaSeq reads for HG003 mapping to the subgraph.
Read coverage is approximately 50x.
File `HG003.kff` contains 29-mer counts for the reads computed using KMC.

## Mapping reads to the original graph

First we need to build the indexes used by Giraffe:

```
vg autoindex --workflow giraffe --prefix micb-kir3dl1 -g micb-kir3dl1.gfa
```

This converts the graph into GBZ format (`micb-kir3dl1.giraffe.gbz`) and builds the distance index (`micb-kir3dl1.dist`) and the minimizer index (`micb-kir3dl1.min`).

Then we can map the reads, specifying the GBZ file, the fastq file, and the output file:

```
vg giraffe -Z micb-kir3dl1.giraffe.gbz -f HG003.fq.gz -p > to-original.gam
```

Option `-p` prints progress information.
We could also use option `-o bam` to project the alignments to reference paths and output them in BAM format.
However, because the graph contains two sets of reference paths (GRCh38 and CHM13), we would also have to specify which reference paths to use with option `--ref-paths`.

We can use the `vg stats` subcommand to compute statistics on the aligned reads:

```
vg stats -a to-original.gam micb-kir3dl1.giraffe.gbz
```

## Haplotype sampling

### Preprocessing the graph

In order to build the haplotype information used for haplotype sampling, we need a GBZ graph, a distance index, and an r-index.
We build them using the following commands:

```
vg gbwt -G micb-kir3dl1.gfa --max-node 32 --gbz-format -g micb-kir3dl1.gbz -r micb-kir3dl1.ri
vg index -j micb-kir3dl1.dist micb-kir3dl1.gbz
```

We used option `--max-node 32` to chop long nodes to 32 bp pieces instead of the default 1024 bp.
This makes the resulting GBZ graph (and hence the distance index) indentical to the one built by `vg autoindex` in the previous section.

Now we can build the haplotype information file `micb-kir3dl1.hapl`:

```
vg haplotypes -v 2 --subchain-length 300 -H micb-kir3dl1.hapl micb-kir3dl1.gbz
```

Option `-v 2` prints extended progress information, while `--subchain-length 300` sets the target length of subchains (blocks) to 300 bp.
The latter makes the example more interesting than the default 10000 bp.

### Manual sampling and read mapping

We can now sample the haplotypes using the haplotype information, the k-mer counts, and the original graph:

```
vg haplotypes -v 2 -i micb-kir3dl1.hapl -k HG003.kff \
    --include-reference --diploid-sampling \
    -g HG003.gbz micb-kir3dl1.gbz
```

Option `--include-reference` includes reference paths in the sampled graph, while option `--diploid-sampling` enables the diploid sampling mode.

Now we have the sampled graph `HG003.gbz`.
We could build the indexes for Giraffe manually.
But because the graph is a temporary object intended only for mapping a specific set of reads, we can let Giraffe handle index construction:

```
# Set the temporary directory as you wish
export TMPDIR=/scratch/tmp

vg giraffe -Z HG003.gbz -f HG003.fq.gz -p --index-basename ${TMPDIR}/HG003 > to-sampled.gam
```

Here we specify that the indexes go to the temporary directory with base name `HG003` using option `--index-basename`.

We can again use `vg stats` to compute statistics:

```
vg stats -a to-sampled.gam HG003.gbz
```

### Giraffe integration

Instead of sampling the haplotypes manually, we can let Giraffe handle everything:

```
# Set the temporary directory as you wish
export TMPDIR=/scratch/tmp

vg giraffe -Z micb-kir3dl1.gbz --haplotype-name micb-kir3dl1.hapl --kff-name HG003.kff \
    -f HG003.fq.gz -p --index-basename ${TMPDIR}/integration > giraffe-integration.gam
```

This uses the best practices for haplotype sampling, which currently means including the reference paths and using diploid sampling.
We put the sampled graph and the automatically built indexes into the temporary directory.
The naming scheme (e.g. `${TMPDIR}/integration.HG003.gbz`) consists of three parts:

* Base name `${TMPDIR}/integration` specified with option `--index-basename`.
* Sample name `HG003` guessed from the name of the KFF file or specified with option `--sample`.
* Extension (e.g. `gbz`) that depends on the file.

The statistics should be similar to those in the previous section:

```
vg stats -a giraffe-integration.gam ${TMPDIR}/integration.HG003.gbz
```

## Mapping reads to a linear reference

Given a GBZ file, we can extract the paths corresponding to a given sample in fasta format:

```
vg paths -x micb-kir3dl1.gbz --extract-fasta --sample GRCh38 > reference.fa
```

The above example extracts the relevant parts of GRCh38.
The resulting fasta file can be used as the reference with any linear aligner.

## Further reading

Documentation for the `vg haplotypes` subcommand can be found in the [vg wiki](https://github.com/vgteam/vg/wiki/Haplotype-Sampling).
You may also want to read about the [best practices](https://github.com/vgteam/vg/wiki/Giraffe-best-practices) for mapping reads with Giraffe.

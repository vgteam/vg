# SAFARI: mapping ancient DNA with `vg safari`

`vg safari` maps ancient DNA reads to a pangenome graph. It is `vg giraffe` with one more way of
finding seeds, made for reads that deamination has damaged. It was first written as a fork of vg
1.44 (<https://github.com/grenaud/SAFARI>) and is now part of vg.

Please cite:

> Joshua Rubin, Jan van Waaij, Louis Kraft, Jouni Sirén, Peter Wad Sackett, Gabriel Renaud.
> SAFARI: pangenome alignment of ancient DNA using purine/pyrimidine encodings.
> *NAR Genomics and Bioinformatics*, Volume 7, Issue 4, December 2025, lqaf170.
> <https://doi.org/10.1093/nargab/lqaf170>

```bibtex
@article{rubin2025safari,
  title   = {{SAFARI}: pangenome alignment of ancient {DNA} using purine/pyrimidine encodings},
  author  = {Rubin, Joshua and van Waaij, Jan and Kraft, Louis and Sir{\'e}n, Jouni and
             Sackett, Peter Wad and Renaud, Gabriel},
  journal = {NAR Genomics and Bioinformatics},
  volume  = {7},
  number  = {4},
  pages   = {lqaf170},
  year    = {2025},
  month   = dec,
  doi     = {10.1093/nargab/lqaf170},
  url     = {https://academic.oup.com/nargab/article/7/4/lqaf170/8376687}
}
```

## The problem

Giraffe seeds a read with its minimizers: k-mers of the read that it looks up exactly in a
minimizer index of the graph. In ancient DNA, cytosine deaminates, so a C in the molecule reads as
T, and on the other strand a G reads as A. The damage is heaviest near the ends of the fragment,
and ancient fragments are short, so a damaged base often falls inside every minimizer of a read,
and the read gets no seed at its true place.

## RYmers

An RYmer is a k-mer written with R for a purine (A or G) and Y for a pyrimidine (C or T). C to T
and G to A each keep a base in its class, so deamination never changes an RYmer.

`vg rymer` builds the RYmer index from a GBZ graph. It chooses minimizers from the DNA of the
haplotypes exactly as `vg minimizer` would with the same k and w, then stores each under its RYmer,
with its graph position and its DNA k-mer. Looking up an RYmer therefore gives every graph k-mer
that a read k-mer with that RYmer could be a damaged copy of.

## Seeding

For each read, `vg safari`:

1. finds the read's minimizers as Giraffe does;
2. finds the read's minimizers under the RYmer index's k and w, and looks up each one's RYmer;
3. for each graph k-mer found, computes the posterior probability that the read k-mer is that
   graph k-mer with damage, rather than an unrelated k-mer with the same RYmer (below);
4. for the first graph k-mer of a minimizer whose posterior is above `--posterior-threshold`,
   makes a stand-in read: all N except for that graph k-mer at the read k-mer's offset;
5. finds minimizers in each stand-in read with the ordinary minimizer index and adds them to the
   read's own, as seeds at the same read offsets.

Everything after seeding (clustering, extension, alignment, mapping quality) is Giraffe's.

### The posterior

Two models explain a read k-mer given a graph k-mer with the same RYmer:

- **Damage.** Each read base comes from the graph base at the same offset, with the probability the
  deamination profiles give. At a position, for each true base, the profile of the end (5' or 3')
  more likely to change that base is used. As in the original SAFARI, the position used for every
  base of the k-mer is that of its first base.
- **Spurious.** The k-mer comes from elsewhere. For k-mers shorter than 22 bases the number of
  mismatches is binomial with rate 1.0014 · k^−0.6628; longer k-mers have likelihood 0.01.

With `--spurious-prior` as the prior of the second model, the posterior of the first is reported.
Equal k-mers have posterior 1.

## Deamination profiles

`--deam-5p` and `--deam-3p` take `.prof` tables: a header, then
one row per position from that end, outermost first, with tab-separated rates for
`A>C A>G A>T C>A C>G C>T G>A G>C G>T T>A T>C T>G`. Positions past the last row use the last row.
The files may be gzip-compressed. The original SAFARI repository has example profiles in `profs/`.

## Usage

```bash
vg autoindex --workflow giraffe -r ref.fa -v variants.vcf.gz -p graph   # or any Giraffe indexes
vg rymer -k 29 -w 11 -o graph.ry graph.giraffe.gbz
vg safari -Z graph.giraffe.gbz -d graph.dist -m graph.shortread.withzip.min -z graph.shortread.zipcodes \
    -q graph.ry --deam-5p sample_5p.prof --deam-3p sample_3p.prof -f merged_reads.fq.gz > out.gam
```

`vg safari` takes every `vg giraffe` option. Unless `--hard-hit-cap` or a parameter preset is
given, it sets the hard hit cap to 50000, the original SAFARI default. RYmer seeds are added to
single-end reads only; merge read pairs before mapping, as is usual for ancient DNA.

## Differences from the original SAFARI

- Mapping quality is Giraffe's. The original SAFARI set every mapped read's MAPQ to 60 before
  Giraffe's caps were applied.
- A minimizer found in a stand-in read is dropped if the read already has the same minimizer at
  the same offset, so an undamaged k-mer is not counted twice.
- The read k-mer's offset is the minimizer's own; the original searched the read for the first
  occurrence of its RYmer, which is the same offset unless the RYmer repeats earlier in the read.
- Damage rates are looked up without a precomputed table, so reads longer than 1000 bases work.
- Everything else in vg is current, not vg 1.44, and `vg giraffe` is unchanged.

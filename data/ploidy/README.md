# Ploidy files for `vg call`

Each file gives the ploidy of one sample sex on one reference assembly, for `vg call --ploidy-bed`.
Pass the same file to every run, whole-genome or one contig at a time: a contig with no line is
diploid, so the autosomes need none.

| File | Assembly | Contig names | Sample |
|---|---|---|---|
| `chm13v2.XY.ploidy.bed` | T2T-CHM13v2.0 | `chrX`, `chrY`, `chrM` | male |
| `chm13v2.XX.ploidy.bed` | T2T-CHM13v2.0 | `chrX`, `chrY`, `chrM` | female |
| `grch38.XY.ploidy.bed` | GRCh38 | `chrX`, `chrY`, `chrM` | male |
| `grch38.XX.ploidy.bed` | GRCh38 | `chrX`, `chrY`, `chrM` | female |
| `grch37.ucsc.XY.ploidy.bed` | GRCh37 (UCSC hg19) | `chrX`, `chrY`, `chrM` | male |
| `grch37.ucsc.XX.ploidy.bed` | GRCh37 (UCSC hg19) | `chrX`, `chrY`, `chrM` | female |
| `grch37.b37.XY.ploidy.bed` | GRCh37 (b37, hs37d5) | `X`, `Y`, `MT` | male |
| `grch37.b37.XX.ploidy.bed` | GRCh37 (b37, hs37d5) | `X`, `Y`, `MT` | female |

In an XY file, chrX is haploid outside its two pseudoautosomal regions (PARs), which are diploid,
and chrY and the mitochondrion are haploid. chrY is haploid along its whole length, PARs included,
since the sample's second PAR copy is counted on chrX. In an XX file, chrX is diploid, chrY has
ploidy 0, which means no calls, and the mitochondrion is haploid.

Each file's header cites the source of its contig lengths and PAR coordinates. The format is
`CHROM START END PLOIDY`, 0-based and half-open, with `CHROM` spelled as in the output VCF (`chrX`,
not `CHM13#0#chrX`). Ploidy is 0, 1 or 2, and intervals must not overlap. For another assembly or
karyotype, copy the nearest file and edit it.

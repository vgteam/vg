# Read-likelihood genotyping (`vg call --read-likelihood`)

`vg call --read-likelihood` genotypes one sample from its reads aligned to a pangenome graph.
Default values of the options named below are listed by `vg call --help`; the few constants that
have no option are listed under [Fixed constants](#fixed-constants).

The method is described in four documents, and its source code in a fifth. This one describes the
caller as a whole, defines the vocabulary that the others use, and can be read on its own.
[read-likelihood-direct-genotyping.md](read-likelihood-direct-genotyping.md) describes how a site is
genotyped from its reads, [read-likelihood-linkage-model.md](read-likelihood-linkage-model.md)
describes the model that re-decides genotypes from the haplotypes stored in the graph and phases
them, and [read-likelihood-read-phasing.md](read-likelihood-read-phasing.md) describes how the reads
then re-decide the phases and correct the likelihoods. Read this document first. Read each of the
other three after it, or when this document reaches its part of the method, under [Settling
genotypes](#settling-genotypes) and [Phasing](#phasing). The source files that implement the method,
and an order in which to read them, are described in
[read-likelihood-architecture.md](read-likelihood-architecture.md).

A minimal run takes a graph in
[GBZ](https://github.com/vgteam/vg/wiki/Extra-details-on-vg-file-formats#gbz-gbwtgraph--gbz) format
and reads aligned to it by `vg giraffe` in
[GAM](https://github.com/vgteam/vg/wiki/Extra-details-on-vg-file-formats#gam-graph-alignment--map-vgs-bam)
format, and writes the calls in [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format:

```
vg call graph.gbz --read-likelihood --gam reads.gam > calls.vcf
```

`vg call` finds the graph's sites itself. It can instead read them, with `--snarls`, from a file
of [snarls](https://github.com/vgteam/vg/wiki/Snarls-and-chains#snarls) made by `vg snarls`.

The caller works in four steps:

1. **Site likelihood computation.** For each genotype the site could have, vg computes a **site
   likelihood**: a number that measures how well the genotype explains the site's reads.
2. **Genotyping.** vg chooses the site's genotype. It starts from the site's **direct call**, the
   genotype with the highest site likelihood. Under the **linkage model**, it then re-decides the
   genotype from the site likelihoods of the site and its neighbours, using the haplotypes stored in
   the graph, which tend to carry the same combinations of alleles at neighbouring sites as the
   sample. With [nested calling](#nested-sites), on by default, sites nested inside other sites are
   genotyped too, and get VCF records of their own.
3. **Phasing.** Where the linkage model runs, vg assigns each genotype's alleles to the sample's
   haplotypes, first from the stored haplotypes and then, optionally, from reads that span several
   sites. With `--regenotype`, the phase is then used to correct the site likelihoods, and
   genotyping and phasing are repeated. The genotype of the last repetition is the site's **settled
   genotype**, the one vg reports.
4. **Output.** vg writes a VCF file. It can also write a **mosaic** file, which describes each of
   the sample's haplotypes as a walk through the graph, and an
   [anchor file](#assembly-anchors---anchors-out), which ties reads to haplotypes.

The caller is built in two parts. **Direct genotyping**, which is step 1 and the direct call,
genotypes each site on its own, from the site's reads. **Linkage-based genotyping**, the linkage
model, works on top of it: it chooses the genotypes of many neighbouring sites together, from
their site likelihoods, and gives the first phase in step 3. How these steps are ordered into
passes over the sites is described under [Passes and rounds](#passes-and-rounds), after the
vocabulary.

## Vocabulary

### Graph

- **Node.** A node holds a DNA sequence and can be traversed in either **orientation**: forward,
  reading its sequence, or reverse, reading its reverse complement. A **walk** is a sequence of
  oriented node visits along the graph's edges.
- **Site.** A place in the graph where the sample's genome may differ from the reference. vg
  identifies sites as **snarls**. A snarl is a subgraph separated from the rest of the graph by two
  **boundary nodes**, a start and an end. On this page a site is a snarl that `vg call` considers
  for genotyping. Its boundary nodes belong to it, and its other nodes, including those of any
  snarls nested in it, are its **interior**. A site is oriented from its start boundary node to its
  end boundary node, and its alleles are written in that direction.
- **Chain.** A series of snarls joined end to end, each snarl's end boundary node being the next
  one's start. Snarls nest: a snarl can contain chains of smaller snarls, its **child chains**. The
  sites in a site's child chains are **nested** in it, and it is their **parent**. A site with no
  parent is a **top-level** site. The linkage chains defined below, and the phase chains of read
  phasing, are instead sequences of sites.
- **Allele.** A walk through a site from its start boundary node to its end boundary node (in the
  code, a *traversal*).
- **Reference.** The **reference paths** are the graph paths on which `vg call` reports positions.
  By default they are the paths that the graph's metadata marks as reference paths or as generic
  paths (named paths with no sample); `--ref-path`, `--path-prefix` or `--ref-sample` selects
  other paths instead, by name, by name prefix or by sample. A site's **reference allele** is the
  walk a reference path takes through it. In the VCF the reference allele is allele 0, and the
  other alleles are numbered from 1. A site's **position** is the first base of whichever of its
  boundary nodes comes first on the reference path. `vg paths --compute-gref` adds a **gRef
  cover** (graph-reference cover) to a graph: a **gRef copy** of each reference path, under a
  name starting `gref_`, and **gRef fragments**, paths named `gref_<reference>_<N>_alt` that run
  through sequence the reference does not cover. Selecting the cover's paths as the reference
  paths, for example with `--ref-sample`, makes the fragments reference paths. The panel (below)
  leaves out the gRef fragments, and leaves out each gRef copy unless the original's sample is
  absent from the GBWT, so that the reference is in the panel once.

### Sample and reads

- **Haplotype.** One copy of a chromosome, or of part of one. A haplotype that passes through a
  site carries one allele there.
- **Ploidy.** A site's ploidy, $P$, is the number of the sample's haplotypes that pass through the
  site. Options set the ploidy of each contig or region, and a nested site's ploidy follows from its
  parent's genotype (see [Ploidy](#ploidy)). The ploidy options accept only 1 or 2, because the
  linkage model, phasing and nested calling are built for at most two haplotypes. A nested site that
  none of its parent's alleles passes through has ploidy 0 (see
  [Which child chains are genotyped](#which-child-chains-are-genotyped)).
- **Panel.** The haplotypes stored in the GBZ file's
  [GBWT](https://github.com/vgteam/vg/wiki/Extra-details-on-vg-file-formats#gbwt-gbwt) index, or in
  a separate GBWT file given with `--gbwt`. A panel haplotype can be stored as several paths that
  each cover part of a chromosome, so it need not pass through every site.
- **Genotype.** The multiset of the $P$ alleles that the sample's haplotypes carry at a site, such
  as $\lbrace 0, 1 \rbrace$.
- **Strand.** One of the sample's haplotypes at a site. At ploidy 2 the two strands are numbered 0
  and 1. At ploidy 1 the strand is 0, except at a site of a nested haploid chain (see
  [Phasing](#phasing)), which is on the strand that carries it. The word names a haplotype of the
  sample, not a strand of DNA. Which strand carries which allele of a genotype is the genotype's
  **phase**, and is decided separately from the genotype.
- **Read placements.** A mapper, such as `vg giraffe`, aligns the reads to the graph. The reads'
  **placements** are the walks their alignments take, together with the **edits** inside each node
  (the runs of matching bases, mismatches, insertions and deletions against the node's sequence)
  and their mapping qualities (MAPQ).
- **Reads of a site.** The reads whose placements visit an interior node of the site, or visit
  both of its boundary nodes. A read that visits both boundary nodes and no interior node supports
  an allele that skips the interior, such as a deletion. A read whose only node in the site is one
  boundary node fits every allele equally, and is not used.

## Passes and rounds

vg goes over the sites in **passes**. The **direct pass** comes first and runs once. It is the only
pass that fetches reads: the steps after it work from the evidence it keeps for each read.
**Rounds** follow, each built around a **linkage pass**. Last, the **render** builds the records,
and vg writes the output files:

```
direct pass                         once, every top-level site
    direct genotyping               the site's likelihoods and direct call, from its reads
    descent                         the same for the sites nested in it, at provisional ploidies
    staging                         keep what the site's records will be built from
round 1, 2, ...                     round 2 on only with --regenotype
    likelihood correction           from round 2 on: correct the likelihoods from the phase
    linkage pass                    one level at a time, level 0 first
        linkage-based genotyping    choose a genotype for each site of the level
        phasing from the panel      put each genotype's alleles on the strands
        ploidy of the next level    from the genotypes just chosen
    read phasing                    re-decide the phases from the reads (--read-phasing)
render                              build each site's records from its settled genotype
write                               the VCF, and the mosaic and anchor files if asked for
```

The order follows from what each step needs. The linkage model chooses a site's genotype from the
likelihoods of every site of the site's [linkage chain](#linkage-based-genotyping), so direct
genotyping has run on every site before any linkage-based genotyping runs. A nested site's ploidy
depends on its parent's genotype, so a linkage pass reaches a site only after its parent. It goes
by **level**: top-level sites are level 0, the sites of their child chains level 1, and so on (see
[Nested sites](#nested-sites)).

A site's **chosen genotype** is the one the latest linkage pass gave it: the genotype the linkage
model chose for it, or its direct call where the model leaves the site out or does not run. Each
linkage pass chooses afresh. The chosen genotype of the last round is the site's **settled
genotype**, the one vg reports (see [Settling genotypes](#settling-genotypes)).

### The direct pass

vg visits the top-level sites, and goes through the reads once. At each site it runs direct
genotyping, which computes the site's likelihoods and its direct call, and then descends into the
site's child chains and does the same for their sites. A nested site gets a **provisional ploidy**:
the number of its parent's direct-call alleles that cross it. Under the linkage model, a nested site
that none of them crosses is genotyped too, at the parent's ploidy, since a linkage pass may choose
for the parent a genotype whose alleles cross it; without the linkage model it is skipped. When the
nested site has more than one candidate allele, vg also computes its likelihoods and direct call at
ploidy 1 or 2, whichever it was not genotyped at, because a linkage pass may give it either. A
nested site with one candidate allele is genotyped at one ploidy only, and keeps it unless a
linkage pass drops the site. Each site is **staged**: vg keeps what the site's records will be
built from, and writes them in the render.

### The linkage pass

A linkage pass goes through the levels in order. At each level, the linkage model chooses a genotype
for every site of the level's linkage chains, and phases each chain from the panel (see
[Phasing](#phasing)). Each site of the next level then takes its ploidy from its parent's chosen
genotype, and vg uses its likelihoods and direct call at that ploidy: those the direct pass
computed, corrected from round 2 on (see [Rounds](#rounds)). A site at ploidy 0 is dropped, with
everything nested in it. At level 0 each linkage chain starts from the whole panel. At a deeper
level it starts from the panel haplotypes that the parent's strands copy.

Each linkage pass decides every ploidy afresh, so a site dropped by one linkage pass can come back
in the next. vg records which of a parent's candidate alleles cross a nested site only when the
parent has at most a fixed number of candidate alleles (see [Fixed constants](#fixed-constants)). A
nested site under a parent with more keeps its provisional ploidy, even where that disagrees with
its parent's chosen genotype. Without the linkage model, the linkage pass changes nothing: every
site keeps its direct call and its provisional ploidy.

### Rounds

Round 1 is a linkage pass followed by [read phasing](#from-the-reads), which runs when
`--read-phasing` is given. With [`--regenotype`](#re-genotyping-from-the-phase), more rounds
follow. Each starts with a likelihood correction: the phase that the previous round ended with
corrects the likelihoods that the direct pass computed, not those of the previous round. The
linkage pass then chooses the genotypes from the corrected likelihoods, and read phasing runs on
them. The new genotypes can change the phase, and with it the next correction.

`--regeno-passes` caps the number of rounds, and so of linkage passes, the first included. No round
follows one whose correction changes no site's direct call, whose linkage pass changes no chosen
genotype, or whose chosen genotypes return to those of an earlier round. Every correction that is
applied is followed by a linkage pass, so the settled genotypes are chosen from the likelihoods
that `GL` reports. With `--regeno-passes 1` there is one round. After it, the correction is
computed and reported, on standard error and in the `--regeno-ledger` file, but not applied, and no
linkage pass follows. `--regeno-ledger` writes one
line for each site whose direct call the last correction changed.

## Settling genotypes

A site's **settled genotype** is the genotype that vg reports for it. vg chooses it among the
genotypes that can be made from the site's [candidate alleles](#candidate-alleles) at the site's
[ploidy](#ploidy). The settled genotype is either the site's direct call or the genotype that the
linkage model chose in the last round:

- Without the linkage model, that is under support enumeration (see
  [Candidate alleles](#candidate-alleles)), with `--linkage-weight 0`, or with a panel of fewer
  than two haplotypes, every site's settled genotype is its direct call.
- Otherwise, under haplotype enumeration, the linkage model chooses the sites' genotypes.

### Direct genotype

Direct genotyping computes the site likelihood $\mathcal{L}(G)$ of every genotype $G$ that can be
made from the site's candidate alleles. It is described in
[read-likelihood-direct-genotyping.md](read-likelihood-direct-genotyping.md). The genotype with the
highest $\mathcal{L}(G)$ is the site's direct call.

A tie, as when every read fits the tied genotypes equally well, goes to the homozygous reference
genotype if it is one of the tied genotypes. Otherwise it goes to the tied genotype that comes
first in the following order. Number the candidate alleles in the order of the site's list of
candidate alleles (see [Candidate alleles](#candidate-alleles)), write each genotype as its allele
numbers in increasing order, and sort the genotypes by their largest number, then by the next
largest, and so on. This is the order in which VCF lists `GL` entries, applied to the candidate
alleles. The reference allele need not come first in that list, so this order can differ from
that of a record's `GL` entries, which follow the record's alleles (see [VCF fields](#vcf-fields)).

Direct genotyping also keeps two values for each read $r$ of the site, which phasing,
re-genotyping and the anchor file use. The read's **mismapping probability** $e_r$ is the
probability that the mapper put it in the wrong place, computed from its MAPQ. Its **relative
likelihood** $p_{ra}$ under each candidate allele $a$ is 1 for the allele that fits the read best,
and smaller for alleles that fit it worse.

### Linkage-based genotyping

The linkage model treats each strand as a mosaic of panel haplotypes, a Li–Stephens copying model
as in PanGenie (Ebler et al., Nature Genetics 2022). It is described in
[read-likelihood-linkage-model.md](read-likelihood-linkage-model.md).

The model runs along **linkage chains**. A linkage chain is a sequence of sites in order of
position. At the top level, a linkage chain holds the top-level sites of one contig, split wherever
the ploidy changes. Below the top level, a linkage chain holds the sites of one child chain that
have the same ploidy, whether or not they are adjacent; at ploidy 1, also on the same strand of the
parent. Strands do not correspond across a ploidy change at the top level, but below it they are the
parent's strands, which its phase names. The model's evidence at each site is the site's
likelihoods. From them, and from the alleles that the panel haplotypes carry at every site of the
chain, it computes each genotype's **posterior probability**, and chooses the genotype with the
highest. A site with more alleles than the model can hold is left out of it (see
[States and emissions](read-likelihood-linkage-model.md#states-and-emissions)), and keeps its direct
call, written unphased.

The linkage model also phases each linkage chain from the panel (see [Phasing](#phasing)).

### Sites with no reads

A site with no reads is left without a genotype, and the linkage model leaves it out. It gets a
record only when `--genotype-snarls` is given, and that record has a missing `GT` and `FILTER`
`noreads`. Every read of a nested site is also a read of its parent, so the sites nested in a site
with no reads have none either, and are not genotyped.

## Candidate alleles

A site's candidate alleles come from the panel or from the reads.

- When the graph is a GBZ with at least two panel haplotypes, the candidate alleles are by default
  the distinct walks that panel haplotypes take through the site (**haplotype enumeration**). A
  GBZ with fewer may hold only the reference, which would offer no other allele. `--gbz` asks for
  haplotype enumeration from the GBZ whatever its number of haplotypes, and `--gbwt` takes the
  panel from a separate GBWT file. So an allele that no panel haplotype takes cannot be called,
  unless it is the reference allele. A panel haplotype that passes through the site more than once
  offers each of its walks through it.
- Otherwise, or with `--enumerate-support`, the candidates are the walks with the most read
  support (**support enumeration**). vg finds them with Yen's k-shortest-paths algorithm, over the
  node and edge coverage in a file made by
  [`vg pack`](https://github.com/vgteam/vg/wiki/vg-manpage#pack) and given with `--pack`. vg keeps
  at most a fixed number of them (see [Fixed constants](#fixed-constants)).
- The reference allele is always a candidate.

vg keeps a site's candidate alleles in a list, in the order in which vg's search finds them, with
the reference allele added at the end if the search did not find it.

`--max-snarl-edges` sets a limit on a site's edges, counting the edges of the sites nested in it.
vg skips a site over the limit, and genotypes the sites of its child chains as top-level sites.
They take the ploidy of their contig or region, and join the contig's top-level sites in the
linkage model. A skipped site with no child chains gets no genotype and no record, even with
`--genotype-snarls`. The limit exists because Yen's search is slow on very large sites. Under
haplotype enumeration, which does not use Yen's search, there is no limit by default.

## Ploidy

A site's ploidy comes from `--ploidy`. `--ploidy-regex` overrides it per contig, with a
comma-separated list of `REGEX:PLOIDY` rules. For each reference path, the first rule whose regular
expression matches the path's name applies. `--ploidy-bed` overrides both per region, with a BED
file of `CHROM START END PLOIDY` lines. `CHROM` is spelled as in the output VCF. Intervals are
0-based and half-open, and must not overlap. A site takes the ploidy of the interval that contains
its position. A site outside every interval keeps the ploidy from `--ploidy` or `--ploidy-regex`.
Every ploidy must be 1 or 2. A site that these options make haploid has a single allele as its
`GT`.

A nested site takes its ploidy from its parent's genotype instead (see
[Which child chains are genotyped](#which-child-chains-are-genotyped)). The exception is a child of
a site that vg could not genotype as a top-level site. (This happens when the site has no reads,
when `--max-snarl-edges` or the allele-length limits `-c` and `-C` skip it, when no reference path
runs through it from one boundary node to the other, or when vg finds no allele through it.) vg
genotypes such a child as a top-level site, at the ploidy of its contig or region.

## Nested sites

A site can contain child chains, and their sites can contain chains in turn. **Nested calling**
genotypes the sites of each child chain and writes them in records of their own. In the direct
pass, vg genotypes the sites of a site's child chains by recursion from the site, which is called
**descent**. A site's **level** counts the descents that reach it.
A site that vg genotypes as a top-level site is level 0, including the child of a site that could
not be genotyped (see [Ploidy](#ploidy)). A site reached by descent is one level below its parent.
The level is not `INFO/LV`, which counts a record's enclosing sites that have records (see
[Nesting tags](#nesting-tags)).

Nested calling is on by default with `--read-likelihood`, and `--nested` turns it on with the other
genotyping methods of `vg call`. `--no-nested` turns it off. With `--no-nested`, each site is
genotyped against its full walks, variation inside nested sites is reported in the enclosing
site's alleles, and a nested site is genotyped on its own only when its parent could not be
genotyped. How nested sites are written as records is described under
[Records](#records).

### Which child chains are genotyped

A strand has a copy of a nested site when the strand's allele at the parent crosses the site. So a
nested site's ploidy is the number of the parent's strands whose allele crosses it, and a strand
whose allele crosses it twice counts once. Each linkage pass decides it from the parent's chosen
genotype (see [The linkage pass](#the-linkage-pass)). Ploidy is decided site by site, so two sites
of one child chain can have different ploidies.

- **0**: the sample has no copy of the site. vg writes no record for it or for the sites nested in
  it.
- **1**: the site is genotyped at ploidy 1. When the genotypes are phased and the site is in a
  [nested haploid chain](#phasing), its `GT` is written `a|.` or `.|a`, and the position of the
  allele says which strand carries the site. Otherwise the `GT` is a single allele.
- **2**: the site is genotyped at ploidy 2.

A child chain that its parent's reference allele does not cross is an **off-reference chain**. Its
sites have no position on the parent's reference path. vg skips an off-reference chain, except in
these cases:

- With `--anchors-out`, unless `--no-off-ref-nesting` is given. The chain's sites then have
  entries in the anchor file, and no VCF records. Phasing [from the reads](#from-the-reads)
  includes these sites, so genotyping these chains can change the phase written for other records
  and, under [`--regenotype`](#re-genotyping-from-the-phase), their genotypes.
- When the reference paths include a [gRef fragment](#graph). A chain on a gRef fragment then gets
  records, with the fragment as their contig.

## Phasing

Phasing decides each genotype's phase: which strand carries which allele. At ploidy 2, a site's
phase is an order of the two alleles of its chosen genotype. The first allele is on strand 0 and
is written to the left of the `|` in `GT`, and the second is on strand 1.

A diploid site is **phaseable** when its chosen genotype holds two different alleles, so
that its two possible phases differ. A phaseable site can still be homozygous in the VCF, when its
two alleles differ only inside a child chain.

A **nested haploid chain** is a linkage chain of ploidy-1 sites whose parent is diploid, or is a
site of another nested haploid chain. All its sites lie on one strand of the nearest diploid
ancestor, the strand that carries them, directly or through ploidy-1 parents (see
[Which child chains are genotyped](#which-child-chains-are-genotyped)). Its sites' `GT` is `a|.` on
strand 0 and `.|a` on strand 1.

The linkage model first phases each linkage chain from the panel. It finds the panel haplotypes
that the strands most probably copy along the chain, given each site's chosen genotype, and orders
each site's alleles to match, keeping the genotype. It is described in
[Phasing from the panel](read-likelihood-linkage-model.md#phasing-from-the-panel). Read phasing,
when on, then re-decides the phases from reads that span several sites (see
[From the reads](#from-the-reads)).

Each record's `FORMAT/PS` names its **phase set**: the reference position of the first site of its
top-level linkage chain. A nested site takes its parent's phase set. Sites of one phase set are
phased relative to one another. So a phase set spans a contig, or the part of one between ploidy
changes, and does not mark where the phase is reliable.

At some phaseable sites, both strands of the panel phase copy the
[wildcard](read-likelihood-linkage-model.md#wildcard-haplotype) or a panel haplotype that does
not pass through the site, so the panel does not order the two alleles. vg still writes them with
`|`, in an order that carries no phase, and keeps the site in its phase set, so that read phasing,
when on, can order it. No field marks these sites.

Phasing is on wherever the linkage model runs. `--phased` makes vg call fail when the linkage model
does not run (see [Settling genotypes](#settling-genotypes)), and `--no-phased` turns phasing off.
Read phasing, re-genotyping and `--anchors-hom-split` start from the panel phase, so an explicit
`--read-phasing`, `--regenotype` or `--anchors-hom-split` is an error with `--no-phased` or where
the linkage model does not run, and a preset's are turned off. Where the linkage model runs,
`--no-phased` also turns nested calling off (an explicit `--nested` is then an error), and variation
inside nested sites is reported in the enclosing site's alleles. Nested calling needs phasing there
because the linkage pass takes each nested site's strand, and the panel haplotypes its linkage chain
starts from, from its parent's phase.

### From the reads

A read that spans two phaseable sites shows directly whether their alleles lie on the same strand.
`--read-phasing` uses such reads to re-decide the phase of each phaseable site, nested and
off-reference sites included, within the phase sets the panel gave. It leaves out the chains that a
[block record](#block-records) spells out. Read phasing changes phases and keeps every genotype. It
is off by default and on under `--preset ont`. It uses its own statistics, not the linkage model's,
and is described in [Read phasing](read-likelihood-read-phasing.md#read-phasing).

### Re-genotyping from the phase

`--regenotype` uses the phase to give each read its own probabilities of having come from each
strand, and so corrects each site's likelihoods, at the start of every round after the first (see
[Rounds](#rounds)). Off-reference chains, and chains that a block record spells out, keep the
likelihoods of the direct pass. It needs `--read-phasing` and, like it, is off by default and on
under `--preset ont`. It cannot be combined with `--top-down` or `--bottom-up`: an explicit
`--regenotype` with either is an error, and a preset's is turned off. The correction is described
in [Re-genotyping from the phase](read-likelihood-read-phasing.md#re-genotyping-from-the-phase).

## Output

### Records

From the settled genotype vg writes the site's **records**, its lines in the VCF. A site usually
has one record. Under nested calling it can have one for each place where it differs from the
reference allele (below). A site has no record when every allele of its settled genotype is
written as the reference allele, unless `--genotype-snarls` is given.

#### Reporting each difference once

Two alleles of a site that take the same route except inside a child chain differ only in that
chain. To report such a difference once, vg compares each called allele with the reference allele
in **symbolic** form: the walk with each crossing of a child chain replaced by one symbol for that
chain. A called allele whose symbolic form equals the reference allele's is written as the
reference allele at this site, and the records of the child chain's sites report the difference.
When every called allele is written as the reference allele, the site has no record of its own,
and the records of its child chains' sites report all its differences.

#### Block records

A called allele can also differ from the reference allele in several places, separated by nodes that
both share. `--atomize-blocks`, on by default with nested calling, then writes one record per
difference. `--no-atomize-blocks` turns it off, and so does `--genotype-snarls`, which writes the
same records for every sample. vg aligns the symbolic form of each called allele to that of the
reference allele, minimising edit distance. Each maximal stretch of nodes and chain symbols that the
alignment does not match is a **block**, and becomes a record. A block record's `GT` gives the
allele that each strand carries over that block. A strand's block allele spells the strand's own
walk inside its blocks and the reference's over the steps it matches. A matched chain symbol says
only that the strand crosses the same chain, perhaps by another walk, and the records of that
chain's sites report that walk.

Blocks of the two strands that overlap or touch on the reference become one record, so that each
stretch of reference appears in at most one record. Where that makes more than one record, the
blocks replace the record for the whole site. Where it makes one record, the block replaces it only
where the site's record would say more than the block does: where a strand's site allele crosses a
matched chain by a walk that spells other bases. That chain's own records report the walk, so the
site's record would report it a second time. The block is still not written where two of its
alleles would stand for one site allele (see below), which happens where two strands' walks spell
one site allele but differ in the block, nor where the reference or a strand crosses a chain more
than once, since the chain is genotyped from each allele's first crossing only. Otherwise, where
the site cannot be split into blocks, and where `-L` merged two of the site's called alleles, vg
writes one record for the whole site.

`INFO/SB` gives each block record's index among the site's block records, counting from 0, and the
number of block records the site writes. A block record's ID (the VCF `ID` column) is the site's ID
followed by `_` and that index, so that no two records share an ID. The block records share the
site's evidence, because the likelihood is computed for the whole site. Their `GQ`, `GQI`, `GQN`,
`GP`, `QUAL`, `DP`, `DR` and `BL` are the site's values (see [VCF fields](#vcf-fields)), except
that where the linkage model moved the site, settling a genotype other than its direct call, each
block record's `GQN` and `lowconf` come from its own `GL` (see
[Linkage and re-genotyping](#linkage-and-re-genotyping)). Each block allele stands for one site
allele: the block's REF for the site's reference allele, and each ALT for the site allele of the
first strand that carries the ALT. The block record copies its `AD` and `GL` entries from the site
record, reading each block allele as the site allele it stands for. So where both strands'
different site alleles carry one ALT, the block's `GT` is homozygous and its `GL` entry is that of
the first strand's site allele, twice. Summing or averaging these fields over a site's records
therefore counts the site's evidence more than once.

When every crossing of a child chain by the alleles of the parent's settled genotype lies inside a
block, the block's ALT spells out the chain. The records of the chain, and of the sites nested in
it, would repeat that ALT, so they are not written. The chain is still genotyped and phased from
the panel, and its sites get anchors. Each linkage pass decides again, from the parent's chosen
genotype, whether a block spells the chain out.

### VCF fields

| Field | Meaning |
|---|---|
| `GT` | the settled genotype, phased where phasing ran |
| `GL` | $\log_{10} \mathcal{L}(G)$ for every genotype of the record's alleles, in the order VCF specifies for the record's alleles |
| `GQ` | the difference between the log-likelihoods of the direct call and the **runner-up**, the genotype with the second-highest $\mathcal{L}(G)$, in phred units. It is multiplied by the **explained share**, the fraction of reads whose best allele is in the direct call (a tied read split as for `AD`), and by the `--depth-quality` factor where that applies |
| `GQI` | the same difference, with neither factor |
| `GQN` | the same difference divided by the **achievable gap** (below), held at 1 or less, and multiplied by the explained share |
| `GP` | one value: the natural log of the posterior probability of the direct call, computed from $\mathcal{L}(G)$ with a uniform prior over genotypes, and still the direct call's on a record the linkage model moved. (In the VCF specification, `GP` is a phred-scaled value per genotype.) |
| `QUAL` | phred-scaled posterior probability, under the same uniform prior, of the genotype whose alleles are all the reference allele; 0 when `GT` is all reference |
| `DP` | number of reads of the site |
| `AD` | for each allele in the record, the number of reads whose best allele it is, rounded; a read tied between alleles counts a fraction to each |
| `DR` | the **depth ratio**: the site's effective read count over the expected read count of the genotype the record is written with (see [Depth term inputs](read-likelihood-direct-genotyping.md#depth-term-inputs)) |
| `BL` | mean over reads of each read's best log-likelihood score at the site, $\max_a \ell_{ra}$, in nats (see [Relative likelihood](read-likelihood-direct-genotyping.md#relative-likelihood)) |
| `FORMAT/PS` | the phase set |
| `INFO/SB` | not strand bias: on a block record, its index counting from 0 among the site's block records, and the number of block records the site writes (see [Reporting each difference once](#reporting-each-difference-once)) |
| `FILTER=noreads` | the site had no reads, so no genotype is called (`GT` is `./.`). Such a record is written only with `--genotype-snarls`, which also writes reference calls |
| `FILTER=lowconf` | `GQN` is below `--min-confidence`, when that is above 0 |

The per-site `GQ` is a difference of log-likelihoods, where the VCF specification's `GQ` is a
posterior. On records that the linkage model moved (below), `GQ` comes from a posterior instead, so
the `GQ` values of the two kinds of record are not comparable. No field flags a moved record, but a
moved whole-site record usually has a negative `GQN`. A read whose best allele is in neither the
direct call nor the runner-up fits both about equally, and adds little to their difference. The
explained share lowers `GQ` where the direct call leaves such reads unexplained. `GL`, `GQ`, `GQI`,
`GP` and `QUAL` are over-confident at high depth, because reads are treated as independent. `AD`
need not sum to `DP`: every candidate allele was scored, but only the alleles written in the record
have an entry.

`DR` is computed for the genotype the record is written with, the settled genotype. A value near 1
means that the site has as many reads as that genotype predicts. `DR` is written whether or not the
depth term is on, and is left out where the site's read-start rate is 0. Where an allele of the
written genotype was not among the site's candidate alleles, or the record's ploidy is not the one
the site was genotyped at, it is the direct call's.

#### Normalised quality

`GQ` depends on depth, since the likelihood difference is a sum over reads. It also depends on
ploidy. At ploidy 1 the runner-up is a different allele, and most reads can tell it apart from the
direct call. At ploidy 2 the runner-up usually differs from the direct call on one strand only, so
fewer reads tell the two apart, or each read tells them apart less.

`GQN` removes both effects by dividing by the **achievable gap**. This is the difference that the
[read term](read-likelihood-direct-genotyping.md#read-term) alone, the part of $\mathcal{L}(G)$ from
the reads' fits to the alleles, would give between the direct call and the runner-up if each of the
site's reads were ideal. Ideal reads are shared among the direct call's haplotypes in proportion to
its [mixture weights](read-likelihood-direct-genotyping.md#mixture-weights). Each ideal read has
$e_r = \epsilon_{\min}$ (`--mismap-min`), and fits its own haplotype's allele with relative
likelihood 1 and every other allele with 0.

The numerator of `GQN` is the observed difference in $\ln \mathcal{L}$,
[depth term](read-likelihood-direct-genotyping.md#depth-term) included. It can exceed the achievable
gap, so the fraction is held at 1. `GQN` is `.` when there is no difference to normalise (no reads,
or a single possible genotype). Otherwise it lies in $[0, 1]$, except on moved records.

#### Linkage and re-genotyping

A record is **moved** when its settled genotype is not its direct call. Its `GQ` and `GQN` are
computed again for the settled genotype, from the direct call's explained share, `GQ` factor and
achievable gap, which the linkage model keeps for each site. The direct call is the one the site
entered the linkage model with: its call in the direct pass, or, for a site of a nested chain whose
ploidy a linkage pass sets from its parent's chosen genotype, its call at that ploidy. Under
re-genotyping it is the genotype with the highest corrected likelihood in the last round.

A moved record's `GQ` is $-10 \log_{10}(1 - \text{posterior})$ times the direct call's `GQ`
factor, then capped at `GQI`. The posterior is the linkage model's posterior probability of the
settled genotype, so it includes the panel's prior. The `GQ` factor is what the per-site `GQ`
multiplies its difference by: the explained share unless `--no-share-quality`, times the
`--depth-quality` factor where that applies. `GQI` is the reads' confidence in their own best
genotype, the most they support any call at the site, and the cap keeps the prior from claiming
more where the reads are weak.

A moved record's `GQN` is the margin of the settled genotype's entry in `GL` over the largest other
entry, divided by the direct call's achievable gap, both in phred units, and multiplied by its
explained share, as the per-site `GQN` is, and held within $[-1, 1]$. It is negative where `GL`
favours another genotype over the settled one, as it does on a moved whole-site record whose direct
call is among the record's genotypes. It is `.` where the direct call had no achievable gap.
`lowconf` is decided again from this `GQN`, so a `--min-confidence` above 0 marks such a negative
record, and is cleared where `GQN` is `.`. Under re-genotyping, a record is moved if the last
round's linkage pass moved it (see [Rounds](#rounds)).

Re-genotyping is applied with `--regenotype` when `--regeno-passes` is above 1 (see
[Rounds](#rounds)). `GL` and `QUAL` are then written from the corrected likelihoods. In a round
whose correction changes the genotype with the highest $\mathcal{L}(G)$, `GQ` is recomputed from
them as the per-site `GQ` is, for that genotype, but with the direct pass's `--depth-quality`
factor. Each round starts again from the direct pass's likelihoods and `GQ`, so `GQ` follows the
last round's correction, and is the direct pass's where that correction left the best genotype
alone. `GP`, `GQI`, `GQN` and `lowconf` keep their values from before the correction, so they
describe the direct call made from the uncorrected likelihoods. `DR` depends only on the reads and
the written genotype. A moved record takes the `GQ`, `GQN` and `lowconf` described above, whether or
not re-genotyping ran.

#### Options that change the fields

`--cluster` merges similar called ALT alleles, as `vg deconstruct --cluster` does, at sites at
least `--cluster-min-len` long. Two alleles are merged when their similarity, as
`vg deconstruct --cluster` measures it, is at least the value of `--cluster`. `GT` is rewritten
for the merged alleles, their `AD` entries are summed, each merged genotype's `GL` entry is the
largest of those it replaces, and `INFO/MAT` records the merge. It is off by default.

`--no-share-quality` leaves the explained share out of the per-site `GQ`; `GQN` keeps it.
`--depth-quality` multiplies the per-site `GQ` by $e^{-D_q \vert \ln \mathrm{DR} \vert}$, where
$D_q$ is its value, at records where an allele of the direct call differs in length from the
reference allele by at least a fixed number of bases (see [Fixed constants](#fixed-constants)).
`--min-confidence` marks records with `FILTER=lowconf` and keeps them in the output. These three
options change only quality fields and `FILTER`.

#### Nesting tags

`vg call` has three other ways of genotyping nested sites, chosen by general options that also work
with `--read-likelihood`: `--all-snarls` genotypes every site on its own, `--top-down` genotypes
children after their parents and takes a child's candidate alleles from its parent's genotype, and
`--bottom-up` genotypes children before their parents. With any of them, and whenever
[off-reference chains](#which-child-chains-are-genotyped) are genotyped, records carry vg's nesting
INFO tags:

- `INFO/LV` counts the record's enclosing sites that have records on the record's own contig.
- `INFO/CH` counts the changes of contig on the way up from the record through the records of its
  enclosing sites, or gives the level of the record's contig if that is larger. A reference path
  that is not a [gRef fragment](#graph) is level 0, a gRef fragment whose ends attach to a level-0
  path is level 1, one whose ends attach to a level-1 fragment is level 2, and so on. A nested
  site that the reference path crosses stays on its parent's contig, even where a called allele
  deletes it. When the reference paths include gRef fragments, a site that only inserted sequence
  holds is reported on one, so its `INFO/CH` is at least 1.
- `INFO/PS` (parent snarl) is the ID of the nearest enclosing site that has a record. Where that
  site is written as [block records](#block-records), their IDs carry it before the `_`. It is
  unrelated to `FORMAT/PS`.
- `INFO/RC`, `INFO/RS` and `INFO/RD` give a contig, start and end at which to look the record up,
  normally those of its outermost enclosing site; the VCF header gives the details.

### Mosaic (`--mosaic-out`)

The mosaic file describes each of the sample's strands as a walk through the graph, and says which
panel haplotype the strand copies along each part of the walk. The file is the phasing in another
form, so `--mosaic-out` is rejected without the linkage model or with `--no-phased`.

The file is tab-separated. Header lines start with `#` and are identified by their first field:

| Key | Meaning |
|---|---|
| `#mosaic-version` | the format version |
| `#graph` | the input graph |
| `#sample` | the sample name, set with `--sample` |
| `#reference` | a reference path that positions refer to, one line for each, except gRef fragments |
| `#gref-fragments` | the number of gRef fragments among the reference paths, when there are any |
| `#decoding` | how the strands were chosen; its only value is `constrained-viterbi`, the Viterbi path restricted to the settled genotypes, as in [Phasing from the panel](read-likelihood-linkage-model.md#phasing-from-the-panel) |
| `#patch`, `#nested`, `#unexplained` | the choices of the three options described below |
| `#haplotype` | a panel haplotype's index and name, one line for each |
| `#note` | text describing the columns |
| `#H` | the column names |

A reader should skip header keys it does not recognise.

#### Segments and rows

A data line, or **row**, starts with `H`. Rows are built from **segments**. A segment is a maximal
stretch of consecutive sites on one strand, over which the strand copies one panel haplotype.
Consecutive segments of a strand join end to end where the walk can continue from one to the next.
Each maximal walk so formed is a **mosaic fragment**, and a strand can have several. A new mosaic
fragment starts after a gap in the walk left unfilled (see [Forming segments](#forming-segments)),
where the direction of travel reverses, as at an inversion, and on each side of a row that cannot
be walked, such as a `*` row.

| Column | Meaning |
|---|---|
| `contig` | reference contig of the row's sites |
| `strand` | `0` or `1`, the sample's strand; strand 0 carries the allele to the left of the `\|` in `GT` |
| `fragment` | the mosaic fragment's number; with `contig` and `strand`, it identifies one walk |
| `ref_start`, `ref_end` | reference positions of the row's first and last sites, or, for a `ref` row that fills a gap, of the sites on either side; approximate, since the nodes define the row |
| `start_node`, `end_node` | oriented node IDs (node ID times 2, plus 1 if reverse) where the row starts and ends |
| `hap_index` | the panel haplotype's index in the `#haplotype` lines; `ref` on a row that follows the reference instead of a haplotype the strand copies; `*` on a row where the strand is on the wildcard |
| `haplotype` | the panel haplotype's name, as its sample name and haplotype number joined by `#`; for a `ref` row, the reference's name in the panel; `*` on a wildcard row |
| `sites` | number of called sites the row covers, or `.` for a `ref` row that fills a gap; a `ref` row that replaces a segment keeps its count |
| `gbwt_offset` | with `start_node`, a position in the graph's GBWT from which the haplotype can be followed to `end_node`; `.` if there is none |

Within a mosaic fragment, each row's `end_node` is the next row's `start_node`, so the fragment
expands to one walk in the graph, counting each shared node once. Each row lies within one of the
paths that store its panel haplotype, so a segment over several such paths gives several rows.
`gbwt_offset` is valid only for the graph named in `#graph`, and `hap_index` only within one file;
`haplotype` is the name to compare across files. A contig whose phased sites are all haploid has
rows for strand 0 only. In a haploid region of a diploid contig, strand 1 is on the wildcard.

#### Forming segments

A strand's sites include the nested sites on its walk, so a change of panel haplotype at a nested
site starts a new segment. Where read phasing reversed a site's order, the two strands' panel
haplotypes are swapped there too, so each strand keeps the haplotype that carries its allele, and
can start a new segment there. Only sites with a record count, whether written as one record or as
[block records](#block-records); a site with no record is left out, so `--genotype-snarls`, which
writes records for more sites, can change the segments.

Between two consecutive segments, the walk follows the first segment's panel haplotype if that
haplotype continues to the second segment. Otherwise it follows the second segment's haplotype, if
that haplotype reaches back to the first. Where neither does, the stretch between them is a **gap in
the walk**. The gap is filled with the reference, on a `ref` row, if the reference is a panel
haplotype that crosses it. The linkage model can have a strand copy a panel haplotype at a site that
the haplotype does not pass through, as an
[unknown allele](read-likelihood-linkage-model.md#wildcard-haplotype), so a segment's haplotype
need not cover the whole segment. Such a segment is replaced by a `ref` row where the reference
crosses it, and is otherwise a row that cannot be walked.

Three options change how the rows are formed, and the header records each choice:

- `--no-mosaic-patch-gaps` turns off both uses of the reference, so gaps stay unfilled, and a new
  mosaic fragment starts after each.
- `--no-mosaic-nested` leaves nested sites out of the segments, so a strand follows its enclosing
  site's haplotype through them.
- Where a strand is on the [wildcard](read-likelihood-linkage-model.md#wildcard-haplotype), the
  panel cannot name a haplotype for it. By default those sites are left out of the segments, and
  the walk crosses them by the rule above. `--mosaic-break-unexplained` writes a `*` row for them
  instead.

### Assembly anchors (`--anchors-out`)

The anchor file records, for each genotyped site, which reads support which of the sample's
strands, for use in pangenome-guided assembly. `--anchors-out` also turns on the genotyping of
[off-reference chains](#which-child-chains-are-genotyped), which can change the VCF's phases
and, under `--regenotype`, its genotypes; `--no-off-ref-nesting` turns that off.

A **pin** is a point between two adjacent bases of the graph, with no sequence of its own. Each site
has two. The **start pin** lies just after the site's start boundary node, where an allele's walk
enters the interior. The **end pin** lies just before its end boundary node, where the walk leaves.

A site's reads are divided among its **slots**. A slot is a strand, numbered as a field of `GT`:
slot 0 the first allele, to the left of the `|`, and slot 1 the second. But a diploid site that is
not phaseable pools both strands in one slot, 0, unless `--anchors-hom-split` divides its reads
between slots 0 and 1, which then carry the same allele. A haploid site also has a single slot, 0,
except at a site of a nested haploid chain whose `GT` is `.|a`, where it is 1. An **anchor** is one
pin together with the reads of one slot that cross it.

#### Placing reads in slots

For each distinct called allele $a$, let $x_{ra} = (1 - e_r) v_a p_{ra}$, where $v_a$ are the
allele-length weights of [What each read says](read-likelihood-read-phasing.md#what-each-read-says),
and $v_a = 1$ at a site that is not phaseable. Each read goes to the slot whose allele has the
largest $x_{ra}$, except where the read phase decides (below). A read's **anchor confidence** is
$-10 \log_{10}\left(1 - x_{ra} / (\sum_b x_{rb} + e_r)\right)$, where $a$ is the allele of the slot
it is placed in. The sum runs over distinct alleles, so the two slots of a split site count their
allele once. A site's **anchor reliability** is the mean anchor confidence of the reads written for
it, after the read filters below, each read counted once: paired mates share a name, and count at
the higher of their confidences.

#### Using the read phase

With `--read-phasing`, read placement also uses each read's tempered strand log-odds $y_{rs}$ (see
[Tempering](read-likelihood-read-phasing.md#tempering)), computed leaving the site out. The
temper $\tau$ is the one re-genotyping used, where it ran and $\tau$ was above 0. Otherwise it is
fitted in the same way from the final phase.

- At a phaseable site, a read's slot is by default chosen from its $x_{ra}$ with $v_a$ replaced by
  the per-read weights $\pi$ of [Likelihood
  correction](read-likelihood-read-phasing.md#likelihood-correction). The read's anchor confidence
  is still computed from the $x_{ra}$. `--no-anchors-phase-hets` chooses the slot from the $x_{ra}$
  alone. `--anchors-strict-hets` chooses it from the sign of $y_{rs}$ alone, slot 0 for a positive
  sign, and a read whose $y_{rs}$ is 0 then keeps the slot its $x_{ra}$ give.
- `--anchors-hom-split`, which needs `--read-phasing`, divides the reads of a diploid site that is
  not phaseable between two slots by the sign of their $y_{rs}$. A site is split only when, for each
  sign, at least `--split-min-side` reads have a $y_{rs}$ of that sign and of absolute value at
  least `--split-min-q`. A split site leaves out a read whose strand is not usable there, because
  the read was used in another phase set or in more than one: its strand belongs to another phase
  set's labelling. Any other read whose $y_{rs}$ is 0 is placed by a coin flip derived from its
  name, so it takes the same slot at every site.

#### Rows of the anchor file

The file is tab-separated. Each anchor is written as an `A` row, followed by one `R` row for each
of its reads. Header lines start with `#`, and a reader should skip keys it does not recognise.
Among them, the `#read` lines list the read names, which `R` rows refer to by index. The `#H` lines
name the columns of each kind of row, and the `#note` lines describe them.

| `A` column | Meaning |
|---|---|
| `node` | the graph's ID of the pin's boundary node: the first node of `snarl` for the start pin, the second for the end pin. Under `--translation` or `--gbz-translation`, `snarl` names the translated segments, while `node` stays the graph's ID |
| `snarl` | the site's ID, as in the VCF `ID` column: its start and end boundary node IDs, each preceded by `>` or `<` for its orientation. It is also written for off-reference sites, which have no VCF record |
| `slot` | the slot |
| `allele` | the slot's allele, as its index in the site's list of candidate alleles rather than a VCF allele number. The list is not written, so the index serves to compare a site's slots |
| `gqn` | the `GQN` that a whole-site record for the site would carry, moved or not; `.` where it has none |
| `explained` | the site's explained share |
| `reliability` | the site's anchor reliability |

| `R` column | Meaning |
|---|---|
| `read_id` | the read's index in the `#read` lines; paired mates share a name, and so an index |
| `strand` | the direction in which the read crosses the pin, not one of the sample's strands: 0 in the site's direction, 1 against it |
| `offset` | the 0-based offset, in the read as sequenced, of the read's last base before the pin, reading in the site's direction |
| `score` | the read's anchor confidence |

#### Selecting sites and reads

- `--anchors-het-only` writes anchors only at phaseable sites, leaving out the other diploid sites
  and haploid sites, including the sites of nested haploid chains. `--anchors-leaf-only` writes
  them only at sites with no child chains.
- `--anchors-min-gqn`, when above 0, skips sites whose `GQN` is below its value or missing.
- `--anchors-min-q` skips reads whose anchor confidence is below its value.
- `--anchors-keep-off-call` keeps reads whose best candidate allele was not called; by default they
  are left out.
- `--anchors-reads` skips anchors with fewer reads than its value.
- With `--anchors-end-new` set to $N$, a slot's end anchor is written only if at least $N$ of the
  slot's reads cross the end pin but not the start pin, which drops end anchors that mostly repeat
  the start anchor's reads.

## Options

`--preset ont` sets several of the options below to values suited to Oxford Nanopore reads, and
`vg call --help` lists which, and each option's default. An option given explicitly overrides the
preset. Where an option has a `--no-` form, the two set the same thing and the one given last wins,
so a `--no-` form also turns off a setting that a preset turned on.

The table includes general `vg call` options: `--pack`, `--gbwt`, `--gbz`, `--max-snarl-edges`,
`--nested`, `--no-nested`, `--atomize-blocks`, `--no-atomize-blocks`, the Ploidy row, `--cluster`
and `--cluster-min-len`. Its other options are rejected without `--read-likelihood`. The options
that modify `--anchors-out` (`--no-off-ref-nesting` among them), `--mosaic-out` and `--regenotype`
are also rejected when those are not in use. Other general options used on this page, such as
`--genotype-snarls`, `--sample`, `--snarls`, `--ref-path`, `--path-prefix`, `--ref-sample`,
`--all-snarls`, `--top-down`, `--bottom-up`, `--translation` and `--gbz-translation`, are described
by `vg call --help`.

| Part | Options |
|---|---|
| [Reads](read-likelihood-direct-genotyping.md#read-input) | `--gam`, `--gaf-reads`, `--gam-index`, `--gaf-index`, `--gaf-base`, `--gbz-base`, `--gaf-base-binary`, `--read-window`, `--read-min-mapq` |
| [Candidate alleles](#candidate-alleles) | `--enumerate-support`, `--pack`, `--gbwt`, `--gbz`, `--max-snarl-edges` |
| [Relative likelihood](read-likelihood-direct-genotyping.md#relative-likelihood) | `--gap-open`, `--gap-extend`, `--insertion-nats`, `--optimal-pairing`, `--no-optimal-pairing` |
| [Mismapping](read-likelihood-direct-genotyping.md#mismapping-probability) | `--mismap-min`, `--mismap-max`, `--no-mismap-term` |
| [Mixture weights](read-likelihood-direct-genotyping.md#mixture-weights) | `--flat-mixture` |
| [Depth term](read-likelihood-direct-genotyping.md#depth-term) | `--depth-term`, `--depth-count-raw` |
| [Linkage](read-likelihood-linkage-model.md) | `--linkage-weight`, `--linkage-scale`, `--linkage-prior`, `--hp-prior`, `--hp-prior-run` |
| [Nested sites](#nested-sites) | `--nested`, `--no-nested`, `--atomize-blocks`, `--no-atomize-blocks`, `--no-off-ref-nesting` |
| [Ploidy](#ploidy) | `--ploidy`, `--ploidy-regex`, `--ploidy-bed` |
| [Phasing](#phasing) | `--phased`, `--no-phased`, `--read-phasing`, `--no-read-phasing`, `--phase-min-q`, `--phase-coherence`, `--phase-coh-rounds`, `--phase-break`, `--phase-relink`, `--phase-hang`, `--phase-prior`, `--phase-cap` |
| [Re-genotyping](#re-genotyping-from-the-phase) | `--regenotype`, `--no-regenotype`, `--regeno-temper`, `--regeno-ceiling`, `--regeno-passes`, `--regeno-haploid`, `--no-regeno-haploid`, `--regeno-ledger` |
| [Quality fields](#options-that-change-the-fields) | `--no-share-quality`, `--depth-quality`, `--min-confidence` |
| [Allele merging](#options-that-change-the-fields) | `--cluster`, `--cluster-min-len` |
| [Mosaic](#mosaic---mosaic-out) | `--mosaic-out`, `--mosaic-patch-gaps`, `--no-mosaic-patch-gaps`, `--no-mosaic-nested`, `--mosaic-break-unexplained` |
| [Anchors](#assembly-anchors---anchors-out) | `--anchors-out`, `--anchors-reads`, `--anchors-min-gqn`, `--anchors-min-q`, `--anchors-het-only`, `--anchors-leaf-only`, `--anchors-keep-off-call`, `--anchors-end-new`, `--anchors-phase-hets`, `--no-anchors-phase-hets`, `--anchors-strict-hets`, `--anchors-hom-split`, `--split-min-q`, `--split-min-side` |
| Debugging and evaluation | `--dump-likelihoods` (writes each read's $e_r$ and $p_{ra}$ at every site to a file, as TSV), `--regeno-shuffle` (randomises the sign of each $\Lambda_{rs}$ before re-genotyping, as a control for how much the phase contributes); `--flat-mixture` and `--anchors-strict-hets`, listed above, also serve to measure the parts they replace |
| Presets | `--preset` |

### Fixed constants

These constants have no option. Each is located by its name, or by the function or member that
holds it.

| Constant | Defined in | What it sets |
|---|---|---|
| `LinkageModel::Params::escape` | `src/linkage_model.hpp` | the escape penalty $\epsilon_{\mathrm{esc}}$ on a strand with an [unknown allele](read-likelihood-linkage-model.md#wildcard-haplotype) |
| `LinkageModel::Params::rho_min` | `src/linkage_model.hpp` | the floor $\rho_{\min}$ in the switch probability (see [Transitions](read-likelihood-linkage-model.md#transitions)) |
| `LinkageModel::Params::window`, `margin` | `src/linkage_model.hpp` | the number of sites each forward–backward window keeps, and the extra sites decoded and discarded on each side |
| `RATE_BUCKET`, `RATE_ID_WINDOW` | `src/allele_likelihood.hpp` | reference bp per rate-window bucket (a rate window is three buckets), and node IDs per rate window when there are no reference positions |
| `banded` and `band`, in `score_by_optimal_pairing` | `src/allele_likelihood.cpp` | the grid size, in read visits times allele visits, above which optimal pairing is restricted to a band, and the band's half-width in allele visits |
| `max_yens_traversals` | `src/subcommand/call_main.cpp` | the most candidate alleles support enumeration keeps |
| the difference limit in `LinkageModel::run_length_site` | `src/linkage_model.cpp` | the largest homopolymer length difference to which `--hp-prior` applies |
| the `int8_t` of `LinkageCollector::allele_arena` | `src/linkage_model.hpp` | the most alleles a site's compact allele set may hold for the linkage model |
| the 64-bit crossing mask of `child_crossing_mask` | `src/graph_caller.cpp` | the most candidate alleles a parent site can have for vg to tell which of them cross a child site |
| the `min_length` argument of `set_depth_quality` | `src/read_likelihood_caller.hpp` | the change in allele length at which `--depth-quality` applies |
| the minimum counted reads in `read_phase_flips` | `src/read_phasing.cpp` | the reads a site needs before low coherence can remove it from the phase chain |
| `RegenotypeParams::fit_bins`, `fit_min_per_bin`, and the grid in `fit_calibration` | `src/regenotype.hpp`, `src/regenotype.cpp` | the bins, the fewest observations per bin, and the candidate values used to fit $\tau$ |

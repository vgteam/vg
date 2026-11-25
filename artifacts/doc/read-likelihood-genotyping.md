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

The caller is built in two parts. **Direct genotyping** computes each site's likelihoods and its
direct call from the site's reads alone. **Linkage-based genotyping**, the linkage model, works on
top of it: it chooses the genotypes of many neighbouring sites together, from their site
likelihoods, and gives them their first phase. How the four steps are ordered into passes over the
sites is described under [Passes and rounds](#passes-and-rounds), after the vocabulary.

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
  linkage model, phasing and nested calling are built for at most two haplotypes. A ploidy BED, a
  table of ploidies by region given with `--ploidy-bed` or `--ploidy-preset`, may also give 0 for a
  region the sample has no copy of. A nested site that none of its parent's alleles passes through
  has ploidy 0 (see [Which child chains are genotyped](#which-child-chains-are-genotyped)).
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

vg goes over the sites in **passes**. The **direct pass** comes first and runs once. It fetches
the reads: the steps after it work from the evidence it keeps for each read, except for the
[read-consensus step](#alleles-built-from-the-reads---read-consensus), which fetches the reads of
the sites it tests.
**Rounds** follow, each built around a **linkage pass**. Last, the **render** builds the records,
and vg writes the output files:

```
direct pass                         once, every top-level site
    direct genotyping               the site's likelihoods and direct call, from its reads
    descent                         the same for the sites nested in it, at provisional ploidies
    staging                         store the site's alleles, direct call and likelihoods
round 1, 2, ...                     round 2 on only with --regenotype
    likelihood correction           from round 2 on: correct the likelihoods from the phase
    linkage pass                    one level at a time, level 0 first
        linkage-based genotyping    choose a genotype for each site of the level
        phasing from the panel      put each genotype's alleles on the strands
        ploidy of the next level    from the genotypes just chosen
    read phasing                    re-decide the phases from the reads (--read-phasing)
    read consensus                  round 1 only, with --read-consensus: add alleles built from
                                    the reads, then a linkage pass and read phasing again
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
linkage pass drops the site. Each site is **staged**: vg stores its candidate alleles, its direct
call and the likelihoods of its genotypes (at both ploidies, for a nested site that has both), and
builds the site's records from them in the render.

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
computed and reported on standard error, but not applied, and no linkage pass follows.

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
Where the site's candidate alleles have groups (see
[Choosing among alleles that differ inside child chains](#choosing-among-alleles-that-differ-inside-child-chains)),
the [grouped likelihood](read-likelihood-direct-genotyping.md#groups-of-alleles)
$\tilde{\mathcal{L}}(G)$ takes the place of $\mathcal{L}(G)$ here, in the [VCF fields](#vcf-fields),
in the linkage model and in [re-genotyping](#re-genotyping-from-the-phase), whose corrections are
added to it. `GQN` then divides a difference in $\ln \tilde{\mathcal{L}}(G)$ by the achievable gap
of $\mathcal{L}(G)$, which takes no account of the groups (see
[Normalised quality](#normalised-quality)).

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

With [`--read-consensus`](#alleles-built-from-the-reads---read-consensus), a snarl can also gain
an allele built from its reads after the direct pass.

`--max-snarl-edges` sets a limit on a site's edges, counting the edges of the sites nested in it.
vg skips a site over the limit, and genotypes the sites of its child chains as top-level sites.
They take the ploidy of their contig or region, and join the contig's top-level sites in the
linkage model. A skipped site with no child chains gets no genotype and no record, even with
`--genotype-snarls`. The limit exists because Yen's search is slow on very large sites. Under
haplotype enumeration, which does not use Yen's search, there is no limit by default.

### Alleles built from the reads (`--read-consensus`)

Haplotype enumeration offers only the alleles of the panel, and
[off-panel detection](read-likelihood-direct-genotyping.md) adds changes inside nodes. With
`--read-consensus`, off by default, a snarl whose reads carry a longer sequence that no walk spells
gains that sequence as an allele. It needs `--phased`. The step runs once, at the end of round 1,
on each snarl that has a phase, at any level, and on each of its strands:

1. The strand's reads are the reads that **span** the snarl, from its start bound's node to its
   end bound's, and that [read phasing](#from-the-reads) puts on the strand: those whose strand
   log-odds, without the snarl's own term, are positive for strand 0 and negative for strand 1.
   Where fewer than half of the spanning reads have a strand, as in a long homozygous stretch or
   without `--read-phasing`, all of them are taken to be on strand 0. A snarl with one strand has
   all of them on it.
2. The strand's **called sequence** is its allele with each site nested in the snarl spelled by
   the allele the strand has there, as a [flat record](#flat-records) spells it. A read's
   **gain** by another sequence is what the read's bases between the bounds score more against
   that sequence than against the called sequence, in nats, each scored as a read is against an
   allele given as its bases. A gain is no more than **one read's say**, the most by which the
   genotyper lets one read favour an allele (the negative log of
   [`--mismap-min`](read-likelihood-direct-genotyping.md#mismapping-probability)), and no less than that below zero. So one or two
   reads cannot carry an allele.
3. One threshold, `--consensus-min-gain` (9 nats), decides what is built and what is kept. A
   strand gets a consensus where it has enough spanning reads to gain that much, each by one
   read's say: 4 reads with the `ont` preset.
4. The consensus is that of the strand's spanning reads. At most 20 are used, those
   nearest the median length. Each is aligned to the called sequence, each
   base becomes whichever of a base or a deletion most reads have there, and an insertion is made
   where more than half the reads insert bases; the result is aligned to again until it stops
   changing. A difference that only changes the length of a run of one base is undone where the
   run has at least nine bases in the sequence now called; a shorter run takes the reads' length.
   The called sequence is the consensus where the first round changes nothing else and the
   reads are neither all longer than the called sequence nor all shorter. The strand is a
   **candidate** where the consensus differs from the called sequence.
5. A consensus by which the strand's reads gain `--consensus-min-gain` nats in all
   becomes an allele of the snarl, unless one of the snarl's alleles spells it already. Two strands
   with the same consensus add one allele. With a consensus, the snarl also gains what each
   strand is called as, where the sites nested in the snarl make that a sequence none of its
   alleles spells, so that a strand can keep what it has. No panel haplotype carries a
   consensus allele.
6. The consensus allele is the walk the strand is called as, with changes to the bases of its
   steps where the consensus differs. The two are aligned and the gaps moved left, and each
   difference changes the steps that hold its bases: a node whose bases are all deleted stays
   in the walk with none, and an insertion goes after the base before it, or before the first
   base after the start bound where it comes first. Where the snarl has sites nested in it, or
   no walk spells what the strand is called as, the allele is instead given as its bases between
   the bounds, as is what a strand is called as. Such an allele visits no node between the
   bounds, so it crosses no site nested in the snarl, as a deletion of everything between the
   bounds does.

A snarl gains its consensus whatever the snarls nested in it gain. The two are not called
together on one strand: the consensus allele of such a snarl crosses no nested site, so a strand that takes it has
no allele at the snarls nested in the snarl, and the linkage pass chooses.

A snarl that gained an allele is genotyped again from its reads, at the ploidy it has. Each read
that spans it is scored by aligning its bases between the bounds to every allele's (see
[Relative likelihood](read-likelihood-direct-genotyping.md#relative-likelihood)), and the other
reads are scored as before and say nothing about the new allele. Where [groups](#choosing-among-alleles-that-differ-inside-child-chains)
choose the genotype, the new allele is a group of its own. The snarl is then filed with the linkage
model in place of its first call, and a linkage pass and read phasing run again before any further
round. So the linkage model, to which the allele is one no haplotype carries, read phasing and
re-genotyping treat the snarl as if the direct pass had found the allele.

#### Runs of neighbouring snarls

In a tandem repeat the mapper can place one difference in any of several neighbouring snarls of a
chain, and it does not choose the same one for every read. Then no one snarl holds
the difference in enough reads, and the bases that two reads have
between the bounds of one snarl cannot be compared. A difference can also lie in the nodes between
two snarls, where no snarl holds it. So the step tests a **run** of snarls as it tests one snarl,
and before it treats any snarl on its own. Runs follow the chain, not reference positions. A run
inside a snarl too large to be treated itself is a run of the snarls of a chain nested in it.

1. Two snarls are **neighbours** where one follows the other in a chain, at any level: they
   share a bound node, or only nodes that every walk of the chain takes lie between them. One
   phase set orders the strands of both, enough reads span each to gain `--consensus-min-gain`
   nats, and each has a called
   sequence of no more bases than the scorer aligns of a read. Their **junction** is those
   nodes' bases, the bound nodes included, and has at most 1,000 bases.
2. No allele **explains** a spanning read that, between the bounds, scores
   `--consensus-min-gain` nats or more below a sequence that spells it against the allele it
   fits best. A snarl's **signal** is one read's say for each such read. A read's **gap** is a
   run of bases it inserts in a node or beside it, or a run of a node's bases it lacks. Where
   enough of the reads that span a snarl to gain `--consensus-min-gain` nats have a gap inside
   the node of its bound that faces a neighbour, or in the nodes of their junction beyond that
   node, and a gap of that length costs a read `--consensus-min-gain` nats, the signal of both
   snarls is the say of that many reads at least. So a difference that lies in the nodes
   between two snarls, in no snarl, makes a run of the two. A **seed** is a snarl with a signal
   of `--consensus-min-gain` or more.
3. A run begins and ends at a seed, and holds every snarl between two of its seeds. Two seeds
   that follow each other in a sequence of neighbours are in one run where at most 1,000 bases
   lie between them, in the nodes of the junctions and in the snarls that are no seeds, and the
   run has no more bases between its outer bounds than the scorer aligns of a read. The bases of
   a snarl are those of its longest called sequence. Where one more seed would pass either
   number, the run ends and the next begins at that seed.
4. A run of two or more is **tested**: the strand's called sequence across the run is its
   called sequence at each member with the junctions between them, a read spans the run where it
   spans from the first member's start bound to the last member's end bound, and steps 1 to 5
   above apply to the run as to one snarl, with each read put on the strand that its strand
   log-odds at the members, summed, favour.
5. The run is **treated** where a strand has a consensus by which its reads gain
   `--consensus-min-gain` nats and the run can become [one site](#a-run-as-one-site). A snarl
   in no treated run is treated on its own where it has a consensus, as above.

#### A run as one site

A treated run becomes one site, the **merged site**, in place of its members, so that its alleles
can differ from what is called anywhere between the run's outer bounds, the nodes between two
snarls included:

- Its bounds are the start bound of the run's first member and the end bound of its last. It is
  named by them as a snarl is, has the first member's place in the nesting tree and on the
  reference and its ploidy, and has no child sites.
- Its alleles are sequences across the run. The first is the walk that takes each member's
  reference allele, with the junctions between them, which is its reference allele. Then comes
  what each strand is called as across the run, as the walk that takes the strand's allele at
  each member where that walk spells it, and otherwise as its bases between the bounds. Last come
  the consensuses, each as that walk of its strand with changes to its steps, as at one snarl,
  and as its bases between the bounds where the walk does not spell what the strand is called
  as. The nodes between two members are steps of the walk like any other. Alleles that spell the
  same bases are one,
  and a run whose consensuses are all spelled by the others stays as it is.
- It is genotyped from its reads as a snarl that gained an allele is, with every read that spans
  it aligned to every allele, and filed with the linkage model in place of its first member. The
  other members, every site nested in a member, and the
  [edit sites](read-likelihood-direct-genotyping.md) inside the nodes of the junctions are no
  longer sites: the merged site's alleles spell what was called at them. Phasing, re-genotyping
  and the records treat the merged site as a snarl. A site that holds the run's chain spells the
  merged site's chosen alleles in its own, as it does those of any site nested in it.

A run cannot become one site, and is left alone, where its members are not placed alike on the
reference, or where a change is called at an edit site inside a node of a junction, which the
strands' called sequences across the run do not have.

#### Which snarls are tested, and limits

The step fetches reads only at the snarls whose reads can gain by an allele. The direct pass
keeps, for each snarl, the number of reads that span it, the number of those that no allele
explains, and the gaps its reads have around each bound. A snarl is tested on its own where its
signal is `--consensus-min-gain` or more, and in a run where the run is tested.
The scorer aligns a read with at most 10,000 bases between the bounds, and makes at most 2,000
alignments at one snarl (see [Fixed constants](#fixed-constants)). A snarl with a called sequence
or a consensus of more bases than that, or whose reads and alleles need more alignments, gains no
allele, and such a run is left alone. A read with more bases is scored as one that does not span
the snarl. An allele can have any number of bases. The step fetches
reads without adding to the read-start counts of
the [depth term](read-likelihood-direct-genotyping.md#depth-term) or to the off-panel candidates,
so both are as the direct pass left them.

A top-level site with such an allele called at it, or at a site nested in it, is written as
[flat records](#flat-records) whether or not `--flatten-records` is given: the allele is a
sequence that no walk spells, so one record of the site would hold every base between its bounds
for a difference of a few. The site is a run of its own, and its records have `INFO/RC` and
`INFO/FL`. Every other site is written as without `--read-consensus`. A site that cannot be
written as flat records is written as one site record with `INFO/RC`, not as
[block records](#block-records): the allele puts its bases after the start bound's node, which
no block covers. Its `AT` walk is the two bound nodes alone. vg reports on standard error how many
snarls it tested on their own and how many strands were candidates at them, how many runs were
treated and how many of those became one site, how many sites gained an allele or were made of a
run and how many of them are called with such an allele, how many consensus alleles of those
sites are a walk with changes and how many are given as bases, because no walk spells what the
strand is called as or because the snarl has sites nested in it, and how many snarls were over a
limit of the scorer.

## Ploidy

A site's ploidy comes from a ploidy BED of `CHROM START END PLOIDY` lines, given as a file with
`--ploidy-bed` or by name with `--ploidy-preset`. `CHROM` is spelled as in the output VCF.
Intervals are 0-based and half-open, and must not overlap. A site takes the ploidy of the interval
that contains its position. A site outside every interval takes the ploidy from `--ploidy`, 2 by
default, or from `--ploidy-regex`, a comma-separated list of `REGEX:PLOIDY` rules of which the
first whose regular expression matches the reference path's name applies. Every ploidy must be 1
or 2, except that an interval of the ploidy BED may have ploidy 0: the sample has no copy of that
region, so a site there is not genotyped and gets no record. A site at ploidy 1 has a single allele
as its `GT`.

`--ploidy-preset` names a ploidy BED built into vg, one per reference assembly and sample karyotype:
`chm13v2`, `grch38`, `grch37.ucsc` and `grch37.b37`, each as `XY` and `XX`, for example
`chm13v2.XY`. An `XY` preset makes chrX haploid outside its pseudoautosomal regions. It also makes
chrY and the mitochondrion haploid. An `XX` preset gives chrY ploidy 0 and makes the mitochondrion
haploid. The `grch37.b37` presets name the contigs `X`, `Y` and `MT`, and the others `chrX`, `chrY`
and `chrM`. A preset has no lines for the autosomes, so one preset serves a whole-genome run and a
run on any one contig alike, with no `--ploidy`. `--ploidy-bed` and `--ploidy-preset` cannot be
given together. The presets are in `artifacts/ploidy_regions.cpp`; for another assembly or
karyotype, copy the nearest one into a BED file, edit it, and pass it with `--ploidy-bed`.

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

Nested calling is part of `--read-likelihood`, and on by default with it. `--no-nested` turns it
off. With `--no-nested`, each site is genotyped against its full walks, variation inside nested
sites is reported in the enclosing site's alleles, and a nested site is genotyped on its own only
when its parent could not be genotyped. How nested sites are written as records is described under
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

### Choosing among alleles that differ inside child chains

Several candidate alleles of a site can take one route and differ only inside a child chain. The
reads that fit the route are then split among those alleles, so a genotype of alleles that take
other routes can fit the reads better. With `--group-candidates`, on by default, vg groups the
candidate alleles of each site that differ only inside child chains: two alleles share a group when
their [symbolic forms](#reporting-each-difference-once) are equal. The groups then choose the site's
genotype, and the alleles choose within the chosen groups, by the
[grouped likelihood](read-likelihood-direct-genotyping.md#groups-of-alleles). Grouping does not
depend on nested calling: under `--no-nested`, the site's records report the differences inside its
child chains from the alleles chosen within the groups. `--no-group-candidates` turns grouping off.

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

At some phaseable sites the panel does not order the two alleles: the two strands' haplotypes carry
the same allele there, or neither carries either chosen allele, so either strand could carry either
(see [Phasing from the panel](read-likelihood-linkage-model.md#phasing-from-the-panel)). vg still
writes them with `|`, in an order that carries no phase, and keeps the site in its phase set, so
that read phasing, when on, can order it. No field marks these sites.

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
is off by default and on under `--preset ont` and `--preset hifi`. It uses its own statistics, not
the linkage model's, and is described in [Read
phasing](read-likelihood-read-phasing.md#read-phasing).

### Re-genotyping from the phase

`--regenotype` uses the phase to give each read its own probabilities of having come from each
strand, and so corrects each site's likelihoods, at the start of every round after the first (see
[Rounds](#rounds)). Off-reference chains, and chains that a block record spells out, keep the
likelihoods of the direct pass. It needs `--read-phasing` and, like it, is off by default and on
under `--preset ont` and `--preset hifi`. The correction is described
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
genotype, whether a block spells the chain out. An ALT that both strands carry is spelled once, as
the first strand has it, so it does not spell a site at which the two strands carry different
alleles. Such a site keeps its record where the record that would spell it out is of a site whose
two strands carry one walk.

#### Flat records

`--flatten-records`, off by default, writes the top-level sites as the base-level differences of
their called strands from the reference, and writes no records for the sites nested in them. It
needs `--phased`, and cannot be combined with `--genotype-snarls`. Only the records change: the
genotypes, the phase, the anchors and the mosaic are what they are without it. Without it, the
only sites written this way are those called with an
[allele built from the reads](#alleles-built-from-the-reads---read-consensus), each as a run of
its own.

The records are cut over a **run** of sites. A run is a longest series of top-level snarls that
follow each other in one chain, where:

- the nodes from one snarl's end boundary node to the next snarl's start boundary node, both
  included, have at most `--flatten-distance` bases (1000 by default). Two snarls that share a
  boundary node have the bases of that node between them;
- one phase set holds all of them, they have the same ploidy, and each heterozygous one has a
  phased GT;
- the reference path reads them one after another; and
- the run's reference has at most 50,000 bases. A snarl that would take a run past that starts
  the next run.

A top-level [edit site](read-likelihood-direct-genotyping.md) in the nodes between two snarls of
a run is in the run, and gets no records of its own. Where such an edit site cannot be, because
it is in another phase set, has an unphased heterozygous GT, or has a node outside those nodes,
the two snarls are in different runs. Every other top-level site is a run of its own, and with
`--flatten-distance 0` every site is.

A strand's **sequence over a run** is spelled from the first site's start boundary node through
the last site's end boundary node: the allele the strand carries at each site of the run, with
each nested site spelled by the allele that strand carries there, at every depth, the changes of
the edit sites between the sites, and the nodes between the sites. The run's **reference** is the
same with every site at its reference allele, which is the reference path through the run. vg
aligns each strand's sequence to the reference by itself, end to end, at the least cost: 40 for a
mismatch, and for a gap of $L$ bases the lesser of $78 + 7L$ and $162 + 3L$. These are the scores
that minimap2 aligns two assemblies of one species with (its `asm5` preset), as costs of an
end-to-end alignment, so that the records are cut as an assembly-to-reference alignment cuts
them. It then moves each gap as far left as it goes within the run.

Each maximal run of alignment columns that are not matches is a **difference**, and the records are
cut from the differences:

- A run of mismatches is one record (an MNP where it is longer than one base).
- A gap is one record, with the reference base before it at the start of REF and of each ALT. A gap
  with mismatches touching it is one record that covers both.
- A difference that both strands have, with the same reference bases and the same bases in their
  place, is one record with GT `1|1`. Any other difference is a record of its own strand, with GT
  `1|0` or `0|1`. The `0` says only that the other strand lacks this difference: what the other
  strand has at those reference bases is in its own records. So two records of a run can cover
  the same reference bases, one for each strand, and no flat record has two ALTs.

`POS` is the position of the record's first reference base on the reference path, as in every other
record. A record belongs to the site of the run it starts in, or to the site before it where it
starts in the nodes between two sites. The `GT` lists the strands in the order of the phase. It is
phased, with that site's `PS`, where the site's own GT is phased, and is otherwise written with
`/` and no `PS`. A run whose strands both spell the reference writes nothing.

The other fields follow the [block records](#block-records). The ID of a record of a run of one
site is the site's ID, and that of a run of several is the ID of a site from the start boundary
node of the run's first site to the end boundary node of its last, in both cases followed by `_`
and the record's index among the run's records. `INFO/SB` gives that index and the number of flat
records of the run. `QUAL`, `FILTER`, `DP`, `DR`, `BL`, `GQ`, `GQI`, `GQN` and `GP` are those of
the site the record belongs to, and `AD` and `GL` are copied from that site's record through the
site allele each strand carries. A difference that comes from a nested site, from an edit site
between two sites, or from another site of the run therefore carries the evidence of the site it
belongs to. `INFO/FL` marks a flat record. `INFO/AT` is `.` for each allele, since no walk of the
graph spells a flat allele.

A site is written as it is without `--flatten-records`, and so are the sites nested in it, where it
has no strands to spell (it has no genotype or no phase, or the strands of its phase do not carry
the genotype that is written), where a gap at the site's first base has no reference base before
it, or where the
aligner reaches its limit. The aligner fills a band of the alignment matrix and widens it until it
can show that the alignment costs the least, keeping a byte for each cell of the band. A band of
more than $2^{27}$ cells is filled for one strand at a time, and the aligner gives up when the
band would hold more than $2^{30}$ cells, or the proof more than $2^{30}$. A run of several sites
whose strands cannot be aligned, or spelled along it, is written as a run for each of its sites
and for each edit site between them. The call reports on stderr, in lines starting `[vg call] flatten:`, the runs and
lines written, the sites in the runs, the runs of several sites, the nested sites and the edit
sites between sites left without a line, the runs of several sites written site by site, by
reason, and the sites not flattened, by reason.

With the environment variable `VG_CHECK_FLAT_RECORDS` set, vg checks every run once all its sites
have written their records: it puts the run's flat records back into the run's reference, taking
for each strand the allele its GT entry names, and compares the result with the strand's sequence
over the run. A `[vg call] flatten check:` line reports the runs checked and the runs that failed.

### VCF fields

| Field | Meaning |
|---|---|
| `GT` | the settled genotype, phased where phasing ran |
| `GL` | $\log_{10} \mathcal{L}(G)$ for every genotype of the record's alleles, in the order VCF specifies for the record's alleles |
| `GQ` | the difference between the log-likelihoods of the direct call and the **runner-up**, the genotype with the second-highest $\mathcal{L}(G)$, in phred units. It is multiplied by the **explained share**, the fraction of reads whose best allele is in the direct call (a read tied between alleles counts a fraction to each), and by the `--depth-quality` factor where that applies |
| `GQI` | the same difference, with neither factor |
| `GQN` | the same difference divided by the **achievable gap** (below), held at 1 or less, and multiplied by the explained share |
| `GP` | one value: the natural log of the posterior probability of the direct call, computed from $\mathcal{L}(G)$ with a uniform prior over genotypes, and still the direct call's on a record the linkage model moved. (In the VCF specification, `GP` is a phred-scaled value per genotype.) |
| `QUAL` | phred-scaled posterior probability, under the same uniform prior, of the genotype whose alleles are all the reference allele; 0 when `GT` is all reference |
| `DP` | number of reads of the site |
| `AD` | for each allele in the record, the number of reads whose best allele it is, rounded; a read tied between alleles counts a fraction to each |
| `DR` | the **depth ratio**: the site's effective read count over the expected read count of the genotype the record is written with (see [Depth term inputs](read-likelihood-direct-genotyping.md#depth-term-inputs)) |
| `BL` | mean over reads of each read's best log-likelihood score at the site, $\max_a \ell_{ra}$, in nats (see [Relative likelihood](read-likelihood-direct-genotyping.md#relative-likelihood)) |
| `FORMAT/PS` | the phase set |
| `INFO/SB` | not strand bias: on a block record, its index counting from 0 among the site's block records, and the number of block records the site writes (see [Reporting each difference once](#reporting-each-difference-once)); likewise on a [flat record](#flat-records) |
| `INFO/FL` | the record is a [flat record](#flat-records), written with `--flatten-records` |
| `INFO/RC` | the called genotype of the record's snarl, or of its [merged site](#a-run-as-one-site), has an [allele built from the reads](#alleles-built-from-the-reads---read-consensus), written with `--read-consensus`; on a flat record, the called genotype of its top-level snarl or of a snarl nested in it has one |
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
least `--cluster-min-len` long. Two alleles are merged when their similarity is at least the value
of `--cluster`. The similarity is the total length of the node visits inside the site that both
alleles make, over the total length of those either allele makes. Against a pure deletion, which
makes no visit inside the site, it is 1 minus the other allele's length inside the site over the
larger of that length and the reference allele's. `GT` is rewritten for the merged alleles, their
`AD` entries are summed, each merged genotype's `GL` entry is the largest of those it replaces,
and `INFO/MAT` records the merge. It is off by default.

`--no-share-quality` leaves the explained share out of the per-site `GQ`; `GQN` keeps it.
`--depth-quality` multiplies the per-site `GQ` by $e^{-D_q \vert \ln \mathrm{DR} \vert}$, where
$D_q$ is its value, at records where an allele of the direct call differs in length from the
reference allele by at least a fixed number of bases (see [Fixed constants](#fixed-constants)).
`--min-confidence` marks records with `FILTER=lowconf` and keeps them in the output. These three
options change only quality fields and `FILTER`.

#### Nesting tags

`vg call` has one other way of genotyping nested sites, chosen by a general option that also works
with `--read-likelihood`: `--all-snarls` genotypes every site on its own. With it, and whenever
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

### Mapping without the reference

The anchor file ties reads to the sample's strands, and the mosaic names the panel haplotypes those
strands copy, so both depend on where the mapper placed each read. `vg haplotypes` makes a graph for
the sample whose haplotypes are pieced together from panel haplotypes chosen to resemble the
sample's. With `--include-reference`, which its presets turn on, the graph also keeps the reference
paths, which were not chosen that way. The mapper can place a read that the sampled haplotypes fit
imperfectly on a reference path instead, including on sequence that only the reference visits, and
so away from the sample's alleles.

To prevent this, map the reads to a copy of the graph without the reference paths, and call on the
graph that has them, since `vg call` reports positions on the reference paths. Removing the
reference paths leaves out the nodes that no other path visits, and keeps the rest with their IDs,
so a placement in the copy is also a placement in the graph. For a `graph.gbz` whose reference
sample is `CHM13`:

```
vg gbwt -Z graph.gbz -o graph.gbwt
vg gbwt -x graph.gbz -R CHM13 -g graph.noref.gbz graph.gbwt
```

The first command copies out the GBZ file's GBWT, the index that holds its paths, because `vg gbwt`
cannot remove paths from a GBWT inside a GBZ file. The second reads the paths from `graph.gbwt` and
the node sequences from `graph.gbz`, removes the paths of sample `CHM13`, and writes the rest as
`graph.noref.gbz`. `vg paths -M -x graph.gbz` lists each path's sample and sense; give `-R` once for
each sample whose paths are marked `REFERENCE`, and when there are several, name the one to report
positions on with `--ref-sample` when calling. Then index the copy, map the reads to it, and call
on the graph. For nanopore reads:

```
vg autoindex -w lr-giraffe -G graph.noref.gbz -p map
vg giraffe -b r10 -Z graph.noref.gbz -d map.dist -m map.longread.withzip.min \
    -z map.longread.zipcodes -f reads.fq.gz -o gaf > reads.gaf
vg call graph.gbz --read-likelihood --preset ont --gaf-reads reads.gaf \
    --mosaic-out mosaic.tsv --anchors-out anchors.tsv > calls.vcf
```

`--gaf-reads` alone loads every read into memory; for a whole genome, sort and index the GAF and
give the index with `--gaf-index` (see [Read input](read-likelihood-direct-genotyping.md#read-input)).

`vg giraffe` can also sample haplotypes itself, with `--haplotype-name` and `--kff-name`. With its
default sampling preset, and with the presets that `-b hifi` and `-b r10` select, the graph it
samples and maps to keeps the reference paths.

### Mosaic (`--mosaic-out`)

The mosaic file describes each of the sample's strands as a walk through the graph, and says which
panel haplotype the strand copies along each part of the walk. The file is the phasing in another
form, so `--mosaic-out` is rejected without the linkage model or with `--no-phased`. For a graph
sampled with its reference paths kept, map the reads as described under [Mapping without the
reference](#mapping-without-the-reference).

The file is tab-separated. Header lines start with `#` and are identified by their first field:

| Key | Meaning |
|---|---|
| `#mosaic-version` | the format version: 5, or 6 for a file with `A` rows (see [Called alleles](#called-alleles---mosaic-exact)) |
| `#graph` | the input graph |
| `#sample` | the sample name, set with `--sample` |
| `#reference` | a reference path that positions refer to, one line for each, except gRef fragments |
| `#gref-fragments` | the number of gRef fragments among the reference paths, when there are any |
| `#decoding` | how the strands were chosen; its only value is `constrained-viterbi`, the Viterbi path restricted to the settled genotypes, as in [Phasing from the panel](read-likelihood-linkage-model.md#phasing-from-the-panel) |
| `#patch`, `#nested`, `#unexplained` | the choices of the three options described below |
| `#haplotype` | a panel haplotype's index and name, one line for each |
| `#note` | text describing the columns |
| `#H` | the column names of `H` rows |
| `#A` | the column names of `A` rows, in a version 6 file |

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
| `hap_index` | the panel haplotype's index in the `#haplotype` lines; `ref` on a row that follows the reference instead of a haplotype the strand copies; `*` on a row where no haplotype is named for the strand |
| `haplotype` | the panel haplotype's name, as its sample name and haplotype number joined by `#`; for a `ref` row, the reference's name in the panel; `*` on a row that names no haplotype |
| `sites` | number of called sites the row covers, or `.` for a `ref` row that fills a gap; a `ref` row that replaces a segment keeps its count |
| `gbwt_offset` | with `start_node`, a position in the graph's GBWT from which the haplotype can be followed to `end_node`; `.` if there is none |

Within a mosaic fragment, each row's `end_node` is the next row's `start_node`, so the fragment
expands to one walk in the graph, counting each shared node once. Each row lies within one of the
paths that store its panel haplotype, so a segment over several such paths gives several rows.
`gbwt_offset` is valid only for the graph named in `#graph`, and `hap_index` only within one file;
`haplotype` is the name to compare across files. A contig whose phased sites are all haploid has
rows for strand 0 only. In a haploid region of a diploid contig, strand 1 names no haplotype.

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
[unknown allele](read-likelihood-linkage-model.md#mutation-and-escape), so a segment's haplotype
need not cover the whole segment. Such a segment is replaced by a `ref` row where the reference
crosses it, and is otherwise a row that cannot be walked. A strand can also carry another allele
than its haplotype at a site, as a
[mutation](read-likelihood-linkage-model.md#mutation-and-escape), and its segment then names that
haplotype, whose walk spells the haplotype's allele there rather than the strand's.
`--mosaic-exact` adds the rows that give the strand's own allele at such sites (see
[Called alleles](#called-alleles---mosaic-exact)).

Three options change how the rows are formed, and the header records each choice:

- `--no-mosaic-patch-gaps` turns off both uses of the reference, so gaps stay unfilled, and a new
  mosaic fragment starts after each.
- `--no-mosaic-nested` leaves nested sites out of the segments, so a strand follows its enclosing
  site's haplotype through them.
- Where the linkage model names no haplotype for a strand, as for a nested haploid chain whose
  strand it could not determine, those sites are left out of the segments by default, and the
  walk crosses them by the rule above. `--mosaic-break-unexplained` writes a `*` row for them
  instead.

#### Called alleles (`--mosaic-exact`)

The walk the `H` rows expand to takes each panel haplotype's own allele at every site. The strand
was called with another allele where it carries a
[mutation](read-likelihood-linkage-model.md#mutation-and-escape) of its haplotype, where its
haplotype is an [unknown allele](read-likelihood-linkage-model.md#mutation-and-escape), where the
allele is one [no panel haplotype carries](read-likelihood-direct-genotyping.md), and at a site
with no record that lies on a segment of a haplotype with another allele. With `--mosaic-exact` the
file also has rows that start with `A`, which put the strand's own alleles on that walk, and its
version is 6. The `H` rows are the same with and without the option.

An `A` row gives one stretch of one strand. There are two kinds:

- A row for a top-level site gives the strand from the site's start bound node to its end bound
  node. It is the strand that the site's record spells: the allele the strand was called with,
  with every site nested in it spelled by the allele the strand was called with there, with or
  without a record of its own. With [`--flatten-records`](#flat-records), it is the strand whose
  differences from the reference the flat records of the site's run are, between those two
  nodes.
- A row for one node of the top-level chain gives the node's bases with the changes that the
  strand's alleles make to them at the [edit sites](read-likelihood-direct-genotyping.md) on the
  chain. Where two of those changes take one base of the node, the one further along the node is
  made, as in the records.

| Column | Meaning |
|---|---|
| `contig`, `strand` | as in an `H` row |
| `fragment` | the mosaic fragment whose walk the row is put on |
| `start_node`, `end_node` | as oriented node IDs: the site's two bound nodes, in the order the strand passes them, or the one node twice |
| `walk` | the strand's walk from `start_node` to `end_node`, both included: `>` and a node ID for a node read forwards, `<` and a node ID for one read backwards |
| `bases` | `.` where `walk` is a walk of the graph and every step of it has the bases of its node; otherwise one entry for each step of `walk`, separated by commas |

An entry of `bases` gives the bases the strand has at its step, read along `walk`. Within an
entry, `*` stands for the bases the visit already has, so `*` alone leaves the step as it is. An
entry can end with `+` and bases, or with `+` alone: those bases, or none, lie after the step on no
node, and the next step follows them without an edge of the graph. So `*+ACG` leaves the step as
it is and puts `ACG` between it and the next step, and `*+` joins the two steps directly. The row
of an [allele built from the reads](#alleles-built-from-the-reads---read-consensus) that is
given as its bases between the bounds has the two bound nodes alone as its `walk`, and the first
entry of its `bases` has the bases between them.

A site's row never changes the bases of its two bound nodes: their entries start with `*`. A
change to a bound node is in that node's own row. So no two `A` rows of a strand change the same
visit of a node, and the rows can be applied in any order.

To rebuild a strand:

1. Expand each fragment's `H` rows to walks, and join them on the nodes they share.
2. For each `A` row of the fragment, find each visit of `start_node` on that walk, and the next
   visit of `end_node` after it, which is the same visit for a row of one node.
3. Replace the visits between the two with the steps of `walk` between its first and last.
4. Give each step of `walk` its entry of `bases`.

A strand's `A` rows follow its `H` rows, in the order of the walk. A row is written only where
it changes the walk of the `H` rows: for a site, where the strand's walk between the bound nodes
is another one, or a step has an entry other than `*`; for a node, where the strand's alleles
change its bases.

Three things limit what the rebuilt strand says, and `vg call` reports the count of each:

- A site gets no row where the strand's `H` rows do not join its two bound nodes in one walk, or
  join them from end to start. The `H` rows run from a strand's first site with a record to its
  last, so a site before or after them is one of these.
- A top-level site with no phase gets no row, and neither do the sites nested in it.
- A change of an edit site gets no row where the strand's walk does not visit the node, or
  visits it between the bound nodes of a top-level site.

`walk` is found by its bound nodes, so a strand that crosses one site several times in one
fragment gets the same walk at each crossing.

`artifacts/scripts/mosaic_to_path.py` rebuilds the strands of a mosaic of either version from the
graph's paths as above. It checks that every `A` row has a stretch to go on, that no two rows
change the same visit, and, with `--gfa`, that every step of a rebuilt strand is an edge of the
graph or is joined by a row's `+`. With `--vcf` it also compares the bases of each strand with
the records:

- For each site with records, top-level or nested, the reference from the first base of the start
  bound node to the last base of the end bound node, with the strand's allele of each of the
  site's records put in, must be what the strand has there. The two bound nodes count with their
  own bases. A site written as [flat records](#flat-records) or as block records is compared with
  all its records together.
- For the nodes that the records of edit sites change, the reference over those nodes, with the
  strand's alleles put in, must be what the strand's visits of those nodes have.
- Where an `A` row makes a strand differ from the reference, a record must lie there.

The script prints the number of strands, nodes and stretches it does not compare, by reason. A
site with the records of other sites inside it is not compared, since its record spells each
allele once, with the alleles that one strand has at the sites inside; the sites inside are
compared. Records with no phase, strands with no allele at a record, sites and nodes the strand
does not cross, and two records of edit sites that take one base are the other reasons.

### Assembly anchors (`--anchors-out`)

The anchor file records, for each genotyped site, which reads support which of the sample's
strands, for use in pangenome-guided assembly. It also records the reads at points spaced along
runs of nodes with no branch, the [fill anchors](#fill-anchors). `--anchors-out` also turns on the genotyping of
[off-reference chains](#which-child-chains-are-genotyped), which can change the VCF's phases
and, under `--regenotype`, its genotypes; `--no-off-ref-nesting` turns that off. For a graph
sampled with its reference paths kept, map the reads as described under [Mapping without the
reference](#mapping-without-the-reference).

A **pin** is a point between two adjacent bases of the graph, with no sequence of its own. Each site
has two. The **start pin** lies just after the site's start boundary node, where an allele's walk
enters the interior. The **end pin** lies just before its end boundary node, where the walk leaves.

An off-panel site, whose ID begins `op`, lies inside nodes, and its pins bracket the bases its
alleles change: the start pin lies just after the last unchanged base before them, and the end pin
just before the first unchanged base after them, an indel's whole repeat tract counting as changed.
Where the changed bases end the site's walk, the end pin lies just after the last of them. An
insertion with no repeat tract changes no base, and its one pin is the point where it goes.
Off-panel changes on one node that the call bar admits whose changed bases
touch or overlap are one site where the reads spell alleles that carry each of them and they are not
all indels at one anchor, with ID `op<node>.<offset>.M<length>`: `offset` is its first
changed base along the node's forward strand and `length` the bases it replaces, and its ALTs are
the spellings of those bases that enough reads give, as the call bar counts them.

No point of a read alignment is in two anchors. A point of a read is the number of read bases
before it, and two alignments of one read name, such as paired-end mates, are told apart. Where an
off-panel site's anchor would hold a point that a snarl's anchor holds, or that an off-panel anchor
earlier in the file holds, it leaves that read out, and an off-panel anchor left with fewer than
`--anchors-reads` reads is not written. A [fill anchor](#fill-anchors) gives way to both kinds of
site anchor in the same way. vg reports on standard error how many read placements and anchors it
left out.

A site's reads are divided among its **slots**. A slot is a strand. Slot 0 is the strand of the
first allele of `GT`, to the left of the `|`, and slot 1 the strand of the second. But a diploid
site that is not phaseable pools both strands in one slot, 0, unless `--anchors-hom-split` divides
its reads between slots 0 and 1, which then carry the same allele. The alleles are those of the
site's line: two called traversals that spell one sequence are one VCF allele, so a site called
with both is homozygous, and a called traversal that is symbolically the reference is allele 0.
Each read's $p_{ra}$ of such a traversal is added to that of the traversal its allele is spelled
from. Alleles that `-L` merges as similar are not pooled. A site with no line keeps its called
traversals. A haploid site also has a single
slot, 0, except at a site of a nested haploid chain whose `GT` is `.|a`, where it is 1. An
**anchor** is one pin together with the reads of one slot that cross it.

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
  correction](read-likelihood-read-phasing.md#likelihood-correction). At a site in a
  [run](read-likelihood-read-phasing.md#runs), $y_{rs}$ comes from $\Lambda_{rs}$ as [Strand
  log-odds](read-likelihood-read-phasing.md#strand-log-odds) defines it, and a read whose strand
  comes from another run has its $y_{rs}$ replaced by
  $\mathrm{logit}(p \sigma(y_{rs}) + (1 - p) \sigma(-y_{rs}))$, with $p$ the probability that read
  phasing orients the two runs as the reads would. At a site in no run, $y_{rs}$ comes from the
  read's whole sum over its other sites of the phase set, whatever their runs. The read's anchor
  confidence is still computed from the $x_{ra}$. `--no-anchors-phase-hets` chooses the slot from
  the $x_{ra}$ alone.
- `--anchors-hom-split`, which needs `--read-phasing`, divides the reads of a diploid site that is
  not phaseable between two slots by the sign of their $y_{rs}$. A site is split only when, for each
  sign, at least `--split-min-side` reads have a $y_{rs}$ of that sign and of absolute value at
  least `--split-min-q`. A split site leaves out a read whose strand is not usable there, because
  the read was used in another phase set or in more than one: its strand belongs to another phase
  set's labelling. Any other read whose $y_{rs}$ is 0 is placed by a coin flip derived from its
  name, so it takes the same slot at every site.

#### Fill anchors

A site's pins lie at its boundaries, so sequence that no site varies, and the inside of a site
whose alleles are long, can lie far from every pin. **Fill anchors** pin points that the graph
alone chooses, so that such sequence is pinned too. `--anchors-fill` sets their spacing $K$, 1000
bp by default, and 0 turns them off.

A **unitig** is a longest run of nodes in which each node's right side has one edge, to the next
node of the run, and that node's left side has that one edge alone. The two nodes of each such join
differ. Inside each unitig of length $L$ greater than $K$, $n = \lceil L/K \rceil - 1$ **fill
points** divide the length evenly: the $i$-th point has $iL/(n+1)$ bases of the unitig before it,
rounded to the nearest base, half up. A unitig is read from its end whose node has the smaller ID. A
unitig of one node is read along the node's forward strand, and so is one that closes on itself,
starting from its node with the smallest ID. A point that falls on a join between two nodes moves
forward one base at a time until it lies inside a node. It is dropped if it reaches the end of the
unitig, and a later point that the move reaches or passes is dropped too. The points depend only on
the graph and $K$, so every run on one graph uses the same ones.

Each read that crosses a fill point is placed there, as at a site's pin, and the reads at a point
make one **fill anchor**, which has one slot. The reads pass the same filters as a site's reads:
`--read-min-mapq`, and no secondary alignments. A read is placed only where it visits the point's
node once, has a read base aligned to the node base before the point, along the node's forward
strand, and has a read base on the point's other side outside its soft clips. A read's anchor
confidence at a fill anchor is $-10 \log_{10} e_r$, which is its confidence at a site with one
allele.

The slot and kind of a fill anchor come from the called walks of its **deciding site**: the
innermost site that holds the point's node and has a called genotype. A nested chain of one copy
whose strand is not known, as where the linkage model gave its parent no phase, is passed over, so
the site around it decides. So is a site that was not genotyped, as under `--no-nested` or below a
nested site that could not be genotyped.

- **fill-shared**, in slot 0, where both of the sample's strands cross the point: in sequence that
  both called alleles of the deciding site share, at a homozygous site, at a site in a region of
  ploidy 1, or where no site holds the node.
- **fill-strand**, in slot $s$, where only the called walk of strand $s$ crosses it: on a node that
  only one called allele visits, or in a nested haploid chain on strand $s$.
- **fill-off-call**, in slot 0, where no called walk of the deciding site crosses it, or where sites
  hold the node but none has a called genotype.

A fill-strand slot is field $s$ of the deciding site's phased `GT`, as a site anchor's slot is. The
`A` row does not name the deciding site, so the slot names a haplotype only where a reader finds
that site, from the snarl decomposition, and it is phased.

A fill anchor gives way at a read point to every snarl's and off-panel site's anchor, and to a fill
anchor written before it. It is written only if at least `--anchors-reads` of its reads remain,
whatever its kind. `--anchors-min-q` applies to its reads. `--anchors-het-only`,
`--anchors-leaf-only` and `--anchors-min-gqn` choose sites, so they do not apply to it.

A fill anchor is only as right as the read placements it rests on. In centromeric arrays and other
large repeats, a read from one copy of the repeat can be placed on another copy in the panel. A
site anchor leaves out a read whose best allele is not one the site was called with, which drops
some of those reads, but a fill anchor keeps every read that crosses its point, so in such repeats a
fill anchor holds misplaced reads more often than the site anchors around it. Where sequence has
close copies elsewhere in the genome, reads from those copies can be placed on anchors of every
kind. Better graphs or haplotype sampling are the remedy, not a rule about anchors.

#### Rows of the anchor file

The file is tab-separated. Each anchor is written as an `A` row, followed by one `R` row for each
of its reads. Header lines start with `#`, and a reader should skip keys it does not recognise.
Among them, `#anchors-version` gives the format version, 10, which changes whenever the columns or
their meaning change. The `#filters` line gives the thresholds the anchors were written with,
`fill=` the spacing of the fill anchors among them. The `#read` lines list the read names, which
`R` rows refer to by index. The `#H` lines name the columns of each kind of row, and the `#note`
lines describe them.

| `A` column | Meaning |
|---|---|
| `node` | the graph's ID of the node the pin lies on. For a snarl it is the pin's boundary node: the first node of `snarl` for the start pin, the second for the end pin. Under `--translation` or `--gbz-translation`, `snarl` names the translated segments, while `node` stays the graph's ID |
| `snarl` | the site's ID, as in the VCF `ID` column: its start and end boundary node IDs, each preceded by `>` or `<` for its orientation. It is also written for off-reference sites, which have no VCF record. A fill anchor's is `fill<node>.<pin>` |
| `slot` | the slot |
| `allele` | the slot's allele, as its index in the site's list of candidate alleles rather than a VCF allele number: at a site with a line, the candidate its VCF allele is spelled from. The list is not written, so the index serves to compare a site's slots. `.` for a fill anchor |
| `gqn` | the `GQN` that a whole-site record for the site would carry, moved or not; `.` where it has none, and for a fill anchor |
| `explained` | the site's explained share; `.` for a fill anchor |
| `reliability` | the site's anchor reliability, counting the reads placed at an off-panel site before any was left out for a point another anchor holds; `.` for a fill anchor |
| `pin` | where the pin lies in `node`: the node's bases before it, along the node's forward strand |
| `kind` | `het` for a strand of a genotype with two different alleles, `hom` for the one slot of a homozygous genotype, `hom-split` for a slot of a homozygous genotype that `--anchors-hom-split` divided, `haploid` for the one slot of a genotype of one allele (a nested haploid chain, or a site in a region of ploidy 1), and `fill-shared`, `fill-strand` or `fill-off-call` for a fill anchor |

| `R` column | Meaning |
|---|---|
| `read_id` | the read's index in the `#read` lines; paired mates share a name, and so an index |
| `strand` | the direction in which the read crosses the pin, not one of the sample's strands: 0 in the site's direction, 1 against it. A fill anchor's direction is its node's forward strand |
| `offset` | the 0-based offset, in the read as sequenced, of the read's last base before the pin, reading in the site's direction. The point of the read the pin is at follows that base on strand 0 and precedes it on strand 1, so it has `offset` + 1 read bases before it on strand 0 and `offset` on strand 1 |
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
- `--anchors-fill` sets the spacing of the [fill anchors](#fill-anchors), and 0 writes none.
- With `--anchors-end-new` set to $N$, a slot's end anchor is written only if at least $N$ of the
  slot's reads cross the end pin but not the start pin, which drops end anchors that mostly repeat
  the start anchor's reads.

## Options

`--preset ont` sets several of the options below to values suited to Oxford Nanopore reads, and
`--preset hifi` to values suited to PacBio HiFi reads; `vg call --help` lists which, and each
option's default. The `hifi` values are the `ont` ones with `--insertion-nats -1.8`, `--hp-prior 0`
and `--linkage-weight 1`, chosen by F1 against the HG002 benchmark on Revio reads aligned to an HPRC
graph sampled with the reads' k-mers. An option given explicitly overrides the preset. Where an
option has a `--no-` form, the two set the same thing and the one given last wins, so a `--no-` form
also turns off a setting that a preset turned on.

The table includes general `vg call` options: `--pack`, `--gbwt`, `--gbz`, `--max-snarl-edges`,
the Ploidy row, `--cluster` and `--cluster-min-len`. Its other options are rejected without `--read-likelihood`. The options
that modify `--anchors-out` (`--no-off-ref-nesting` among them), `--mosaic-out` and `--regenotype`
are also rejected when those are not in use. Other general options used on this page, such as
`--genotype-snarls`, `--sample`, `--snarls`, `--ref-path`, `--path-prefix`, `--ref-sample`,
`--all-snarls`, `--translation` and `--gbz-translation`, are described by `vg call --help`.
`--read-likelihood` cannot be combined with `--bottom-up`, `--top-down`, `--star-allele`,
`--chains`, `--gaf`, `--traversals`, `--vcf` or `--legacy`.

| Part | Options |
|---|---|
| [Reads](read-likelihood-direct-genotyping.md#read-input) | `--gam`, `--gaf-reads`, `--gam-index`, `--gaf-index`, `--gaf-base`, `--gaf-base-binary`, `--read-window`, `--read-min-mapq` |
| [Candidate alleles](#candidate-alleles) | `--enumerate-support`, `--pack`, `--gbwt`, `--gbz`, `--max-snarl-edges` |
| [Relative likelihood](read-likelihood-direct-genotyping.md#relative-likelihood) | `--gap-open`, `--gap-extend`, `--insertion-nats`, `--optimal-pairing`, `--no-optimal-pairing` |
| [Mismapping](read-likelihood-direct-genotyping.md#mismapping-probability) | `--mismap-min`, `--mismap-max`, `--no-mismap-term` |
| [Mixture weights](read-likelihood-direct-genotyping.md#mixture-weights) | `--flat-mixture` |
| [Depth term](read-likelihood-direct-genotyping.md#depth-term) | `--depth-term`, `--depth-count-raw` |
| [Alleles built from the reads](#alleles-built-from-the-reads---read-consensus) | `--read-consensus`, `--consensus-min-gain` |
| [Linkage](read-likelihood-linkage-model.md) | `--linkage-weight`, `--linkage-scale`, `--linkage-prior`, `--linkage-mutation`, `--hp-prior`, `--hp-prior-run` |
| [Nested sites](#nested-sites) | `--nested`, `--no-nested`, `--group-candidates`, `--no-group-candidates`, `--atomize-blocks`, `--no-atomize-blocks`, `--flatten-records`, `--flatten-distance`, `--no-off-ref-nesting` |
| [Ploidy](#ploidy) | `--ploidy`, `--ploidy-regex`, `--ploidy-bed`, `--ploidy-preset` |
| [Phasing](#phasing) | `--phased`, `--no-phased`, `--read-phasing`, `--no-read-phasing`, `--phase-min-q`, `--phase-coherence`, `--phase-coh-rounds`, `--phase-break`, `--phase-relink`, `--phase-hang`, `--phase-prior`, `--phase-cap` |
| [Re-genotyping](#re-genotyping-from-the-phase) | `--regenotype`, `--no-regenotype`, `--regeno-temper`, `--regeno-ceiling`, `--regeno-passes`, `--regeno-haploid`, `--no-regeno-haploid`, `--regeno-ledger` |
| [Quality fields](#options-that-change-the-fields) | `--no-share-quality`, `--depth-quality`, `--min-confidence` |
| [Allele merging](#options-that-change-the-fields) | `--cluster`, `--cluster-min-len` |
| [Mosaic](#mosaic---mosaic-out) | `--mosaic-out`, `--mosaic-patch-gaps`, `--no-mosaic-patch-gaps`, `--no-mosaic-nested`, `--mosaic-break-unexplained`, `--mosaic-exact` |
| [Anchors](#assembly-anchors---anchors-out) | `--anchors-out`, `--anchors-fill`, `--anchors-reads`, `--anchors-min-gqn`, `--anchors-min-q`, `--anchors-het-only`, `--anchors-leaf-only`, `--anchors-keep-off-call`, `--anchors-end-new`, `--anchors-phase-hets`, `--no-anchors-phase-hets`, `--anchors-hom-split`, `--split-min-q`, `--split-min-side` |
| Presets | `--preset` |

### Fixed constants

These constants have no option. Each is located by its name, or by the function or member that
holds it.

| Constant | Defined in | What it sets |
|---|---|---|
| `LinkageModel::Params::escape` | `artifacts/linkage_model.hpp` | the escape penalty $\epsilon_{\mathrm{esc}}$ on a strand with an [unknown allele](read-likelihood-linkage-model.md#mutation-and-escape) |
| `LinkageModel::Params::rho_min` | `artifacts/linkage_model.hpp` | the floor $\rho_{\min}$ in the switch probability (see [Transitions](read-likelihood-linkage-model.md#transitions)) |
| `LinkageModel::Params::window`, `margin` | `artifacts/linkage_model.hpp` | the number of sites each forward–backward window keeps, and the extra sites decoded and discarded on each side |
| `RATE_BUCKET`, `RATE_ID_WINDOW` | `artifacts/allele_likelihood.hpp` | reference bp per rate-window bucket (a rate window is three buckets), and node IDs per rate window when there are no reference positions |
| `banded` and `band`, in `score_by_optimal_pairing` | `artifacts/allele_likelihood.cpp` | the grid size, in read visits times allele visits, above which optimal pairing is restricted to a band, and the band's half-width in allele visits |
| `AlleleLikelihoodParams::sequence_allele_max_bases`, `sequence_allele_max_alignments` | `artifacts/allele_likelihood.hpp` | the most bases of a read between a snarl's bounds that are aligned, and the most alignments of a read to an allele at one snarl, at a snarl with an [allele built from the reads](#alleles-built-from-the-reads---read-consensus) |
| `ReadConsensusParams::min_reads`, `max_reads`, `max_rounds`, `max_medoid` | `artifacts/read_consensus.hpp` | for a consensus of reads: the fewest and the most reads it is built from, the most times they are aligned to it, and the most inserted strings compared to choose the bases of an insertion |
| `ReadConsensusOptions::run_max_junction` | `artifacts/read_consensus_step.hpp` | for a [run](#runs-of-neighbouring-snarls): the most bases in a junction, and between two seeds that follow each other |
| `max_yens_traversals` | `src/subcommand/call_main.cpp` | the most candidate alleles support enumeration keeps |
| the difference limit in `LinkageModel::run_length_site` | `artifacts/linkage_model.cpp` | the largest homopolymer length difference to which `--hp-prior` applies |
| the `int8_t` of `LinkageCollector::allele_arena` | `artifacts/linkage_model.hpp` | the most alleles a site's compact allele set may hold for the linkage model |
| the 64-bit crossing mask of `crossing_mask` | `artifacts/child_crossing_finder.hpp` | the most candidate alleles a parent site can have for vg to tell which of them cross a child site |
| the `min_length` argument of `set_depth_quality` | `artifacts/read_likelihood_caller.hpp` | the change in allele length at which `--depth-quality` applies |
| the minimum counted reads in `read_phase_flips` | `artifacts/read_phasing.cpp` | the reads a site needs before low coherence can remove it from the phase chain |
| `RegenotypeParams::fit_bins`, `fit_min_per_bin`, and the grid in `fit_calibration` | `artifacts/regenotype.hpp`, `artifacts/regenotype.cpp` | the bins, the fewest observations per bin, and the candidate values used to fit $\tau$ |
| `PANEL_FIT_MIN_PAIRS`, `PANEL_FIT_BINS`, and the grid in `fit_panel_errors` | `artifacts/read_phasing.cpp` | the fewest pairs of sites the panel's error rates $\lambda$ and $f$ are fitted to, the bins of the fit's report, and the candidate rates |

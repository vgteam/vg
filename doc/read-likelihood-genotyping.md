# Read-likelihood genotyping (`vg call --read-likelihood`)

`vg call --read-likelihood` genotypes one sample from its reads aligned to a pangenome graph.
Default values of the options named below are listed by `vg call --help`; the few constants that
have no option are listed under [Fixed constants](#fixed-constants). The source files that
implement the method, and an order in which to read them, are described in
[read-likelihood-architecture.md](read-likelihood-architecture.md).

A minimal run takes a graph in
[GBZ](https://github.com/vgteam/vg/wiki/Extra-details-on-vg-file-formats#gbz-gbwtgraph--gbz) format
and reads aligned to it by `vg giraffe` in
[GAM](https://github.com/vgteam/vg/wiki/Extra-details-on-vg-file-formats#gam-graph-alignment--map-vgs-bam)
format, and writes the calls in [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format:

```
vg call graph.gbz --read-likelihood --gam reads.gam > calls.vcf
```

`vg call` finds the graph's sites itself. It can instead read them, with `-r`, from a file of
[snarls](https://github.com/vgteam/vg/wiki/Snarls-and-chains#snarls) made by `vg snarls`.

Each site goes through four steps:

1. **Site likelihood computation.** For each genotype the site could have, vg computes a
   likelihood: a number that measures how well the genotype explains the site's reads.
2. **Genotyping.** vg chooses the site's genotype from its likelihoods. It can take the genotype
   with the highest likelihood, the [direct call](#direct-call-or-linkage), or use a *linkage
   model*. The linkage model also uses the haplotypes stored in the graph, which tend to carry the
   same combinations of alleles at neighbouring sites as the sample. With
   [nested calling](#nested-sites), on by default, sites nested inside other sites are genotyped
   too, and get VCF records of their own.
3. **Phasing.** vg assigns each genotype's alleles to the sample's haplotypes, first from the
   stored haplotypes and then, optionally, from reads that span several sites. With
   `--regenotype`, the phase is then used to correct the likelihoods, and genotyping and phasing
   are repeated.
4. **Output.** vg writes a VCF file. It can also write a *mosaic* file, which describes each of
   the sample's haplotypes as a walk through the graph, and an
   [anchor file](#assembly-anchors---anchors-out), which ties reads to haplotypes.

The linkage model decides many neighbouring sites together, and a nested site's genotype depends
on its parent's. So when the linkage model or nested calling is on, as both are by default, vg
first computes the likelihoods of every site. It then genotypes and phases the sites one
[generation](#nested-sites) at a time, parents before children.

## Vocabulary

### Graph

- **Node.** A node holds a DNA sequence and can be traversed in either *orientation*: forward,
  reading its sequence, or reverse, reading its reverse complement. A *walk* is a sequence of
  oriented node visits along the graph's edges.
- **Site.** A place in the graph where the sample's genome may differ from the reference. vg
  identifies sites as *snarls*. A snarl is a subgraph separated from the rest of the graph by two
  *boundary nodes*, a start and an end. On this page a site is a snarl that `vg call` genotypes.
  Its boundary nodes belong to it, and its other nodes, including those of any snarls nested in
  it, are its *interior*. A site is oriented from its start boundary node to its end boundary
  node, and its alleles are written in that direction.
- **Chain.** A series of snarls joined end to end, each snarl's end boundary node being the next
  one's start. Snarls nest: a snarl can contain chains of smaller snarls, its *child chains*. The
  sites in a site's child chains are *nested* in it, and it is their *parent*.
- **Allele.** A walk through a site from its start boundary node to its end boundary node, also
  called a *traversal*.
- **Reference.** The *reference paths* are the graph paths on which `vg call` reports positions.
  By default they are the paths that the graph's metadata marks as reference or generic paths;
  `-p`, `-P` or `-S` selects other paths instead, by name, by name prefix or by sample. A site's
  *reference allele* is the walk a reference path takes through it. In the VCF the reference
  allele is allele 0, and the other alleles are numbered from 1.

### Sample and reads

- **Haplotype.** One copy of a chromosome, or of part of one. A haplotype that passes through a
  site carries one allele there.
- **Ploidy.** A site's ploidy, $P$, is the number of the sample's haplotypes that pass through
  the site. Options set the ploidy of each contig or region, and a nested site's ploidy follows
  from its parent's genotype (see [Ploidy](#ploidy)). The ploidy options accept only 1 or 2, with
  or without `--read-likelihood`, because the linkage model, phasing and nested calling are built
  for at most two haplotypes. A nested site that none of its parent's alleles passes through has
  ploidy 0 (see
  [Which child chains are genotyped](#which-child-chains-are-genotyped)).
- **Panel.** The haplotypes stored in the GBZ file's
  [GBWT](https://github.com/vgteam/vg/wiki/Extra-details-on-vg-file-formats#gbwt-gbwt) index, or in
  a separate GBWT file given with `-g`. A panel haplotype can be stored as several paths that each
  cover part of a chromosome, so it does not necessarily pass through every site.
- **Genotype.** The multiset of the $P$ alleles that the sample's haplotypes carry at a site, such
  as $\lbrace 0, 1 \rbrace$. Which haplotype carries which allele is the genotype's *phase*, and
  is decided separately.
- **Read placements.** A mapper, such as `vg giraffe`, aligns the reads to the graph. The reads'
  *placements* are the walks their alignments take, together with the *edits* inside each node
  (the runs of matching bases, mismatches, insertions and deletions against the node's sequence)
  and their mapping qualities (MAPQ).
- **Reads of a site.** The reads whose placements visit an interior node of the site, or visit
  both of its boundary nodes. A read that visits both boundary nodes and no interior node supports
  an allele that skips the interior, such as a deletion. A read whose only node in the site is one
  boundary node fits every allele equally, and is not used.

## Site likelihood computation

The site likelihood $\mathcal{L}(G)$ measures how well a genotype $G$ explains the reads of one
site. It is built from a model of how a site's reads arise from the sample's genotype. The model
can be written as a formula in a few quantities, and each of those quantities can be computed from
the read placements.

### Model of how reads arise

We treat the reads of a site as generated from the sample's genotype $G$ by this process:

1. **How many reads.** The number of reads that come from the site's haplotypes is Poisson
   distributed, with a mean that grows with the length of $G$'s alleles. Reads from elsewhere in
   the genome can also be placed at the site by mistake; such a read is *mismapped*.
2. **Where each read came from.** Independently for each read placed at the site, with a
   probability set from its MAPQ, the read is mismapped, and then its bases say nothing about $G$.
   Otherwise it came from one of $G$'s haplotypes.
3. **What the read says.** The part of the read inside the site is a copy of part of the allele
   on this haplotype for the site, with sequencing errors from an error model.

The model takes each read's MAPQ as given, and does not model how many reads are mismapped.

Each step leaves a different kind of evidence about $G$ in the reads. The number of reads reflects
the alleles' lengths: a homozygous deletion shows mainly as reads that are missing. The chance of
mismapping means that a read with a doubtful placement counts for less. Most of the evidence comes
from the copying step, because a read fits the allele it was copied from better than it fits the
others.

The model assumes four things that real data do not always satisfy:

- Reads are independent given the genotype, although paired mates, and reads that share a
  systematic error, are not.
- Whether a read is mismapped does not depend on how well it fits the site's alleles, so a read
  that fits every allele badly is no more likely to be mismapped than one that fits well.
- The number of reads depends only on the genotype, and not on mappability or base composition.
- Each haplotype passes through the site once. A duplication can make a haplotype's path pass
  through the site twice. The model still gives that haplotype one allele, and explains the reads
  of both copies with it, so the site has more reads than $G$ predicts. A walk that loops inside
  the site, between its boundary nodes, is different: the whole loop is part of one allele.

vg also does not align the reads' bases again. The probability of a read given an allele comes
from the alignment that the mapper made (see [Relative likelihood](#relative-likelihood)).

### Notation

| Symbol | Meaning |
|---|---|
| $A$ | the site's candidate alleles: the alleles vg considers at the site (see [Candidate alleles](#candidate-alleles)) |
| $a, b$ | alleles in $A$ |
| $P$ | the ploidy at the site |
| $G = \lbrace g_1, \dots, g_P \rbrace$ | a genotype; $g_1, \dots, g_P$ list its alleles in any order, and nothing below depends on that order. A homozygous genotype lists the same allele $P$ times |
| $\mathcal{L}(G)$ | the site likelihood of $G$ |
| $R$ | the reads of the site |
| $r$ | a single read in $R$ |
| $\Pr(r \mid a)$ | the probability of read $r$'s bases, if it was copied from allele $a$, starting where its placement puts it |
| $\Pr(r \mid \text{mismapped})$ | the probability of read $r$'s bases, if it was mismapped |
| $p_{ra}$ | the relative likelihood of read $r$ under allele $a$, defined below |
| $e_r$ | the probability that read $r$ is mismapped |
| $\epsilon_{\min}, \epsilon_{\max}$ | the floor and ceiling on $e_r$ |
| $w_i(G)$ | the mixture weight of the haplotype carrying $g_i$ |
| $f_{\mathrm{Pois}}(n ; \lambda)$ | $\lambda^{n} \exp(-\lambda) / \Gamma(n + 1)$ for a real $n \geq 0$, where the gamma function $\Gamma$ extends the factorial, $\Gamma(n + 1) = n!$; for a whole number $n$ it is the Poisson probability of $n$ events when $\lambda$ are expected |
| $N_{\mathrm{eff}}$ | the *effective read count*, $\sum_{r \in R} (1 - e_r)$: the expected number of the site's reads that were not mismapped, given their MAPQs |
| $\mu_G$ | the *expected read count*: the number of reads from the site's haplotypes that $G$ predicts |
| $\beta$ | the exponent of the depth term, `--depth-term` |
| $T_a$ | the length in bases of allele $a$, excluding the site's two boundary nodes; a node the allele visits twice counts twice |
| $U_i(G)$ | at ploidy 2, the length in bases of the nodes that $g_i$ visits and the other allele of $G$ does not, each node counted once |
| $\bar L$ | the mean read length near the site |
| $\kappa$ | the rate of read starts near the site, each counted as $1 - e_r$, per base and per haplotype of the region |

### Likelihood formula

vg computes the site likelihood as

$$
\mathcal{L}(G) = \prod_{r \in R} \left[ (1 - e_r) \sum_{i=1}^{P} w_i(G) p_{r g_i} + e_r \right] \times f_{\mathrm{Pois}}\left(N_{\mathrm{eff}} ; \mu_G\right)^{\beta}
$$

The product over reads, the *read term*, follows steps 2 and 3 of the model. The last factor, the
*depth term*, follows step 1. Neither is exactly the model's probability, and each is explained in
turn below, with the places where it departs from the model. vg works with $\ln \mathcal{L}(G)$,
the sum of the logarithms of the factors.

#### Read term

Take one read, and suppose its start on each haplotype is where its placement puts it. With
probability $1 - e_r$ the read was not mismapped, and came from one of $G$'s haplotypes. vg gives
the haplotype carrying $g_i$ a weight $w_i(G)$, its *mixture weight*, and in that case the read's
probability is $\Pr(r \mid g_i)$. With probability $e_r$ the read was mismapped, and its
probability is $\Pr(r \mid \text{mismapped})$. The read's probability is the sum of these cases:

$$
(1 - e_r) \sum_{i=1}^{P} w_i(G) \Pr(r \mid g_i) + e_r \Pr(r \mid \text{mismapped})
$$

The model does not say what $\Pr(r \mid \text{mismapped})$ is, since a mismapped read came from
somewhere else in the genome. vg uses the probability of the read under its best allele at the
site, $\max_{b \in A} \Pr(r \mid b)$. This value is the same under every genotype, so the
mismapped case favours no genotype.

Dividing one read's probability by a number that is the same under every genotype divides every
genotype's likelihood by that number. It changes neither which genotype is most likely nor the
ratio of two genotypes' likelihoods. vg divides each read's probability by
$\max_{b \in A} \Pr(r \mid b)$. The result depends on the read only through $e_r$ and the read's
*relative likelihoods*,

$$
p_{ra} = \frac{\Pr(r \mid a)}{\max_{b \in A} \Pr(r \mid b)}
$$

and the mismapped case becomes $e_r$. A relative likelihood is 1 for the read's best allele, and
smaller for alleles that explain the read less well. Each read's factor in the read term therefore
lies between $e_r$ and 1: it is $1 - e_r$ times the average of the read's relative likelihoods
under $G$'s alleles, weighted by the mixture weights, plus $e_r$.

The read term departs from the model in two ways. First, vg chooses the mixture weights, and
$\Pr(r \mid a)$ leaves out the chance of the read's start position; the two are related, as
[Mixture weights](#mixture-weights) explains. Second, each read's factor sums over its own two
cases, mismapped or not, as if that did not affect how many reads came from the site's haplotypes.
The model's exact probability would sum over every assignment of the reads to the two cases at
once, together with the count.

#### Depth term

The depth term asks whether the number of reads fits $G$. Which reads came from the site's
haplotypes is not known, because any read may have been mismapped, so vg uses the effective read
count $N_{\mathrm{eff}}$, which counts each read by its probability of not being mismapped. That
count is usually not a whole number, and the Poisson distribution gives probabilities only to
whole numbers. vg therefore models $N_{\mathrm{eff}}$ with a continuous analogue of the Poisson
distribution: a distribution on $n \geq 0$ whose density is proportional to
$f_{\mathrm{Pois}}(n ; \lambda)$.

The depth term evaluates this function at $N_{\mathrm{eff}}$, as a likelihood of the expected read
count $\lambda = \mu_G$. The density itself would be the function divided by its integral over
$n$. That integral depends on $\lambda$, so leaving it out changes the ratio of two genotypes'
likelihoods, but by less than 1% when both expected read counts are at least 4 and $\beta$ is at
most 1. vg leaves it out.

The depth term is also raised to the power $\beta$, below 1 by default, because read depth varies
between sites for reasons the model leaves out.

With these departures, $\mathcal{L}(G)$ as a whole is not proportional to the model's probability
of the reads. It is a likelihood built from the model's parts.

The expected read count is

$$
\mu_G = \kappa \sum_{i=1}^{P} \left(T_{g_i} + \bar L - 1\right)
$$

$T_{g_i} + \bar L - 1$ is the number of positions at which a read of length $\bar L$ can start on
the haplotype and still include some of the interior of $g_i$. When $g_i$ has no interior
($T_{g_i} = 0$), it is the number of positions from which a read reaches across the junction
between the two boundary nodes.

### Computing the terms

The read term needs, for each read of the site, its relative likelihood under each allele and its
mismapping probability, and, for each genotype, the mixture weights. The depth term needs
$N_{\mathrm{eff}}$, which also comes from the reads of the site, and $\mu_G$, whose $\kappa$ and
$\bar L$ come from the reads near the site.

#### Read input

Reads with MAPQ below `--read-min-mapq`, secondary alignments, and alignments with no placement
are discarded as they are read in. The reads come from one of three sources:

- `--gam` or `--gaf-reads`: a GAM or [GAF](static/GAF.md) file, loaded into memory.
- `--gam` with `--gam-index`: a GAM file sorted and indexed by `vg gamsort -i`. Reads are fetched
  from it as they are needed, one range of consecutive node IDs at a time; node IDs in these
  graphs increase roughly along the genome.
- `--gaf-base`: a [GAF-base](https://github.com/jltsiren/gbz-base) database of alignments,
  queried one range of node IDs at a time by the `gbz-base` program. `--gaf-base-binary` gives the
  path to that program. It reads the graph from the GBZ-base database given with `--gbz-base`, or
  else from the input graph. This source also drops an alignment returned twice by its queries,
  identified by read name and start position.

`--read-window` sets the size of those ranges, in node IDs, for the two indexed sources. It changes
which reads are fetched together, but not which reads a site uses. It can change a likelihood in
its last digits, since the reads' terms are then added in another order.

A read whose alignment crosses the site against the direction of the alleles is
reverse-complemented before it is scored. The direction is decided by a vote over the nodes the
alleles visit, leaving out any node that two alleles visit in opposite orientations. The read is
reversed if more of its visits to those nodes are in the opposite orientation than in the same
one. A tie, including a read that visits none of them, leaves it as it is.

#### Relative likelihood

A read's relative likelihood under an allele comes from scoring the read against the allele. vg
pairs the read's node visits with the allele's, scores the pairing, and compares the score with the
score of the read's best allele. Two visits are the same when they go to the same node in the same
orientation.

A *pairing* lines up the read's node visits inside the site (boundary nodes included) with the
allele's node visits, in order, and may leave some of either unpaired. The rules for scoring a
pairing are:

| Read visit paired with allele visit | Score |
|---|---|
| the same visit | the score of the read's own edits in that node |
| a visit the allele never makes, paired with a visit the read never makes | a *substitution*: the read's bases in its node compared one by one with the allele node's sequence, from the first base of each, over the shorter length, plus a gap for the difference in length |
| any other pair | not allowed |

- A run of allele visits left unpaired between two pairs is a *deletion*, scored as one gap as
  long as their bases. Unpaired allele visits before the read's first same-visit pair, or after
  its last pair, lie outside the read and score nothing.
- A run of consecutive read visits left unpaired after the first same-visit pair is an
  *insertion*, scored as one gap as long as their bases. A visit that the read and the allele
  share can be left unpaired in this way. Before the first same-visit pair, each unpaired read
  visit is a gap of its own.

So every read base inside the site is scored under every allele, and all alleles are scored over
the same read bases.

Scores use vg's alignment scoring, in which a better fit scores higher. Matches and mismatches
score from a substitution matrix adjusted by base quality, in which a mismatch at a low-quality
base costs less; a read without base qualities uses the plain matrix. A gap of length $m$ scores
$-(o + (m - 1) x)$, where $o$ is `--gap-open` and $x$ is `--gap-extend`. The mapper's edits say
which bases match, mismatch, or are inserted or deleted, and their scores are computed here with
these settings. The total score $s_{ra}$ of the pairing is converted to a log-likelihood score in
nats:

$$
\ell_{ra} = \alpha s_{ra} + \iota I_{ra}
$$

$\alpha$, which vg's code calls the scorer's *log base*, is a scale factor that vg computes from
the substitution matrix, so that $\alpha$ times a substitution score is a natural-log likelihood
ratio. $I_{ra}$ is the number of insertions in which the read has bases the allele lacks, counting
those inside the mapper's edits, and $\iota$ is `--insertion-nats`. Under optimal pairing (below),
each unpaired read visit counts as one insertion here, even where a run of them is scored as one
gap. A positive $\iota$ raises $\ell_{ra}$ for each insertion, so extra read bases count against an
allele less than missing ones.

vg's substitution scores are log-odds of the read's bases against a random-sequence background.
Since every allele is scored over the same read bases, $\ell_{ra}$ stands for
$\ln \Pr(r \mid a)$ plus a term that depends on the read and not on the allele. That term cancels in
the relative likelihood:

$$
p_{ra} = \exp\left(\ell_{ra} - \max_{b \in A} \ell_{rb}\right)
$$

$\ell_{ra}$ is only an approximation to a log probability: its gap penalties are fixed costs, not
probabilities from an error model.

vg can choose the pairing in two ways. Optimal pairing finds the best pairing that the rules
allow, and greedy pairing, the default, approximates it. Inside a node that the read and the
allele share, both score the mapper's edits.

##### Greedy pairing

Greedy pairing works through the read's visits in order. For each one, it looks in the allele for
the same visit, after the last allele visit already paired. If it finds one, the two are paired,
and once a first such pair exists, any allele visits skipped over become a deletion. Otherwise, if
the allele's next unpaired visit occurs later in the read, the read visit is left unpaired as an
insertion. Otherwise the read visit is paired with the allele's next unpaired visit as a
substitution. Once the allele's visits are used up, the remaining read visits are insertions. A
pair, once made, is never revised.

Greedy pairing therefore departs from the rules in two ways. It can make a substitution pair that
the rules forbid, and it scores each unpaired read visit as a gap of its own, rather than one gap
per run.

##### Optimal pairing

With `--realign`, optimal pairing finds the pairing that the rules allow with the highest score
$s_{ra}$, by dynamic programming over the read's visits against the allele's, and then adds
$\iota$ times that pairing's insertions. `--realign` chooses the pairing of node visits again; it
does not align bases.

When the read and the allele both have many visits, so that the product of their numbers exceeds a
fixed limit, the search is restricted to a band. For each read visit, the band is centred where
the next visit that the read and the allele share (or, past the last one, the last) puts the
matching allele visit, and it holds the allele visits within a fixed distance of that centre. At
the largest sites, optimal pairing is therefore approximate too.

#### Mismapping probability

MAPQ is the mapper's estimate, on the phred scale, of the probability that the read belongs
somewhere else. vg uses it as the mismapping probability, held between a floor and a ceiling:

$$
e_r = \min\left(\max\left(10^{-\mathrm{MAPQ}_r / 10}, \epsilon_{\min}\right), \epsilon_{\max}\right)
$$

$\epsilon_{\min}$ is `--mismap-min` and $\epsilon_{\max}$ is `--mismap-max`. Each read's factor in
the read term lies between $e_r$ and 1. So one read changes the ratio of two genotypes' read terms
by at most a factor of $1 / e_r$, which is at most $1 / \epsilon_{\min}$. The floor therefore sets
how strongly a single read's fit can count against a genotype. (The read also adds $1 - e_r$ to
$N_{\mathrm{eff}}$, which moves the depth term.)

MAPQ does not cover one kind of error: a read placed at the right locus can still be aligned
wrongly through this particular site. The floor also stands for that case.

The ceiling applies to reads with MAPQ 0 or close to it, whose unclamped $e_r$ is near 1, and keeps
them contributing a little rather than nothing. Such reads are used unless `--read-min-mapq`
excludes them. `--no-mismap-term` sets every $e_r$ to $\epsilon_{\min}$, in $N_{\mathrm{eff}}$ and
$\kappa$ as well as in the read term.

#### Mixture weights

vg chooses the mixture weights; they do not follow from the model. In the model, the chance that a
read came from haplotype $i$ is proportional to the number of positions at which it can start
there, $T_{g_i} + \bar L - 1$. The model's probability of the read then also includes the chance
of its particular start position, $1 / (T_{g_i} + \bar L - 1)$, which $\Pr(r \mid a)$ leaves out.
Between the haplotypes of one genotype, the two cancel, and every haplotype gets the same weight.
What remains is a factor of $\kappa / \mu_G$ for each read, the same for every haplotype of $G$ but
not for every genotype; in the model it combines with the count of reads. The relative
likelihoods leave the start position out, so vg sets the weights by another rule: which reads can
tell the alleles apart.

A read's factor depends on the weights only when the read fits the alleles of $G$ differently: if
$p_{r g_1} = p_{r g_2}$, the weighted average is the same whatever the weights, because they sum
to 1. Call a read *informative* for $G$ if it fits one of $G$'s alleles better than another. vg
gives each haplotype a weight equal to its expected share of the informative reads. For a true
heterozygote, whose informative reads split between its alleles in that proportion, these weights
approximately maximise the expected read term.

Where the two alleles visit different nodes, a read can favour $g_i$ only where it overlaps nodes
that $g_i$ visits and the other allele does not. For an allele with no such nodes, it can favour
$g_i$ only where it spans the junction at which the other allele's extra sequence would be.
Treating those nodes as one stretch of length $U_i(G)$ (0 for a junction), such a read can start at
$U_i(G) + \bar L - 1$ positions, so each haplotype's weight is proportional to that number:

$$
w_i(G) = \frac{U_i(G) + \bar L - 1}{\sum_{j=1}^{P} \left(U_j(G) + \bar L - 1\right)}
$$

Two alleles can also visit the same nodes: an inversion visits them in the other orientation, and
two alleles of a repeat can go round a loop a different number of times. Then $U_i(G) = 0$ for
both, and the weights are equal. For an inversion that is right, by symmetry. For alleles that
differ in how many times they go round a loop, equal weights are an approximation.

Take a heterozygote of the reference allele and a 1,000-base insertion. The only reads that favour
the reference allele are those spanning its junction, which start at about $\bar L - 1$ positions,
while reads that favour the insertion can start at about $1{,}000 + \bar L - 1$. Equal weights
would read this imbalance as evidence against the heterozygote. Making more of one allele's
sequence differ from the other's, without changing its length, does not make its haplotype produce
more reads; it makes more of them informative, and the weight follows.

A homozygous genotype has $U_i(G) = 0$ for both haplotypes, which therefore get equal weights, as
do two alleles whose unique sequence is equally long, such as two alleles that differ by one base.
At ploidy 1, $w_1(G) = 1$. `--flat-mixture` uses $w_i(G) = 1/P$ instead, so that the effect of the
weighting can be measured. It also flattens the per-allele weights that read phasing uses (see
[From the reads](#from-the-reads)).

#### Depth term inputs

- $\kappa$ is measured over a *rate window* on the reference path. The reference path is cut into
  buckets of a fixed length. A site belongs to the bucket that holds the reference position of its
  start boundary node. If that node has none, vg uses the end boundary node, and failing that, the
  nearest enclosing site. The site's rate window is its bucket and one bucket on each side.
- The window's *read rate* is the number of reads whose alignment begins on a reference node in
  the window, each counted as $1 - e_r$, divided by the reference length of the window. Only
  reference nodes count, in both the number and the length, so the rate does not depend on how the
  graph's nodes are numbered or on how many non-reference nodes the window holds. Variation that
  the sample carries in the window still changes it: a deletion leaves reference bases on which no
  reads start.
- The counts are computed once per bucket and shared, and each site divides the rate by the
  ploidy of its region (set by `-d`, `-R` or `--ploidy-bed`) to give $\kappa$, a rate per haplotype.
  That is the site's own ploidy, except at a nested site that only some of its parent's alleles
  cross: the site is genotyped at a lower ploidy, but the window's reads come from every haplotype.
  When no read begins in the window, $\kappa = 0$, and the site has no depth term.
- A site with no reference position anywhere among its enclosing sites, as in a graph without
  reference path positions, uses a fixed block of consecutive node IDs in place of the rate
  window: the block that contains the site's lowest node ID.
- $\bar L$ is the mean length of the reads that begin in the rate window, or of the site's own
  reads if none does. The mixture weights use the same $\bar L$.
- $N_{\mathrm{eff}}$ counts each read of the site as $1 - e_r$. With `--depth-count-raw`, each read
  counts as 1, in both $N_{\mathrm{eff}}$ and $\kappa$, and $N_{\mathrm{eff}}$ is a whole number.
- $\beta$ is `--depth-term`; 0 turns the depth term off.

The ratio $N_{\mathrm{eff}} / \mu_G$ at the direct call is reported in the VCF as `DR`, whether or
not the depth term is on. A site with $\kappa = 0$ has no `DR`.

## Genotyping

Genotyping chooses each site's genotype from its likelihoods, among the genotypes that can be made
from the site's candidate alleles at the site's ploidy. The genotype chosen is the site's *settled*
genotype. From it vg writes the site's *records*, its lines in the VCF. A site usually has one
record. It can have one for each place where it differs from the reference allele (see
[Reporting each difference once](#reporting-each-difference-once)), or none, as when it is called
homozygous for the reference allele and `-a` is not given.

### Candidate alleles

A site's candidate alleles come from the panel or from the reads.

- When the graph is a GBZ with at least two panel haplotypes, the candidate alleles are by default
  the distinct walks that panel haplotypes take through the site (*haplotype enumeration*). `-z`
  asks for this explicitly, and `-g` gives the panel as a separate GBWT file. Haplotype
  enumeration only offers alleles that some panel haplotype takes. A panel haplotype that passes
  through the site more than once offers each of its walks through it.
- Otherwise, or with `--enumerate-support`, the candidates are the walks with the most read
  support (*support enumeration*). vg finds them with Yen's k-shortest-paths algorithm, over the
  node and edge coverage in a file made by
  [`vg pack`](https://github.com/vgteam/vg/wiki/vg-manpage#pack) and given with `-k`. vg keeps at
  most a fixed number of them.
- The reference allele is always a candidate.

`--max-snarl-edges` sets a limit on a site's edges, counting the edges of the sites nested in it.
vg skips a site over the limit, and genotypes the sites of its child chains as if they were
*top-level* sites, those nested in no other site. They take the ploidy of their contig or region,
and join the contig's top-level sites in the linkage model. A skipped site with no child chains
gets no call. The limit exists because Yen's search is slow on very large sites. It is lifted by
default under haplotype enumeration, which does not use Yen's search.

### Ploidy

A site's ploidy comes from `-d`. `-R` (`--ploidy-regex`) overrides it per contig, with a
comma-separated list of `REGEX:PLOIDY` rules. For each reference path, the first rule whose
regular expression matches the path's name applies. `--ploidy-bed` overrides both per region, with
a BED file of `CHROM START END PLOIDY` lines. `CHROM` is spelled as in the output VCF. Intervals
are 0-based and half-open, and must not overlap. A site takes the ploidy of the interval that
contains its first base on the reference path, the first base of whichever of its boundary nodes
comes first there. A site outside every interval keeps the ploidy from `-d` or `-R`. Every ploidy
must be 1 or 2. A site that these options make haploid has a single allele as its `GT`.

A nested site takes its ploidy from its parent's genotype instead (see
[Which child chains are genotyped](#which-child-chains-are-genotyped)). The exception is a child
of a site that vg could not genotype as a top-level site: the site had no reads,
`--max-snarl-edges` skipped it, or it failed in another way, such as having no reference path
through both of its boundary nodes. vg genotypes such a child as a top-level site, at the ploidy
of its contig or region.

### Direct call or linkage

Every genotype that can be made from the site's candidate alleles, that is every multiset of $P$ of
them, is scored. The genotype with the highest $\mathcal{L}(G)$ is the site's *direct call*. An
exact tie goes to the homozygous reference genotype if it is one of the tied genotypes. Otherwise
it goes to the tied genotype that comes first in VCF genotype order (the order of `GL`), with the
candidate alleles numbered in the order in which they were found, not as the record numbers them.

- Under *direct calling*, that is under support enumeration or with the linkage model off
  (`--linkage-weight 0`), the direct call is the settled genotype.
- Otherwise, under haplotype enumeration, the likelihoods become the evidence of the linkage
  model, described next, and the settled genotype is the one with the highest posterior
  probability under it. A site that the model leaves out (see
  [States and emissions](#states-and-emissions)) keeps its direct call.

A site with no reads gets no genotype and has no place in the linkage model. It gets no record,
unless `-a` is given; its record then has a missing `GT` and `FILTER` `noreads`. If such a site is
top-level, or is genotyped as one, the sites of its child chains are genotyped as top-level sites
(see [Ploidy](#ploidy)). If its ploidy came from its parent's genotype, its child chains are not
genotyped, and nothing nested in it is called.

### Linkage model

Alleles at nearby sites are inherited together, so the sample's alleles at neighbouring sites tend
to form combinations that panel haplotypes also carry. The linkage model uses this by treating each
of the sample's haplotypes as a mosaic of panel haplotypes (the Li–Stephens model), following
PanGenie (Ebler et al. 2022).

The model is a hidden Markov model that runs along a *linkage chain*: a sequence of sites in
reference order. At the top level, a linkage chain holds the top-level sites of one contig, split
wherever the ploidy changes. Below the top level, a linkage chain holds the sites of one child
chain that have the same ploidy, whether or not they are adjacent; at ploidy 1, also on the same
strand of the parent. Where it can (see
[Forward–backward windows](#forwardbackward-windows)), the model starts such a chain with all its
probability on the panel haplotypes that the sample's haplotypes copy at the parent site, as the
parent was phased from the panel (see [From the panel](#from-the-panel)).

The model's hidden state at a site says which panel haplotype each of the sample's haplotypes
copies there. Between sites, a haplotype can switch to copying another panel haplotype, and does so
more readily between distant sites. From the site likelihoods of the genotypes that the states
imply, the model computes each genotype's posterior probability at each site.

#### States and emissions

We call each of the sample's haplotypes a *strand*, numbered 0 and 1. The word names a haplotype of
the sample, not a strand of DNA. At ploidy 2 the hidden state at a site is an ordered pair
$(h_0, h_1)$ of panel haplotypes: strand 0 copies $h_0$ there, and strand 1 copies $h_1$. At
ploidy 1 the state is a single panel haplotype.

A panel haplotype carries at most one allele at a site. It can take the walks of more than one
candidate allele there, when it passes through the site twice or is stored as several paths. vg
then takes it to carry the one that comes last in the list of candidate alleles. That list is in
the order in which the search of the panel's GBWT index finds the alleles, with the reference
allele added at the end if the search did not find it.

A state implies a genotype, made of the alleles that its panel haplotypes carry at the site. The
state's *emission*, the evidence the reads give for it, is that genotype's $\mathcal{L}(G)$.

At each site, the model works over the site's *compact allele set*: the alleles of the site's
direct call, and every allele that some panel haplotype carries there. The model can settle a site
only on a genotype made of these alleles. A site whose compact allele set is larger than a fixed
limit is left out of the model. It keeps its direct call, written unphased, and the model links
the sites on either side of it directly.

#### Wildcard haplotype

The model adds a *wildcard* haplotype to the panel, so that it can call a genotype that no pair of
panel haplotypes carries. The wildcard can carry any allele of the compact allele set.

A strand has an *unknown allele* at a site when it copies the wildcard, or a panel haplotype that
does not pass through the site. The emission of a state with one unknown strand is the mean of
$\mathcal{L}(G)$ over the alleles of the compact allele set, each taken in turn as the unknown
strand's allele. With two unknown strands, the mean is over ordered pairs of these alleles. The
mean is of the likelihoods themselves, and it is then multiplied by a fixed *escape* probability
for each unknown strand.

#### Transitions

Between consecutive sites $d$ bases apart, each strand independently keeps copying the same
haplotype or, with probability $\rho(d)$, draws one uniformly from the panel and the wildcard,
possibly the same one:

$$
\rho(d) = \left[ \rho_{\min} + (1 - \rho_{\min}) \left(1 - e^{-d/D}\right) \right]^{\omega}
$$

A site's position is the first base of whichever of its boundary nodes comes first on the
reference path, and $d$ is the difference between two sites' positions, so distances run from
start to start. $D$ is `--linkage-scale`, the distance over which linkage decays. $\omega$ is
`--linkage-weight`, an exponent on the switch probability: a larger $\omega$ makes switches rarer
and linkage stronger. $\omega = 0$ is a special value: vg does not run the model, and every site
keeps its direct call. $\rho_{\min}$ is a small fixed floor, so that a switch is never impossible.

A chain that no reference path passes through is an *off-reference chain*, and its sites have no
reference position of their own. vg takes the first allele of the parent's settled genotype that
crosses such a site, and places the site at the parent's position plus the site's offset along
that allele. Two sites of one off-reference chain are then as far apart as they are along that
allele. No distance is known between such a site and a site with a reference position, and there
$\rho = 1$.

#### Genotype posterior

The forward–backward algorithm gives each state's posterior probability at each site. A
genotype's posterior collects the probability of the states that imply it:

- A state with no unknown strand implies one genotype.
- A state with one unknown strand, the other carrying allele $a$, shares its probability among the
  genotypes $\lbrace a, b \rbrace$, for each allele $b$ of the compact allele set, in proportion to
  their likelihoods.
- A state with two unknown strands shares its probability among all genotypes, in proportion to
  their likelihoods.

At ploidy 1 a state is a single haplotype, and a state with an unknown strand shares its
probability among the alleles in the same way.

#### Allele-frequency prior

A genotype that many pairs of panel haplotypes carry is implied by many states, so it collects more
probability. This acts as a prior from the panel's allele frequencies. The exponent $F$,
`--linkage-prior`, sets the prior's strength:

- The probability from states of the first kind is multiplied by $c_G^{F-1}$, where $c_G$ is the
  number of such states that imply $G$.
- The probability from states of the second kind is multiplied by $n_a^{F-1}$, where $n_a$ is the
  number of panel haplotypes, the wildcard aside, that carry $a$.
- States of the third kind are not rescaled.

At ploidy 1, the probability of the states that carry allele $a$ is multiplied by $n_a^{F-1}$. The
posteriors are then normalised to sum to 1. $F = 1$ keeps the prior that the states imply, $F = 0$
removes it, and $F > 1$ strengthens it. The model settles the site on the genotype with the highest
posterior.

#### Homopolymer sites

Sequencing errors in a long homopolymer run tend to recur in many reads at the same site, and the
read term counts each such read as independent evidence. A larger exponent $F$ at a site like this
gives the panel's allele frequencies more weight against these reads.

`--hp-prior`, when not 0, replaces $F$ at *homopolymer sites*. At a homopolymer site, the reference
allele and another candidate allele differ only in the length of one homopolymer run, by 1 to 49
copies of its base. The run must also be long: in the longer of the two alleles it has at least
`--hp-prior-run` bases, or it reaches an end of the allele. A run that reaches an end of the allele
can continue into the neighbouring site, so its full length is not known, and it counts as long.

#### Forward–backward windows

vg runs the forward–backward algorithm over overlapping windows of a linkage chain. Each window
*keeps* a fixed number of consecutive sites, and the kept sites of successive windows cover the
chain without overlapping. A window also decodes a fixed margin of sites on each side of its kept
sites, fewer at the ends of the chain, and discards their posteriors. The windows are decoded one
after another, each on its own. Every kept posterior therefore has the margin's number of sites on
each side, except near the ends of the linkage chain, and approximates the posterior from decoding
the whole chain at once. The first window of a chain below the top level starts from the parent's
panel haplotypes where the parent was phased and, for a ploidy-1 chain, the panel haplotype on its
strand passes through the chain's first site. Otherwise, and for every later window, the start is
a uniform distribution. A linkage chain at ploidy 2 that starts from its parent's panel haplotypes
is decoded as one window, whatever its length.

### Nested sites

A site can contain child chains, and their sites can contain chains in turn. *Nested calling*
genotypes the sites of each child chain and writes them in records of their own. Genotyping the
child chains of a site is called *descent*. A site's *generation* counts the descents that reach
it. A site that vg genotypes as a top-level site is generation 0, including the child of a site
that could not be genotyped (see [Ploidy](#ploidy)). A site reached by descent is one generation
after its parent.

Nested calling is on by default with `--read-likelihood`, and `--nested` turns it on with the other
genotyping methods of `vg call`. `--no-nested` turns it off. Each site is then genotyped against
its full walks, and variation inside nested sites is reported in the enclosing site's alleles. A
nested site is genotyped on its own only when its parent could not be genotyped.

#### Reporting each difference once

Two alleles of a site that take the same route except inside a child chain differ only in that
chain. To report such a difference once, vg compares each called allele with the reference allele
in *symbolic* form: the walk with each crossing of a child chain replaced by one symbol for that
chain. A called allele whose symbolic form equals the reference allele's is written as the
reference allele at this site, and the records of the child chain's sites report the difference.

##### Block records

A called allele can also differ from the reference allele in several places, separated by nodes
that both share. `--atomize-blocks`, on by default with nested calling, then writes one record per
difference; `--no-atomize-blocks` turns it off. vg aligns the symbolic form of each called allele
to that of the reference allele, minimising edit distance. Each maximal stretch of nodes and chain
symbols that the alignment does not match is a *block*, and becomes a record. A block record's
`GT` gives the allele that each strand carries over that block.

Blocks of the two strands that overlap or touch on the reference become one record, so that two
records never describe the same stretch of reference differently. If the blocks make only one
record, with as many alleles as one record for the whole site, or if the site cannot be split into
blocks, vg writes one record for the whole site.

A site's block records share the site's ID (the VCF `ID` column) and its evidence, because the
likelihood is computed for the whole site. Their `GQ`, `GQI`, `GQN`, `GP`, `QUAL`, `DP`, `DR` and
`BL` are the site's values (see VCF fields, under Output), except that where the linkage model
moved the site, each block record's `GQN` and `lowconf` come from its own `GL` (see
[Linkage and re-genotyping](#linkage-and-re-genotyping)). Each block allele stands for one site
allele: the block's REF for the site's reference allele, and each ALT for the site allele of the
first strand that carries the ALT. The block record copies its `AD` and `GL` entries from the site
record, reading each block allele as the site allele it stands for. Summing or averaging these
fields over a site's records therefore counts the site's evidence more than once. `INFO/SB` gives
each block record's index, counting from 0, and the number of block records the site writes.

When every crossing of a child chain by the alleles of the parent's settled genotype lies inside a
block, the block's ALT spells out the chain. The records of the chain, and of the sites nested in
it, would repeat that ALT, so they are not written. The child chain is still genotyped and phased
from the panel, because the sites nested in it depend on its genotype and phase.

#### Which child chains are genotyped

A nested site has a copy on a strand exactly when the strand's allele at the parent crosses it. So
a nested site's ploidy is the number of the parent's settled alleles that cross it, and an allele
that crosses it twice counts once. Ploidy is decided site by site, so two sites of one child chain
can have different ploidies.

- **0**: the sample has no copy of the site, and no record is written for it or for anything
  nested inside it.
- **1**: the site is genotyped at ploidy 1. When the genotypes are phased and the site is in a
  [nested haploid chain](#phasing), its `GT` is written `a|.` or `.|a`, and the position of the
  allele says which strand carries the site. Otherwise the `GT` is a single allele.
- **2**: the site is genotyped at ploidy 2.

vg skips a child chain that its parent's reference allele does not cross, except in these cases:

- With `--anchors-out`, unless `--no-off-ref-nesting` is given. The chain's sites then have
  entries in the anchor file. Phasing [from the reads](#from-the-reads) includes these sites, so
  genotyping these chains can change the phase written for other records and, under
  [`--regenotype`](#re-genotyping-from-the-phase), their genotypes.
- When the reference paths include a gRef fragment. `vg paths --compute-gref` adds a *gRef cover*
  (graph-reference cover) to a graph: a copy of each reference path, under a name starting
  `gref_`, and *gRef fragments*. A gRef fragment is a path named `gref_<reference>_<N>_alt` that
  runs through sequence the reference does not cover. A chain on such a gRef fragment gets
  records, with the gRef fragment as their contig. The panel leaves out the gRef fragments, and
  leaves out each gRef copy of a reference path unless the original's sample is absent from the
  GBWT, so that the reference is in the panel once.
- With the environment variable `VG_CALL_NO_REF_NESTED` set, for testing.

An off-reference chain gets no VCF records in any of these cases.

#### Genotyping in two passes

Under the linkage model, a site's genotype is settled only after the model has run over the site's
whole linkage chain. A child site's ploidy depends on its parent's settled genotype. So vg
genotypes in two passes, and writes the records after both.

- In the *sweep*, vg goes through the reads once and computes every site's likelihoods and direct
  call. A child site gets a *provisional ploidy*: the number of the parent's direct-call alleles
  that cross it. Under the linkage model, a child that none of them crosses is genotyped too, at
  the parent's ploidy, since the model may move the parent onto an allele that crosses it; under
  direct calling it is not genotyped. When the child has more than one candidate
  allele, vg also computes its likelihoods and direct call at the other ploidy, because the
  barrier may choose either; a child with one candidate allele keeps its provisional ploidy unless
  it is dropped. Each site is *staged*: vg keeps what the site's records will be built from, and
  writes nothing.
- In the *barrier*, vg settles and phases the sites one generation at a time, parents before
  children. The linkage model settles generation 0, and vg phases it from the panel (see
  [From the panel](#from-the-panel)). Each generation-1 site then takes its ploidy from its
  parent's settled genotype, and vg uses the likelihoods and direct call that the sweep computed
  at that ploidy. A site at ploidy 0 is dropped with everything nested in it. The model then
  settles generation 1, each linkage chain starting from the panel haplotypes that the parent's
  strands copy. Each later generation follows in the same way.

After the last barrier pass, vg *renders* each staged site: it builds and writes the site's records
from its settled genotype and phase.

The barrier can run more than once (see [Rounds](#rounds)). Each time, it decides every site's
ploidy afresh, so a site dropped once can come back the next time. vg records which of a parent's
candidate alleles cross a child site only when the parent has at most a fixed number of
candidates. A child of a parent with more keeps its provisional ploidy. Under direct calling there
is no linkage model, so every site keeps its direct call and its provisional ploidy.

## Phasing

Phasing decides each genotype's [phase](#sample-and-reads): which of the sample's
[strands](#states-and-emissions) carries which allele. At ploidy 2, a site's phase is an order of
the two alleles of its settled genotype. The first allele is on strand 0 and is written to the left
of the `|` in `GT`, and the second is on strand 1.

A diploid site is *phaseable* when its settled genotype holds two different candidate alleles, so
that its two possible phases differ. A phaseable site can still be homozygous in the VCF, when its
two alleles differ only inside a child chain.

A *nested haploid chain* is a linkage chain of ploidy-1 sites whose parent is diploid, or is a site
of another nested haploid chain. All its sites lie on one strand of the nearest diploid ancestor,
the strand that carries them, directly or through ploidy-1 parents (see
[Which child chains are genotyped](#which-child-chains-are-genotyped)). Its sites' `GT` is `a|.`
on strand 0 and `.|a` on strand 1.

### From the panel

vg phases each linkage chain from its *Viterbi path*: the most probable sequence of the linkage
model's states along the chain, among those that imply the settled genotype at every site. A state
in which one strand has an unknown allele qualifies when the other strand carries one of the
settled genotype's alleles, and a state in which both strands do always qualifies. So phasing
orders a genotype without changing it. For the path, each qualifying state's emission is the
settled genotype's $\mathcal{L}(G)$, times the escape probability for each unknown strand. A site
settled in an earlier generation, such as a parent decoded with its child chain, is held where it
can be to the state it was phased with. On every linkage chain, the path is decoded over
overlapping windows of the sizes given in [Forward–backward windows](#forwardbackward-windows).
Each window after the first is held, at the previous window's last kept site, to the state that
the previous window chose there.

Each record's `FORMAT/PS` names its *phase set*: the reference position of the first site of its
top-level linkage chain. A nested site takes its parent's phase set. Sites of one phase set are
phased relative to one another.

At some phaseable sites, neither strand of the Viterbi path copies a panel haplotype that carries
either called allele. Nothing then orders the pair. vg still writes it with `|`, in an arbitrary
order, and the site stays in its phase set. Read phasing, when on, can still order it.

Phasing is on wherever the linkage model runs. `--phased` makes vg call fail when the linkage model
does not run (see [Direct call or linkage](#direct-call-or-linkage)), and `--no-phased` turns
phasing off. Where the linkage model runs, `--no-phased` also turns nested calling off (an
explicit `--nested` is then an error), and variation inside nested sites is reported in the
enclosing site's alleles. Nested calling needs phasing there because the barrier takes each nested
site's strand, and the panel haplotypes its linkage chain starts from, from its parent's phase.

### From the reads

A read that spans two phaseable sites shows directly whether their alleles lie on the same strand.
`--read-phasing` uses such reads to re-decide the phase of each phaseable site, nested and
off-reference sites included, within the phase sets the panel gave. It leaves out the chains that a
[block record](#block-records) spells out. Read phasing never changes a genotype. It is off by
default and on under `--preset ont`.

Read phasing takes one phase set at a time, and its phaseable sites in order of position (see
[Transitions](#transitions)). A read is identified by its name, so the two mates of a pair count as
one read when they reach different sites. Where both mates reach the same site, the site has two
entries under one name. Each counts as a separate read in the site's reliability (defined below).
A link pairs the two sites' entries of one name one to one, so a mate that reaches only one of the
sites adds nothing. In coherence, a mate's vote at a site includes its mate's entry there, so a
pair counts even when both mates reach only that site.

#### What each read says

Take a read $r$ of a phaseable site $s$, and let $a_0$ and $a_1$ be the site's alleles on strand 0
and strand 1 in its current phase. For $k \in \lbrace 0, 1 \rbrace$, let
$x_{rsk} = (1 - e_r) v_{a_k} p_{r a_k}$, where $e_r$ is the read's
[mismapping probability](#mismapping-probability) and $p_{r a_k}$ its
[relative likelihood](#likelihood-formula) under $a_k$ at $s$. The *allele-length weights*
$v_{a_0}, v_{a_1}$ are the [mixture weights](#mixture-weights) with $\lvert a_k \rvert$ in place
of $U_i(G)$. $\lvert a \rvert$ is the length of all of allele $a$'s node visits, boundary nodes
included, while $U_i(G)$ counts only the nodes that one allele visits and the other does not. A
read with $x_{rs0} + x_{rs1} = 0$ is not used at $s$. For each read used, we keep:

- $q_{rs} = x_{rs0} / (x_{rs0} + x_{rs1})$, the probability that the read carries $a_0$, given
  that it came from one of the two strands;
- $c_{rs} = (x_{rs0} + x_{rs1}) / (x_{rs0} + x_{rs1} + e_r)$, the probability that it did come from
  one of them, rather than being mismapped;
- its *confidence*,
  $-10 \log_{10}\left(1 - \max(x_{rs0}, x_{rs1}) / (x_{rs0} + x_{rs1} + e_r)\right)$, the
  phred-scaled probability that its better allele is wrong.

#### Links between sites

Two sites $s$ and $t$ are *linked* by the reads they share. For a shared read $r$, let

$$
m_r = q_{rs} q_{rt} + (1 - q_{rs})(1 - q_{rt})
$$

This is the probability that the read's alleles at the two sites lie on one strand, given the
sites' current phases. The read is assumed to report this relation truly with probability
$\gamma_r = c_{rs} c_{rt}$, and to be a coin flip otherwise. The link is

$$
\mathrm{link}(s, t) = \sum_{r} \log_{10} \frac{\gamma_r m_r + (1 - \gamma_r)/2}{\gamma_r (1 - m_r) + (1 - \gamma_r)/2}
$$

A positive link favours the two sites' current phases. A link's *size* is its absolute value.
`--phase-cap`, when not 0, limits the size of each link.

#### Reliable sites

A site's *reliability* is the mean confidence of its reads, and the site is *reliable* if this is
at least `--phase-min-q`. Consider a read with $e_r = \epsilon_{\min}$, the floor on $e_r$, at a
site whose two alleles have equal length. If the read fits one allele perfectly and the other not
at all, its confidence is
$-10 \log_{10}\left(\epsilon_{\min} / (\epsilon_{\min} + (1 - \epsilon_{\min})/2)\right)$, the
*heterozygous score ceiling*. A read that favours the longer of two unequal alleles can have a
higher confidence. vg rejects a `--phase-min-q` above the ceiling, since sites whose alleles are of
similar length could not reach it.

#### Deciding the phases

To *flip* a site is to reverse its phase. For each site, read phasing decides whether to flip it
against the phase the panel gave. It phases the reliable sites first, each relative to the one
before it, and then phases every other site on its own. A wrong link between reliable sites flips
every site after it, while a wrong decision for another site flips only that site. Each phase set
goes through four stages:

1. **Phase chain.** The reliable sites of the phase set, in order of position, form its *phase
   chain*, and each is linked to the next. The phase chain breaks wherever the size of a link is
   below `--phase-break` $\log_{10}$ units, and the breaks divide it into *pieces*. The first site
   of each piece keeps the panel's phase. Each later site of the piece is flipped relative to the
   site before it when their link is negative.
2. **Relink.** At each break, read phasing decides from the links across it whether to flip the
   later piece. Take the last `--phase-relink` sites of the piece before the break and the first
   `--phase-relink` sites of the piece after it, or all of a piece's sites if it has fewer. Sum the
   links between every site on one side and every site on the other, each in the current phases.
   If the sum is negative, every site of the later piece is flipped. If no read links the two
   sides, the later piece keeps the phases stage 1 gave it. Breaks are decided from left to right,
   so each piece is oriented against the piece before it as that piece now stands.
3. **Coherence.** This stage removes phase-chain sites whose reads disagree with the rest of the
   phase chain. Take a read that spans two or more phase-chain sites, and one of those sites, $s$.
   The read's other phase-chain sites $t$ vote for the strand it came from. Strand 0 scores
   $\sum_t \log_{10}(c_{rt} q_{rt} + (1 - c_{rt})/2)$, strand 1 scores
   $\sum_t \log_{10}(c_{rt}(1 - q_{rt}) + (1 - c_{rt})/2)$, and the higher score wins. The read
   *agrees* at $s$ if $q_{rs}$ lies on the winning strand's side of $1/2$. Every $q$ here is taken
   in the current phase. A site's *coherence* is the fraction of its reads that agree, counting
   only reads that span another phase-chain site. A site with at least a fixed number of counted
   reads and a coherence below `--phase-coherence` is removed from the phase chain. If any site is
   removed, stages 1 and 2 run again on the sites left, starting again from the panel's phases,
   and then this stage runs again. This stops when the stage removes no site, or after it has
   removed sites `--phase-coh-rounds` times (once if that is 0). The stage runs only on a phase
   chain of at least 3 sites, and removes nothing if fewer than 2 would be left.
   `--phase-coherence 0` skips this stage.
4. **Hang.** Each phaseable site of the phase set that is not in the phase chain, because it is
   unreliable or stage 3 removed it, is then phased on its own. Its nearest phase-chain sites are
   used, up to $\lfloor H/2 \rfloor + 1$ on each side, where $H$ is `--phase-hang`. The link to each
   of them votes for the phase it implies, with a weight equal to its size. A further vote of
   weight `--phase-prior`, in the same $\log_{10}$ units, favours flipping the site exactly when
   the nearest phase-chain site is flipped. The site takes the phase with the larger total vote,
   and a tie flips it.

When read phasing flips a site, its alleles change strands, and so does every nested haploid chain
that takes its strand from them, directly or through another nested haploid chain. A diploid
nested site keeps its own phase, which read phasing decides like any other site's.

### Re-genotyping from the phase

In the site likelihood, every read has the same [mixture weights](#mixture-weights), because
nothing at one site shows which strand a read came from. Once sites are phased, a read's alleles at
the other phaseable sites it spans do show it. `--regenotype` uses this to give each read its own
weights, tilted towards the strand its other sites place it on, and so corrects each site's
likelihoods. Off-reference chains, and chains that a block record spells out, keep their sweep
likelihoods. It needs `--read-phasing` and, like it, is off by default and on under
`--preset ont`. The quantities $e_r$, $p_{ra}$, $q_{rt}$, $c_{rt}$ and the allele-length weights
$v$ are those of [From the reads](#from-the-reads). $\sigma$ is the logistic function, and
$\mathrm{logit}$ is its inverse.

#### Strand log-odds

A read's *strand log-odds* at site $s$, $\Lambda_{rs}$, measures how strongly its other sites
place it on strand 0. Each other phaseable site $t$ at which $r$ is used adds one term: the
natural-log odds that the read came from strand 0, with $q_{rt}$ taken in the phase read phasing
chose. A site adds one term per read name, even where both mates of a pair are used there.

$$
\Lambda_{rs} = \sum_{t \neq s} \ln \frac{c_{rt} q_{rt} + (1 - c_{rt})/2}{c_{rt}(1 - q_{rt}) + (1 - c_{rt})/2}
$$

Leaving out $s$ keeps a site from confirming its own genotype. A read used in more than one phase
set has no usable strand, because each phase set labels its strands independently, and its
$\Lambda_{rs}$ is 0.

#### Tempering

The terms of $\Lambda_{rs}$ are not independent evidence. A wrong phase at some of the read's
sites, or an error that the read makes the same way at several sites, enters several terms at once.
So $\Lambda_{rs}$ overstates how sure the strand is. The *tempered* strand log-odds is

$$
y_{rs} = \mathrm{logit}\left(C \sigma(\tau \Lambda_{rs}) + (1 - C)/2\right)
$$

where $\tau$ is the *temper* (`--regeno-temper`) and $C$ is `--regeno-ceiling`. With $C = 1$,
$y_{rs} = \tau \Lambda_{rs}$. A smaller $C$ keeps the probability of either strand between
$(1 - C)/2$ and $(1 + C)/2$.

Unless `--regeno-temper` is given, $\tau$ is fitted once, from the phase that read phasing gave
before the first correction. Each read at each phaseable site $s$ gives an *observation* when
$\Lambda_{rs} \neq 0$ and $q_{rs} \neq 1/2$. Each of the two points to a strand: strand 0 when
$\Lambda_{rs} > 0$ or $q_{rs} > 1/2$, and strand 1 otherwise. The observation records whether they
point to the same strand. The observations are sorted by $\vert \Lambda_{rs} \vert$ and grouped
into bins. Every bin but the last holds the larger of a fixed number and a fixed fraction of the
observations, and the last holds the remainder. For each bin, the tempered probability of the
strand, $C \sigma(\tau \vert \Lambda \vert) + (1 - C)/2$ at the bin's mean $\vert \Lambda \vert$,
predicts its rate of agreement. $\tau$ is chosen from a fixed grid, with $C$ held at its given
value, to minimise the squared difference between the predicted and observed rates, with each bin
weighted by its size. With too few observations, $\tau$ is 0 and re-genotyping changes nothing.

#### Likelihood correction

Take a genotype of two different alleles, $a$ on strand 0 and $b$ on strand 1, whose allele-length
weights are $v_a$ and $v_b$. Each read gets its own weights

$$
\pi_{ra} = \frac{v_a e^{y_{rs}}}{v_a e^{y_{rs}} + v_b}, \qquad \pi_{rb} = 1 - \pi_{ra}
$$

The correction added to $\ln \mathcal{L}(G)$ sums over the site's reads $R$:

$$
\sum_{r \in R} \left[ \ln\left((1 - e_r)(\pi_{ra} p_{ra} + \pi_{rb} p_{rb}) + e_r\right) - \ln\left((1 - e_r)(v_a p_{ra} + v_b p_{rb}) + e_r\right) \right]
$$

The other assignment of the two alleles to the strands is scored too, and the larger correction is
kept. Homozygous genotypes are unchanged and have no assignment to choose, so this choice can only
favour heterozygous genotypes.

Both terms of the correction use the allele-length weights $v$, while the site likelihood's read
term used the mixture weights $w$. The two agree when the alleles have the same length and neither
visits a node twice.

#### Nested haploid chains

At a site of a nested haploid chain, each genotype is a single allele on one strand, so there are no
weights to tilt. Instead, a read is made less informative when the phase places it on the other
strand. Let $y$ be the read's tempered strand log-odds towards the chain's strand, and
$\eta = \min(1, e^{y})$. Each of the read's relative likelihoods $p$ becomes $\eta p + 1 - \eta$.
$\eta$ is 1 when $\tau = 0$ or when the read points to the chain's strand, so those reads are
unchanged. `--no-regeno-haploid` turns this off.

#### Rounds

After a correction, the barrier runs again on the corrected likelihoods, and read phasing runs again
on the genotypes it settles. A correction, the barrier and read phasing together make a *round*. The
new genotypes can change the phase, and with it $\Lambda$, so another round can follow. Each round
corrects the likelihoods from the sweep, not those of the previous round. `--regeno-passes` caps the
number of times the barrier runs, the first run included. The rounds stop sooner when the
correction changes no site's direct call, when the barrier changes no settled genotype, or when the
settled genotypes return to an earlier state. Every round's correction is settled by the barrier,
including the round that stops, so the genotypes are settled from the likelihoods that `GL` reports.
With `--regeno-passes 1` the correction is computed and reported but not applied. `--regeno-ledger` writes one line for each site whose direct call the
last correction changed.

## Output

### VCF fields

| Field | Meaning |
|---|---|
| `GT` | the settled genotype, phased where phasing ran |
| `GL` | $\log_{10} \mathcal{L}(G)$ for every genotype of the record's alleles, in VCF order |
| `GQ` | the difference between the log-likelihoods of the direct call and the *runner-up*, the genotype with the second-highest $\mathcal{L}(G)$, in phred units. It is multiplied by the *explained share*, the fraction of reads whose best allele is in the direct call (a tied read split as for `AD`), and by the `--depth-quality` factor where that applies |
| `GQI` | the same difference, with neither factor |
| `GQN` | the same difference divided by the *achievable gap* (below), held at 1 or less, and multiplied by the explained share |
| `GP` | one value: the natural log of the posterior probability of the direct call, computed from $\mathcal{L}(G)$ with a uniform prior over genotypes. (In the VCF specification, `GP` is a phred-scaled value per genotype.) |
| `QUAL` | phred-scaled posterior probability, under the same uniform prior, of the genotype whose alleles are all the reference allele; 0 when `GT` is all reference |
| `DP` | number of reads of the site |
| `AD` | for each allele in the record, the number of reads whose best allele it is, rounded; a read tied between alleles counts a fraction to each |
| `DR` | $N_{\mathrm{eff}} / \mu_G$ at the direct call: the effective read count over the count that the direct call predicts |
| `BL` | mean over reads of $\max_a \ell_{ra}$, each read's best log-likelihood score at the site (see [Relative likelihood](#relative-likelihood)) |
| `FORMAT/PS` | the phase set |
| `INFO/SB` | on a block record, its index counting from 0 among the site's block records, and the number of block records the site writes (see [Reporting each difference once](#reporting-each-difference-once)) |
| `FILTER=noreads` | the site had no reads, so no genotype is called (`GT` is `./.`). Such a record is written only with `-a`, which also writes reference calls and sites with no reads |
| `FILTER=lowconf` | `GQN` is below `--min-confidence` |

The per-site `GQ` is a difference of log-likelihoods, where the VCF specification's `GQ` is a
posterior. On records that the linkage model moved (below), `GQ` comes from a posterior instead, so
the `GQ` values of the two kinds of record are not comparable. A read whose best allele is in
neither the direct call nor the runner-up fits both about equally, and adds almost nothing to their
difference. The explained share lowers `GQ` where the direct call leaves such reads unexplained.
`GL`, `GQ`, `GQI`, `GP` and `QUAL` are over-confident at high depth, because reads are treated as
independent. `AD` need not sum to `DP`: every candidate allele was scored, but only the alleles
written in the record have an entry.

#### Normalised quality

`GQ` depends on depth, since the likelihood difference is a sum over reads. It also depends on
ploidy. At ploidy 1 the runner-up is a different allele, and most reads can tell it apart from the
direct call. At ploidy 2 the runner-up usually differs from the direct call on one strand only, so
fewer reads tell the two apart, or each read tells them apart less.

`GQN` removes both effects by dividing by the *achievable gap*. This is the difference that the
read term alone would give between the direct call and the runner-up if each of the site's reads
were ideal. Ideal reads are shared among the direct call's haplotypes in proportion to its
[mixture weights](#mixture-weights). Each ideal read has $e_r = \epsilon_{\min}$ (`--mismap-min`),
and fits its own haplotype's allele with [relative likelihood](#relative-likelihood) 1 and every
other allele with 0.

The numerator of `GQN` is the observed difference in $\ln \mathcal{L}$, depth term included. It can
exceed the achievable gap, so the fraction is held at 1. `GQN` is `.` when there is no difference to
normalise (no reads, or a single possible genotype). Otherwise it lies in $[0, 1]$, except on moved
records.

#### Linkage and re-genotyping

A record is *moved* when the linkage model settles a genotype other than the direct call. Its
`GQ` and `GQN` are computed again for the settled genotype, from the direct call's explained share,
`GQ` factor and achievable gap, which the linkage model keeps for each site. The direct call is the
one the site entered the linkage model with: its call in the sweep, or, for a nested chain the
barrier records at the ploidy its settled parent implies, its call at that ploidy.

A moved record's `GQ` is $-10 \log_{10}(1 - \text{posterior})$ times the direct call's `GQ`
factor, then capped at `GQI`. The posterior is the linkage model's posterior probability of the
settled genotype, so it includes the panel's prior. The `GQ` factor is what the per-site `GQ`
multiplies its difference by: the explained share unless `--no-share-quality`, times the
`--depth-quality` factor where that applies. The cap keeps `GQ` at or below the reads' own
confidence in the direct call.

A moved record's `GQN` is the margin, in phred units, of the settled genotype's entry in `GL` over
the largest other entry, divided by the direct call's achievable gap and multiplied by its
explained share, as the per-site `GQN` is, and held within $[-1, 1]$. It is negative where `GL`
favours another genotype over the settled one, as it does on a moved whole-site record whose direct
call is among the record's genotypes. It is `.` where the direct call had no achievable gap.
`lowconf` is decided again from this `GQN`, and cleared where it is `.`. Under re-genotyping, a
record is moved if the last barrier run moved it (see [Rounds](#rounds)).

Re-genotyping is applied with `--regenotype` when `--regeno-passes` is above 1 (see
[Rounds](#rounds)). `GL` and `QUAL` are then written from the corrected likelihoods. In a round
whose correction changes the genotype with the highest $\mathcal{L}(G)$, `GQ` is recomputed from
them as the per-site `GQ` is, for that genotype, but with the sweep's `--depth-quality` factor.
Each round starts again from the sweep's likelihoods and `GQ`, so `GQ` follows the last round's
correction, and is the sweep's where that correction left the best genotype alone. `GP`, `GQI`,
`GQN`, `DR` and `lowconf` keep their values from before the correction, so they describe the direct
call made from the uncorrected likelihoods. A moved record takes the `GQ`, `GQN` and `lowconf` described
above, whether or not re-genotyping ran.

#### Options that change the fields

`-L` merges similar called ALT alleles, as `vg deconstruct -L` does, at sites at least
`--cluster-min-len` long. Two alleles are merged when their length-weighted similarity is at least
the value of `-L`. `GT`, `AD` and `GL` are rewritten for the merged alleles, and `INFO/MAT` records
the merge. It is off by default.

`--no-share-quality` leaves the explained share out of the per-site `GQ`; `GQN` keeps it.
`--depth-quality` multiplies the per-site `GQ` by $e^{-D_q \vert \ln \mathrm{DR} \vert}$, where
$D_q$ is its value, at records where an allele of the direct call differs in length from the
reference allele by at least a fixed number of bases. `--min-confidence` marks records with
`FILTER=lowconf` and does not remove them. These three options change only quality fields and
`FILTER`, never a genotype.

#### Nesting tags

`vg call` has three other ways of genotyping nested sites, chosen by general options: `-A`
genotypes every site on its own, `--top-down` genotypes children after their parents and takes a
child's candidate alleles from its parent's genotype, and `--bottom-up` genotypes children before
their parents. With any of them, and whenever [off-reference chains](#transitions) are genotyped,
records carry vg's nesting INFO tags:

- `INFO/LV` counts the record's enclosing sites that have records on the record's own contig.
- `INFO/CH` counts the changes of contig on the way up from the record through the records of its
  enclosing sites, or gives the level of the record's contig if that is larger. A reference path
  that is not a [gRef fragment](#which-child-chains-are-genotyped) is level 0, a gRef fragment
  whose ends attach to a level-0 path is level 1, one whose ends attach to a level-1 fragment is
  level 2, and so on. A site inside a deletion stays on its parent's contig. When the reference
  paths include gRef fragments, a site inside an insertion is reported on one, so its `INFO/CH` is
  at least 1.
- `INFO/PS` is the ID of the record of the nearest enclosing site that has one. It is unrelated to
  `FORMAT/PS`.
- `INFO/RC`, `INFO/RS` and `INFO/RD` give a contig, start and end at which to look the record up,
  normally those of its outermost enclosing site; the VCF header gives the details.

### Mosaic (`--mosaic-out`)

The mosaic file describes each of the sample's strands (its haplotypes) as a walk through the
graph, and says which panel haplotype the strand copies along each part of the walk. The file is the
phasing in another form, so `--mosaic-out` is rejected without the linkage model or with
`--no-phased`.

The file is tab-separated. Header lines start with `#` and are identified by their first field:

| Key | Meaning |
|---|---|
| `#mosaic-version` | the format version |
| `#graph` | the input graph |
| `#sample` | the sample name, set with `-s` |
| `#reference` | a reference path that positions refer to, one line for each, except gRef fragments |
| `#gref-fragments` | the number of gRef fragments among the reference paths, when there are any |
| `#decoding` | how the strands were chosen: `constrained-viterbi`, the Viterbi path restricted to the settled genotypes, as in [From the panel](#from-the-panel) |
| `#patch`, `#nested`, `#unexplained` | the choices of the three options described below |
| `#haplotype` | a panel haplotype's index and name, one line for each |
| `#note` | text describing the columns |
| `#H` | the column names |

A reader should skip header keys it does not recognise. Where read phasing reversed a site's order,
the two strands' panel haplotypes are swapped there too, so each strand keeps the haplotype that
carries its allele.

#### Segments and rows

A data line, or *row*, starts with `H`. Rows are built from *segments*. A segment is a maximal
stretch of consecutive sites on one strand, over which the strand copies one panel haplotype.
Consecutive segments of a strand join end to end where the walk can continue from one to the next.
Each maximal walk so formed is a *mosaic fragment*, and a strand can have several. A new mosaic
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
expands to one walk in the graph, counting each shared node once. A row never crosses the end of
one of the paths that store its panel haplotype, so a segment over several such paths gives several
rows. `gbwt_offset` is valid only for the graph named in `#graph`, and `hap_index` only within one
file; `haplotype` is the name to compare across files. A contig whose phased sites are all haploid
has rows for strand 0 only. In a haploid region of a diploid contig, strand 1 is on the wildcard.

#### Forming segments

A strand's sites include the nested sites on its walk, so a change of panel haplotype at a nested
site starts a new segment. Only sites with a record count, whether written as one record or as
[block records](#block-records); a site with no record is left out.

Between two consecutive segments, the walk follows the first segment's panel haplotype if that
haplotype continues to the second segment. Otherwise it follows the second segment's haplotype, if
that haplotype reaches back to the first. Where neither does, the stretch between them is a *gap in
the walk*. The gap is filled with the reference, on a `ref` row, if the reference is a panel
haplotype that crosses it. The linkage model can have a strand copy a panel haplotype at a site
that the haplotype does not pass through, so a segment's haplotype need not cover the whole
segment. Such a segment is replaced by a `ref` row where the reference crosses it.

Three options change how the rows are formed, and the header records each choice:

- `--no-mosaic-patch-gaps` turns off both uses of the reference, so gaps stay unfilled, and a new
  mosaic fragment starts after each.
- `--no-mosaic-nested` leaves nested sites out of the segments, so a strand follows its enclosing
  site's haplotype through them.
- Where a strand is on the [wildcard](#wildcard-haplotype), the panel cannot name a haplotype for
  it. By default those sites are left out of the segments, and the walk crosses them by the rule
  above. `--mosaic-break-unexplained` writes a `*` row for them instead.

### Assembly anchors (`--anchors-out`)

The anchor file records, for each genotyped site, which reads support which of the sample's
strands, for use in pangenome-guided assembly.

A *pin* is a point between two adjacent bases of the graph, with no sequence of its own. Each site
has two. The *start pin* lies just after the site's start boundary node, where an allele's walk
enters the interior. The *end pin* lies just before its end boundary node, where the walk leaves.

A site's reads are divided among its *slots*. A slot's number names a field of `GT`: slot 0 the
first allele, to the left of the `|`, and slot 1 the second. A phaseable site has one slot for each
allele. Any other diploid site has a single slot, 0, unless `--anchors-hom-split` divides its reads
between slots 0 and 1, which then carry the same allele. A haploid site also has a single slot, 0,
except at a site of a nested haploid chain whose `GT` is `.|a`, where it is 1. An *anchor* is one
pin together with the reads of one slot that cross it.

#### Placing reads in slots

For each distinct called allele $a$, let $x_{ra} = (1 - e_r) v_a p_{ra}$, where $v_a$ are the
allele-length weights of [From the reads](#from-the-reads), and $v_a = 1$ at a site that is not
phaseable. Each read goes to the slot whose allele has the largest $x_{ra}$, except where the read
phase decides (below). A read's *anchor confidence* is
$-10 \log_{10}\left(1 - x_{ra} / (\sum_b x_{rb} + e_r)\right)$, where $a$ is the allele of the slot
it is placed in. The sum runs over distinct alleles, so the two slots of a split site count their
allele once. A site's *anchor reliability* is the mean anchor confidence of the reads written for
it, each read counted once, after the read filters below.

#### Using the read phase

With `--read-phasing`, read placement also uses each read's tempered strand log-odds $y_{rs}$ (see
[Re-genotyping from the phase](#re-genotyping-from-the-phase)), computed leaving the site out. The
temper $\tau$ is the one re-genotyping used, where it ran and $\tau$ was above 0. Otherwise it is
fitted in the same way from the final phase.

- At a phaseable site, a read's slot is by default chosen from its $x_{ra}$ with $v_a$ replaced by
  the per-read weights $\pi$ of [Likelihood correction](#likelihood-correction). The read's anchor
  confidence is still computed from the $x_{ra}$. `--no-anchors-phase-hets` chooses the slot from
  the $x_{ra}$ alone. `--anchors-strict-hets` chooses it from the sign of $y_{rs}$ alone, slot 0 for
  a positive sign, and a read whose $y_{rs}$ is 0 then keeps the slot its $x_{ra}$ give.
- `--anchors-hom-split` divides the reads of a diploid site that is not phaseable between two slots
  by the sign of their $y_{rs}$. A site is split only when, for each sign, at least
  `--split-min-side` reads have a $y_{rs}$ of that sign and of absolute value at least
  `--split-min-q`. At a split site, a read whose $y_{rs}$ is 0, such as a read seen in more than one
  phase set, is placed by a coin flip derived from its name, so it takes the same slot at every
  site.

#### Rows of the anchor file

The file is tab-separated. Each anchor is written as an `A` row, followed by one `R` row for each
of its reads. Header lines start with `#`, and a reader should skip keys it does not recognise.
Among them, the `#read` lines list the read names, which `R` rows refer to by index. The `#H` lines
name the columns of each kind of row, and the `#note` lines describe them.

| `A` column | Meaning |
|---|---|
| `node` | the ID of the pin's boundary node: the first node of `snarl` for the start pin, the second for the end pin |
| `snarl` | the site's ID, as in the VCF `ID` column: its start and end boundary node IDs, each preceded by `>` or `<` for its orientation. It is also written for off-reference sites, which have no VCF record |
| `slot` | the slot |
| `allele` | the slot's allele, as its index in the site's list of candidate alleles rather than a VCF allele number. The list is not written, so the index serves to compare a site's slots |
| `gqn` | the site's `GQN`, computed as for the VCF; `.` where it has none |
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
- `--anchors-min-gqn`, when above 0, skips sites whose `GQN` is below its value or missing. At 0
  or below it skips nothing.
- `--anchors-min-q` skips reads whose anchor confidence is below its value.
- `--anchors-keep-off-call` keeps reads whose best candidate allele was not called; by default they
  are left out.
- `--anchors-reads` skips anchors with fewer reads than its value.
- With `--anchors-end-new` set to $N$, a slot's end anchor is written only if at least $N$ of the
  slot's reads cross the end pin but not the start pin.

## Options

`vg call --help` gives each option's default. `--preset ont` sets several of the options below to
values suited to Oxford Nanopore reads, and `vg call --help` lists which. An option given explicitly
overrides the preset. Where an option has a `--no-` form, the two set the same thing and the one
given last wins, so a `--no-` form also turns off a setting that a preset turned on.

The table includes general `vg call` options: `-k`, `-g`, `-z`, `--max-snarl-edges`, `--nested`,
`--no-nested`, `--atomize-blocks`, `--no-atomize-blocks`, the Ploidy row, `-L` and
`--cluster-min-len`. Its other options are rejected without `--read-likelihood`. The options that
modify `--anchors-out` (`--no-off-ref-nesting` among them), `--mosaic-out` and `--regenotype` are
also rejected when those are not in use. Other general options used on this page, such as `-a`,
`-s`, `-r`, `-p`, `-P`, `-S`, `-A`, `--top-down` and `--bottom-up`, are described by
`vg call --help`.

| Part | Options |
|---|---|
| Reads | `--gam`, `--gaf-reads`, `--gam-index`, `--gaf-base`, `--gbz-base`, `--gaf-base-binary`, `--read-window`, `--read-min-mapq` |
| Candidate alleles | `--enumerate-support`, `-k`, `-g`, `-z`, `--max-snarl-edges` |
| Relative likelihood | `--gap-open`, `--gap-extend`, `--insertion-nats`, `--realign`, `--no-realign` |
| Mismapping | `--mismap-min`, `--mismap-max`, `--no-mismap-term` |
| Mixture weights | `--flat-mixture` |
| Depth term | `--depth-term`, `--depth-count-raw` |
| Linkage | `--linkage-weight`, `--linkage-scale`, `--linkage-prior`, `--hp-prior`, `--hp-prior-run` |
| Nested sites | `--nested`, `--no-nested`, `--atomize-blocks`, `--no-atomize-blocks`, `--no-off-ref-nesting` |
| Ploidy | `-d`, `-R`, `--ploidy-bed` |
| Phasing | `--phased`, `--no-phased`, `--read-phasing`, `--no-read-phasing`, `--phase-min-q`, `--phase-coherence`, `--phase-coh-rounds`, `--phase-break`, `--phase-relink`, `--phase-hang`, `--phase-prior`, `--phase-cap` |
| Re-genotyping | `--regenotype`, `--no-regenotype`, `--regeno-temper`, `--regeno-ceiling`, `--regeno-passes`, `--regeno-haploid`, `--no-regeno-haploid`, `--regeno-ledger` |
| Quality fields | `--no-share-quality`, `--depth-quality`, `--min-confidence` |
| Allele merging | `-L`, `--cluster-min-len` |
| Mosaic | `--mosaic-out`, `--mosaic-patch-gaps`, `--no-mosaic-patch-gaps`, `--no-mosaic-nested`, `--mosaic-break-unexplained` |
| Anchors | `--anchors-out`, `--anchors-reads`, `--anchors-min-gqn`, `--anchors-min-q`, `--anchors-het-only`, `--anchors-leaf-only`, `--anchors-keep-off-call`, `--anchors-end-new`, `--anchors-phase-hets`, `--no-anchors-phase-hets`, `--anchors-strict-hets`, `--anchors-hom-split`, `--split-min-q`, `--split-min-side` |
| Debugging and evaluation | `--dump-likelihoods` (writes each read's $e_r$ and $p_{ra}$ at every site to a file, as TSV), `--regeno-shuffle` (randomises the sign of each $\Lambda_{rs}$ before re-genotyping, as a control for how much the phase contributes); `--flat-mixture` and `--anchors-strict-hets`, listed above, also serve to measure the parts they replace |
| Presets | `--preset` |

### Fixed constants

These constants have no option. Each is located by its name, or by the function or member that
holds it.

| Constant | Defined in | What it sets |
|---|---|---|
| `LinkageModel::Params::escape` | `src/linkage_model.hpp` | the escape probability for a strand with an unknown allele |
| `LinkageModel::Params::rho_min` | `src/linkage_model.hpp` | the floor $\rho_{\min}$ on the switch probability |
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

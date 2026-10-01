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
   with the highest likelihood, the [direct call](#direct-call-or-linkage), or use a **linkage
   model**. The linkage model also uses the haplotypes stored in the graph, which tend to carry the
   same combinations of alleles at neighbouring sites as the sample. With
   [nested calling](#nested-sites), on by default, sites nested inside other sites are genotyped
   too, and get VCF records of their own.
3. **Phasing.** vg assigns each genotype's alleles to the sample's haplotypes, first from the
   stored haplotypes and then, optionally, from reads that span several sites. With
   `--regenotype`, the phase is then used to correct the likelihoods, and genotyping and phasing
   are repeated.
4. **Output.** vg writes a VCF file. It can also write a **mosaic** file, which describes each of
   the sample's haplotypes as a walk through the graph, and an
   [anchor file](#assembly-anchors---anchors-out), which ties reads to haplotypes.

The linkage model decides many neighbouring sites together, and a nested site's genotype depends
on its parent's. So when the linkage model or nested calling is on, as both are by default, vg
first computes the likelihoods of every site. It then genotypes and phases the sites one
[generation](#nested-sites) at a time, parents before children.

## Vocabulary

### Graph

- **Node.** A node holds a DNA sequence and can be traversed in either **orientation**: forward,
  reading its sequence, or reverse, reading its reverse complement. A **walk** is a sequence of
  oriented node visits along the graph's edges.
- **Site.** A place in the graph where the sample's genome may differ from the reference. vg
  identifies sites as **snarls**. A snarl is a subgraph separated from the rest of the graph by two
  **boundary nodes**, a start and an end. On this page a site is a snarl that `vg call` genotypes.
  Its boundary nodes belong to it, and its other nodes, including those of any snarls nested in
  it, are its **interior**. A site is oriented from its start boundary node to its end boundary
  node, and its alleles are written in that direction.
- **Chain.** A series of snarls joined end to end, each snarl's end boundary node being the next
  one's start. Snarls nest: a snarl can contain chains of smaller snarls, its **child chains**. The
  sites in a site's child chains are **nested** in it, and it is their **parent**.
- **Allele.** A walk through a site from its start boundary node to its end boundary node, also
  called a **traversal**.
- **Reference.** The **reference paths** are the graph paths on which `vg call` reports positions.
  By default they are the paths that the graph's metadata marks as reference or generic paths;
  `-p`, `-P` or `-S` selects other paths instead, by name, by name prefix or by sample. A site's
  **reference allele** is the walk a reference path takes through it. In the VCF the reference
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
  as $\lbrace 0, 1 \rbrace$. Which haplotype carries which allele is the genotype's **phase**, and
  is decided separately.
- **Read placements.** A mapper, such as `vg giraffe`, aligns the reads to the graph. The reads'
  **placements** are the walks their alignments take, together with the **edits** inside each node
  (the runs of matching bases, mismatches, insertions and deletions against the node's sequence)
  and their mapping qualities (MAPQ).
- **Reads of a site.** The reads whose placements visit an interior node of the site, or visit
  both of its boundary nodes. A read that visits both boundary nodes and no interior node supports
  an allele that skips the interior, such as a deletion. A read whose only node in the site is one
  boundary node fits every allele equally, and is not used.

## Site likelihood computation

The site likelihood $\mathcal{L}(G)$ measures how well a genotype $G$ explains the reads of one
site. We build it from a model of how a site's reads arise from the sample's genotype. We write
down the model's probability of the reads, and then approximate each of its parts by a quantity
that can be computed from the read placements.

### Model of how reads arise

We treat the reads of a site as generated from the sample's genotype $G$ in four steps:

1. **Draw the reads' mapping qualities.** The site receives some number of reads, each with a
   MAPQ, and so with a probability $e_r$ of being mismapped (see the next step). Only their
   **effective read count** $N_{\mathrm{eff}}$, the sum of $1 - e_r$ over the reads, depends on
   $G$. It is drawn from a distribution centred near the **expected read count** $\mu_G$, which
   grows with the length of $G$'s alleles (see [Depth term](#depth-term)). Every set of MAPQs
   with that effective read count is equally likely.
2. **Decide whether each read is mismapped.** Independently for each read, with probability
   $e_r$, the read is **mismapped**: the mapper put it in the wrong place. Either it came from
   elsewhere in the genome, or it came from this site but its placement gives it the wrong start
   or the wrong local alignment. Otherwise the read is correctly mapped.
3. **Choose the haplotype that each correctly mapped read came from.** A correctly mapped read
   came from $G$'s $i$-th haplotype, which carries allele $g_i$, with probability $w_i(G)$, the
   haplotype's **mixture weight**.
4. **Copy the read from that haplotype's allele.** The read's bases inside the site are a copy of
   part of the allele, starting where the read's placement puts it, with sequencing errors. The
   error model is the one behind vg's alignment scores: base substitutions have probabilities set
   by base quality, and insertions and deletions have fixed costs (see
   [Relative likelihood](#relative-likelihood)). A mismapped read's bases come from wherever the
   read really came from.

Each step leaves a different kind of evidence about $G$ in the reads. The number of reads reflects
the alleles' lengths: a homozygous deletion shows mainly as reads that are missing. The chance of
mismapping means that a read with a doubtful placement counts for less. Most of the evidence comes
from the copying step, because a read fits the allele it was copied from better than it fits the
others.

The model makes the following assumptions. Where the data break one, expect wrong answers.

- **Reads are independent given the genotype.** This does not mean that the reads are unrelated:
  the reads from one haplotype share its allele, and the genotype accounts for that. The
  assumption fails for paired mates, and for reads that share a systematic error. Their evidence is
  then counted more than once, and the likelihoods grow over-confident with depth.
- **Whether a read is mismapped depends only on its MAPQ.** This does not mean that MAPQ is taken
  at face value: $e_r$ is held between a floor and a ceiling (see
  [Mismapping probability](#mismapping-probability)). It does mean that a read that fits every
  allele badly is no more likely to be mismapped than one that fits well. The assumption fails for
  a read from another copy of a repeat that the mapper placed here with a high MAPQ, which then
  counts as evidence against the alleles it fits worst.
- **The effective read count depends only on the genotype.** This does not mean that depth is
  uniform along the genome, since the read rate in $\mu_G$ is measured near the site, or that the
  count follows a Poisson distribution exactly, since the distribution is widened (see
  [Depth term](#depth-term)). The assumption fails where mappability, base composition or a change
  in copy number moves the depth at one site away from that of its neighbourhood. The depth term
  then favours genotypes whose allele lengths fit the wrong depth.
- **Each haplotype passes through the site once.** This does not mean that an allele visits each
  node once: a walk that loops inside the site, between its boundary nodes, is one allele. The
  assumption fails at a duplication, where a haplotype's path passes through the site twice. The
  model still gives that haplotype one allele, and explains the reads of both copies with it, so
  the site has more reads than $G$ predicts.

### Notation

| Symbol | Meaning |
|---|---|
| $A$ | the site's candidate alleles: the alleles vg considers at the site (see [Candidate alleles](#candidate-alleles)) |
| $a, b$ | alleles in $A$ |
| $P$ | the ploidy at the site |
| $G = \lbrace g_1, \dots, g_P \rbrace$ | a genotype; $g_i$ is the allele of its $i$-th haplotype. The haplotypes are listed in any order, and nothing below depends on that order. A homozygous genotype lists the same allele $P$ times |
| $\mathcal{L}(G)$ | the site likelihood of $G$ |
| $R$ | the reads of the site |
| $r$ | a single read in $R$ |
| $\Pr(r \mid G)$ | the probability of read $r$'s bases under genotype $G$ |
| $\Pr(r \mid a)$ | the probability of read $r$'s bases if it was correctly mapped and copied from allele $a$ |
| $\Pr(r \mid \text{mismapped})$ | the probability of read $r$'s bases if it was mismapped |
| $p_{ra}$ | the relative likelihood of read $r$ under allele $a$ |
| $e_r$ | the mismapping probability of read $r$ |
| $\epsilon_{\min}, \epsilon_{\max}$ | the floor and ceiling on $e_r$ |
| $w_i(G)$ | the mixture weight of $G$'s $i$-th haplotype |
| $N_{\mathrm{eff}}$ | the effective read count: the expected number of the site's reads that were not mismapped, given their MAPQs |
| $\mu_G$ | the expected read count: the number of correctly mapped reads that $G$ predicts |
| $f_{\mathrm{Pois}}$ | the Poisson probability, extended to counts that are not whole numbers |
| $h_\beta$ | the distribution of the effective read count |
| $\beta$ | the parameter that widens $h_\beta$, `--depth-term` |
| $T_a$ | the length in bases of allele $a$, excluding the site's two boundary nodes; a node the allele visits twice counts twice |
| $U_i(G)$ | at ploidy 2, the length in bases of the nodes that $g_i$ visits and the other allele of $G$ does not, each node counted once |
| $\bar L$ | the mean read length near the site |
| $\kappa$ | the read-start rate: the expected number of correctly mapped reads that begin at each base of one haplotype near the site |

### Likelihood formula

Under the model, the probability of the site's reads given $G$ is

$$
\Pr(R \mid G) = c_R \, h_\beta\left(N_{\mathrm{eff}} ; \mu_G\right) \prod_{r \in R} \Pr(r \mid G)
$$

The first two factors come from step 1. $h_\beta$ is the density of the effective read count, and
$c_R$, the probability of the reads' particular MAPQs among the sets with the same effective read
count, does not depend on $G$. The product over reads comes from steps 2 to 4, which happen
independently for each read.

We compute the site likelihood as

$$
\mathcal{L}(G) = \prod_{r \in R} \left[ (1 - e_r) \sum_{i=1}^{P} w_i(G) p_{r g_i} + e_r \right] \times \frac{f_{\mathrm{Pois}}\left(N_{\mathrm{eff}} ; \mu_G\right)^{\beta}}{\hat Z_\beta(\mu_G)}
$$

The product over reads, the **read term**, represents $\prod_{r \in R} \Pr(r \mid G)$. The last
factor, the **depth term**, represents $h_\beta(N_{\mathrm{eff}} ; \mu_G)$, with
$\hat Z_\beta$ an approximation of its normaliser $Z_\beta$ (see [Depth term](#depth-term)). Each is approximately
proportional to its counterpart in the model, with a factor that does not depend on $G$, so

$$
\mathcal{L}(G) \approx \frac{\Pr(R \mid G)}{C_R}
$$

for some $C_R$ that depends on the reads but not on $G$. Such a factor changes neither which
genotype is most likely nor the ratio of two genotypes' likelihoods. The sections below derive each
part, and state each approximation where it is made. vg works with $\ln \mathcal{L}(G)$, the sum of
the logarithms of the factors.

#### Read term

Take one read $r$ of the site. By steps 2 to 4 of the model, it arose in one of $P + 1$ mutually
exclusive ways: it was correctly mapped and came from $G$'s $i$-th haplotype, for one $i$ from 1
to $P$, or it was mismapped. Its probability is the sum over these cases:

$$
\Pr(r \mid G) = \sum_{i=1}^{P} (1 - e_r) \, w_i(G) \Pr(r \mid g_i) + e_r \Pr(r \mid \text{mismapped})
$$

$\Pr(r \mid a)$ is the probability of the read's bases inside the site if the read was copied from
allele $a$, starting where its placement puts it. We compute it, up to a factor that depends only
on the read, from an alignment of the read to the allele (see
[Relative likelihood](#relative-likelihood)). [Mixture weights](#mixture-weights) and
[Mismapping probability](#mismapping-probability) give $w_i(G)$ and $e_r$.

The model does not say how probable a mismapped read's bases are, because it does not model the
rest of the genome. We approximate that probability by the read's probability under its best
allele at the site:

$$
\Pr(r \mid \text{mismapped}) \approx \max_{b \in A} \Pr(r \mid b)
$$

The mapper placed the read here because it resembles this site, so the sequence that the read
really came from probably explains it about as well as the site's best allele does. This choice
has a second use: it lets us divide every read's probability by the same value.

The likelihood multiplies the probabilities of all the site's reads together. If we divide one
read's probability by a value that does not depend on $G$, we divide every genotype's likelihood by
that value, and $C_R$ absorbs it. We divide each read's probability by
$\max_{b \in A} \Pr(r \mid b)$, which leaves

$$
\frac{\Pr(r \mid G)}{\max_{b \in A} \Pr(r \mid b)} \approx (1 - e_r) \sum_{i=1}^{P} w_i(G) \, p_{r g_i} + e_r
$$

where $p_{ra}$ is the read's **relative likelihood** under allele $a$:

$$
p_{ra} = \frac{\Pr(r \mid a)}{\max_{b \in A} \Pr(r \mid b)}
$$

This is the read's factor in $\mathcal{L}(G)$. A relative likelihood is 1 for the read's best
allele, and smaller for alleles that explain the read less well. The mixture weights are
probabilities that sum to 1, so $\sum_{i} w_i(G) p_{r g_i}$ is the expected value of the read's
relative likelihood over which of $G$'s haplotypes it came from. The read's factor therefore lies
between $e_r$ and 1.

Dividing by the best allele's probability does four things. We compute $\Pr(r \mid a)$ only up
to a factor that depends on the read, and that factor cancels in $p_{ra}$, so we never need it. It
lets reads scored with and without base qualities, by scorers with different scales, share one
matrix of relative likelihoods. Each read's factor lies between $e_r$ and 1, so its logarithm is
finite and does not underflow, however long the read. And a read's relative likelihoods can be read
directly: 1 for its best allele, and near 0 for an allele that fits it far worse.
`--dump-likelihoods` writes them.

#### Mixture weights

The mixture weight $w_i(G)$ is the prior probability, before we look at the read's bases, that a
correctly mapped read of the site came from $G$'s $i$-th haplotype. The weights of a genotype sum
to 1. We do not compute this prior exactly. Instead we approximate it where it matters.

It matters only for some reads. If a read fits $G$'s alleles equally well, so that
$p_{r g_1} = p_{r g_2}$, its expected relative likelihood is the same whatever the weights, because
they sum to 1. Call a read **informative** for $G$ if it fits one of $G$'s alleles better than
another. Only the informative reads depend on the weights, so we approximate the prior for an
informative read: each haplotype's expected share of the informative reads.

We approximate that share in two steps. First, we take a read to be informative for $g_i$ only
where it overlaps nodes that $g_i$ visits and the other allele of $G$ does not. Where $g_i$ has no
such nodes, as the reference allele has none against an insertion, a read is informative for
$g_i$ only where it spans the junction at which the other allele's extra sequence would be.
Second, we treat $g_i$'s unique nodes as one stretch, of length $U_i(G)$, which is 0 at a
junction. An informative read can then start at $U_i(G) + \bar L - 1$ positions on the haplotype,
and the weights are proportional to that number:

$$
w_i(G) = \frac{U_i(G) + \bar L - 1}{\sum_{j=1}^{P} \left(U_j(G) + \bar L - 1\right)}
$$

A homozygous genotype has $U_i(G) = 0$ for both haplotypes, which therefore get equal weights, as
do two alleles whose unique sequence is equally long, such as two alleles that differ by one base.
At ploidy 1, $w_1(G) = 1$.

This is a rough approximation, and it could be improved. An allele's unique nodes are seldom one
stretch. $U_i(G)$ counts each node once, in either orientation, so two alleles that visit the same
nodes get equal weights. For an inversion, which visits them in the other orientation, that is
right by symmetry. For two alleles of a repeat that go round a loop a different number of times,
it is not: when the loop is longer than a read, the allele with more turns of the loop has
informative reads that the other lacks.

`--flat-mixture` sets $w_i(G) = 1/P$ instead, so that the effect of the weighting can be measured.
It also flattens the per-allele weights that read phasing uses (see
[From the reads](#from-the-reads)).

#### Mismapping probability

Read $r$ is mismapped if it came from elsewhere in the genome, or if it came from this site but its
placement gives it the wrong start or the wrong local alignment. MAPQ is the mapper's estimate, on
the phred scale, of the probability of the first case. vg's mappers compute it from the scores of
the read's distinct placements in the graph, so it measures whether the read belongs at another
locus, not whether it is aligned correctly through this site. The mismapping probability $e_r$ is
the sum of the probabilities of the two cases. We approximate the second by a constant, the floor
$\epsilon_{\min}$ (`--mismap-min`), and the sum by the larger of its two terms. We also hold $e_r$
below a ceiling $\epsilon_{\max}$ (`--mismap-max`):

$$
e_r = \min\left(\max\left(10^{-\mathrm{MAPQ}_r / 10}, \epsilon_{\min}\right), \epsilon_{\max}\right)
$$

The larger of two probabilities is at least half their sum, so this approximation is at least
half the probability it stands for.

Every read is mismapped with probability at least $\epsilon_{\min}$, so the model can always
explain a read that fits $G$'s alleles badly as mismapped. One read can therefore change the ratio
of two genotypes' likelihoods by at most a factor of $1 / \epsilon_{\min}$, and the higher the
floor, the less any single read's fit counts. (The read also adds $1 - e_r$ to
$N_{\mathrm{eff}}$, which moves the depth term.)

The ceiling applies to reads with MAPQ 0 or close to it. Many mappers give MAPQ 0 to a read with
several equally good placements. Its unclamped $e_r$ would be 1, which would make the read's
factor the same under every genotype, so that the read counted for nothing. The ceiling decides
how much such a read still counts. Such reads are used unless `--read-min-mapq` excludes them.

`--no-mismap-term` sets every $e_r$ to $\epsilon_{\min}$, whatever the read's MAPQ, as if every
read were well mapped, so that the contribution of the mismapping case can be measured. The
effective read count and $\kappa$ are computed from $e_r$ too (see
[Depth term inputs](#depth-term-inputs)), so under this option every read counts as
$1 - \epsilon_{\min}$ in them.

#### Depth term

The depth term represents step 1 of the model: how well the effective read count fits $G$.
$N_{\mathrm{eff}}$ is usually not a whole number, so it cannot follow a Poisson distribution, and
we give it a continuous analogue of one. For a real $n \geq 0$ and $\lambda > 0$, let

$$
f_{\mathrm{Pois}}(n ; \lambda) = \frac{\lambda^{n} e^{-\lambda}}{\Gamma(n + 1)}
$$

where the gamma function $\Gamma$ extends the factorial to real numbers: $\Gamma(n + 1) = n!$ for
a whole number $n$. For a whole number $n$, $f_{\mathrm{Pois}}(n ; \lambda)$ is the Poisson
probability of $n$ events when $\lambda$ are expected. The model draws $N_{\mathrm{eff}}$ from the
distribution on $n \geq 0$ with density

$$
h_\beta(n ; \lambda) = \frac{f_{\mathrm{Pois}}(n ; \lambda)^{\beta}}{Z_\beta(\lambda)}, \qquad Z_\beta(\lambda) = \int_0^\infty f_{\mathrm{Pois}}(x ; \lambda)^{\beta} \, dx
$$

with $\lambda = \mu_G$. At $\beta = 1$ this is a continuous analogue of the Poisson distribution.
The parameter $\beta$, `--depth-term`, widens it: for large $\lambda$ it is close to a normal
distribution with mean $\lambda$ and variance $\lambda / \beta$, where a Poisson distribution has
variance $\lambda$. We widen it because read depth varies between sites more than a Poisson count
would, for reasons the model leaves out. At $\beta = 0$, $f_{\mathrm{Pois}}^{\beta} = 1$ under every
genotype, so the count says nothing about $G$, and the depth term is off.

Evaluated at the observed $N_{\mathrm{eff}}$, $h_\beta(N_{\mathrm{eff}} ; \lambda)$ is the
likelihood of $\lambda$. We compute $Z_\beta(\lambda)$ from the normal approximation: for large
$\lambda$, $f_{\mathrm{Pois}}(n ; \lambda)$ is close to a normal density in $n$ with mean and
variance $\lambda$, so

$$
\ln Z_\beta(\lambda) \approx \frac{1 - \beta}{2} \ln (2 \pi \lambda) - \frac{1}{2} \ln \beta
$$

Against numerical integration at $\beta = 0.1$, $0.5$ and $1$, this is within 0.13 of
$\ln Z_\beta(\lambda)$ for every $\lambda \geq 2$, and within 0.02 for $\lambda \geq 30$.
Below $\lambda = 2$ it falls away from the integral, so there we use its value at $\lambda = 2$.
At $\beta = 1$ the approximation is 0. At $\beta < 1$, $Z_\beta(\lambda)$ grows with $\lambda$,
about as $\lambda^{(1 - \beta)/2}$, so the normaliser counts against the genotype with the larger
expected read count, by about $\frac{1 - \beta}{2} \ln (\mu_G / \mu_{G'})$ in
$\ln \mathcal{L}$ between genotypes $G$ and $G'$. That is 0 between genotypes whose alleles have
equal lengths, such as at a SNP. At $\beta = 0.1$, when one genotype expects twice as many reads as
the other, it is 0.31.

The expected read count is

$$
\mu_G = \kappa \sum_{i=1}^{P} \left(T_{g_i} + \bar L - 1\right)
$$

$T_{g_i} + \bar L - 1$ is the number of positions at which a read of length $\bar L$ can start on
the haplotype and still include some of the interior of $g_i$. When $g_i$ has no interior
($T_{g_i} = 0$), it is the number of positions from which a read reaches across the junction
between the two boundary nodes.

### Computing the terms

The read term needs, for each read of the site, its relative likelihoods and its mismapping
probability. The depth term needs the effective read count, and the read-start rate and mean read
length that set the expected read count. Those two come from the reads that begin near the site,
not only from the site's own reads.

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
  else from the input graph.

`--read-window` sets the size of those ranges, in node IDs, for the two indexed sources. It changes
which reads are fetched together and the order in which they arrive, but not which reads a site
uses, and no result depends on it. A site's reads are put in a fixed order, by read name and then
by where the alignment begins, before anything is summed over them, and the reads that begin in a
rate window (see [Depth term inputs](#depth-term-inputs)) are counted by MAPQ.

A read whose alignment crosses the site against the direction of the alleles is
reverse-complemented before it is scored. We decide the direction by a vote over the read's visits
to the nodes that the alleles visit, leaving out any node that two alleles visit in opposite
orientations. The read is reversed if more of those visits are in the opposite orientation to the
alleles' than in the same one. A tie, including a read that visits none of those nodes, leaves it
as it is. We do not rely on the boundary nodes alone, because a read that lies inside the site's
interior visits neither of them. The boundary nodes do take part in the vote, and they decide it
at an inversion, whose inverted nodes the alleles visit in both orientations.

#### Relative likelihood

We compute a read's relative likelihoods by aligning the read to each candidate allele, scoring each
alignment, and comparing each score with the best one. The alignment is of node visits, not of
bases. The read's placement and the allele are both sequences of node visits, and two visits are the
**same visit** when they go to the same node in the same orientation. We align the read's visits
inside the site, boundary nodes included, to the allele's visits in much the way that dynamic
programming aligns two DNA sequences with affine gap costs, with node visits in place of bases.

Such an alignment is a **pairing**. Each read visit is either paired with one allele visit or left
unpaired, and the pairs come in the same order in both sequences. Two kinds of pair are
considered:

- **A read visit paired with the same allele visit.** It scores the read's own edits in that node,
  as the mapper aligned them.
- **A substitution.** It pairs a read visit that the allele never makes with an allele visit that
  the read never makes. It compares the read's bases in its node with the allele node's sequence,
  base by base from the first base of each, over the shorter of the two lengths, and adds a gap
  for the difference in length.

A visit that the read and the allele share is never paired with a different visit; such pairs are
not considered. Visits left unpaired are gaps, measured in bases:

- A run of consecutive read visits left unpaired is an **insertion**, scored as one gap as long as
  their bases. A visit that the read and the allele share may be left unpaired, which can explain
  a read with poor edits in that node better. Before the read's first pair of same visits, each
  unpaired read visit is a gap of its own.
- A run of allele visits left unpaired between two pairs is a **deletion**, scored as one gap as
  long as their bases. Unpaired allele visits before the read's first pair of same visits, or
  after its last pair, lie outside the read and score nothing, so that a read is not penalised for
  being short.

Every read base inside the site is therefore scored under every allele, and all alleles are scored
over the same read bases.

Scores use vg's quality-adjusted alignment scoring, in which a mismatch at a low-quality base costs
less, with gap scores set by `--gap-open` and `--gap-extend`. A read without base qualities is
scored without the quality adjustment. As in any alignment scoring, a better fit scores higher. A
pairing determines a base-level alignment of the read's bases in the site to the allele's sequence,
and the pairing's score $s_{ra}$ is the score of that alignment, except that each pair and each gap
is scored on its own, so a gap at the edge of one is not joined to a gap in the next.

##### Optimal pairing

With `--optimal-pairing`, we use optimal pairing: the pairing with the highest score $s_{ra}$,
found by dynamic programming over the read's visits against the allele's, as in affine-gap
alignment. Only the pairing of visits is searched; bases are not aligned again. When the
product of the numbers of read visits and allele visits exceeds a fixed limit, the dynamic
programming is banded (see [Fixed constants](#fixed-constants)), so at the largest sites optimal
pairing can miss the best pairing.

##### Greedy pairing

Without `--optimal-pairing`, we use greedy pairing, which approximates optimal pairing in one
pass along the read's visits. It pairs each read visit with the next matching allele visit: the
next occurrence of the same visit after the last pair. Allele visits skipped over between two pairs
are a deletion. When there is no matching visit, it looks for a simple insertion: if the allele's
next unpaired visit is one that the read makes later, the read visit is left unpaired. Otherwise it
pairs the read visit with the allele's next unpaired visit as a substitution, if neither of the two
visits is made elsewhere by the other sequence, and leaves the read visit unpaired if one is. Once
the allele's visits are used up, the remaining read visits are left unpaired. A pair, once made, is
never revised.

Greedy pairing considers the same pairs as optimal pairing and scores the pairing it finds by
the same rules, including one gap for each run of unpaired read visits, so the two differ only in
which pairing they find.

##### From scores to relative likelihoods

We convert the pairing's score to a log-likelihood score in nats:

$$
\ell_{ra} = \alpha s_{ra} + \iota I_{ra}
$$

$\alpha$ is the scorer's log base. vg computes it from the substitution matrix so that $\alpha$
times a substitution score is the natural log of a likelihood ratio: the probability of the aligned
bases under the alignment model, over their probability as unrelated random sequence. This is how
vg interprets alignment scores as probabilities elsewhere. The scorers with and without the
quality adjustment have different log bases. $I_{ra}$ is the number of insertions in the pairing
in which the read has bases the allele lacks: those inside the mapper's edits, each unpaired read
visit, and each substitution whose read node is the longer. Each unpaired read visit counts once
here, even where a run of them is scored as one gap. $\iota$ is
`--insertion-nats`. A positive $\iota$ makes an insertion, where the read has bases the allele
lacks, cost less than a deletion of the same length, where the allele has bases the read lacks.
Optimal pairing chooses the pairing by $s_{ra}$ alone, and then adds $\iota$ for its insertions.

These scores define the error model of the copying step. Given the alignment that the pairing
describes, the probability of the read's bases in the site is $\exp(\ell_{ra})$ times $B_r$, their
probability as unrelated random sequence. $B_r$ is the same for every allele, because every allele
is scored over the same read bases. We take the probability of one alignment, the best pairing
found, in place of the sum over all alignments of the read to the allele, so

$$
\Pr(r \mid a) \approx B_r \exp\left(\ell_{ra}\right)
$$

We score the read against every candidate allele, and then divide by the best. Substituting this
approximation into the definition of $p_{ra}$, $B_r$ cancels:

$$
p_{ra} = \frac{\Pr(r \mid a)}{\max_{b \in A} \Pr(r \mid b)} \approx \frac{B_r \exp\left(\ell_{ra}\right)}{B_r \exp\left(\max_{b \in A} \ell_{rb}\right)} = \exp\left(\ell_{ra} - \max_{b \in A} \ell_{rb}\right)
$$

#### Depth term inputs

The depth term compares the effective read count $N_{\mathrm{eff}}$ with the expected read count
$\mu_G = \kappa \sum_{i} (T_{g_i} + \bar L - 1)$, under a distribution of width set by $\beta$. Its
inputs are computed in this order:

1. **The read-start rate $\kappa$**, the expected number of correctly mapped reads that begin at
   each base of one haplotype near the site. It is measured over a **rate window** on the
   reference path. The reference path is cut into buckets of a fixed length. A site belongs to the
   bucket that holds the reference position of its start boundary node. If that node has none, vg
   uses the end boundary node, and failing that, the nearest enclosing site. The site's rate window
   is its bucket and one bucket on each side. The window's **read rate** is the number of reads
   whose alignment begins on a reference node in the window, each counted as $1 - e_r$, divided by
   the reference length of the window. Only reference nodes count, in both the number and the
   length, so the rate does not depend on how the graph's nodes are numbered or on how many
   non-reference nodes the window holds. Variation that the sample carries in the window still
   changes it: a deletion leaves reference bases on which no reads start. The counts are computed
   once per bucket and shared, and each site divides the read rate by the ploidy of its region (set
   by `-d`, `-R` or `--ploidy-bed`) to give $\kappa$, a rate per haplotype. That is the site's own
   ploidy, except at a nested site that only some of its parent's alleles cross: the site is
   genotyped at a lower ploidy, but the window's reads come from every haplotype. When no read
   begins in the window, $\kappa = 0$, and the site has no depth term. A site with no reference
   position anywhere among its enclosing sites, as in a graph without reference path positions,
   uses a fixed block of consecutive node IDs in place of the rate window: the block that contains
   the site's lowest node ID.
2. **The mean read length $\bar L$**, the mean length of the reads that begin in the rate window,
   or of the site's own reads if none does. The mixture weights use the same $\bar L$.
3. **The effective read count** $N_{\mathrm{eff}} = \sum_{r \in R} (1 - e_r)$, over the site's
   reads.
4. **The width parameter $\beta$**, `--depth-term`.

The expected read count $\mu_G$ is then computed for each genotype from $\kappa$, $\bar L$ and the
lengths of the genotype's alleles.

`--depth-count-raw` counts each read as 1 in place of $1 - e_r$, in both $N_{\mathrm{eff}}$ and
$\kappa$, so that the contribution of the mismapping probabilities to the depth term can be
measured. $N_{\mathrm{eff}}$ is then the number of the site's reads, a whole number, and step 1 of
the model draws that number.

The VCF's `DR` field, the **depth ratio**, is the effective read count divided by the expected
read count, $N_{\mathrm{eff}} / \mu_G$, for the genotype the record is written with: the site's
[direct call](#direct-call-or-linkage), the genotype with the highest $\mathcal{L}(G)$, or the
genotype the linkage model settles. A value near 1 means that the site has as many reads as that
genotype predicts. `DR` is written whether or not the depth term is on, and is left out where
$\kappa = 0$. Where an allele of the written genotype was not among the site's candidate alleles,
or the record's ploidy is not the one the site was genotyped at, it is the direct call's.

## Genotyping

Genotyping chooses each site's genotype from its likelihoods, among the genotypes that can be made
from the site's candidate alleles at the site's ploidy. The genotype chosen is the site's
**settled** genotype. From it vg writes the site's **records**, its lines in the VCF. A site usually
has one record. It can have one for each place where it differs from the reference allele (see
[Reporting each difference once](#reporting-each-difference-once)), or none, as when it is called
homozygous for the reference allele and `-a` is not given.

### Candidate alleles

A site's candidate alleles come from the panel or from the reads.

- When the graph is a GBZ with at least two panel haplotypes, the candidate alleles are by default
  the distinct walks that panel haplotypes take through the site (**haplotype enumeration**). `-z`
  asks for this explicitly, and `-g` gives the panel as a separate GBWT file. Haplotype
  enumeration only offers alleles that some panel haplotype takes. A panel haplotype that passes
  through the site more than once offers each of its walks through it.
- Otherwise, or with `--enumerate-support`, the candidates are the walks with the most read
  support (**support enumeration**). vg finds them with Yen's k-shortest-paths algorithm, over the
  node and edge coverage in a file made by
  [`vg pack`](https://github.com/vgteam/vg/wiki/vg-manpage#pack) and given with `-k`. vg keeps at
  most a fixed number of them.
- The reference allele is always a candidate.

`--max-snarl-edges` sets a limit on a site's edges, counting the edges of the sites nested in it.
vg skips a site over the limit, and genotypes the sites of its child chains as if they were
**top-level** sites, those nested in no other site. They take the ploidy of their contig or region,
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
them, is scored. The genotype with the highest $\mathcal{L}(G)$ is the site's **direct call**. An
exact tie goes to the homozygous reference genotype if it is one of the tied genotypes. Otherwise
it goes to the tied genotype that comes first in VCF genotype order (the order of `GL`), with the
candidate alleles numbered in the order in which they were found, not as the record numbers them.

- Under **direct calling**, that is under support enumeration or with the linkage model off
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

The model is a hidden Markov model that runs along a **linkage chain**: a sequence of sites in
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

We call each of the sample's haplotypes a **strand**, numbered 0 and 1. The word names a haplotype
of the sample, not a strand of DNA. At ploidy 2 the hidden state at a site is an ordered pair
$(h_0, h_1)$ of panel haplotypes: strand 0 copies $h_0$ there, and strand 1 copies $h_1$. At ploidy
1 the state is a single panel haplotype.

A panel haplotype carries at most one allele at a site. It can take the walks of more than one
candidate allele there, when it passes through the site twice or is stored as several paths. vg
then takes it to carry the one that comes last in the list of candidate alleles. That list is in
the order in which the search of the panel's GBWT index finds the alleles, with the reference
allele added at the end if the search did not find it.

A state implies a genotype, made of the alleles that its panel haplotypes carry at the site. The
state's **emission**, the evidence the reads give for it, is that genotype's $\mathcal{L}(G)$.

At each site, the model works over the site's **compact allele set**: the alleles of the site's
direct call, and every allele that some panel haplotype carries there. The model can settle a site
only on a genotype made of these alleles. A site whose compact allele set is larger than a fixed
limit is left out of the model. It keeps its direct call, written unphased, and the model links
the sites on either side of it directly.

#### Wildcard haplotype

The model adds a **wildcard** haplotype to the panel, so that it can call a genotype that no pair of
panel haplotypes carries. The wildcard can carry any allele of the compact allele set.

A strand has an **unknown allele** at a site when it copies the wildcard, or a panel haplotype that
does not pass through the site. The emission of a state with one unknown strand is the mean of
$\mathcal{L}(G)$ over the alleles of the compact allele set, each taken in turn as the unknown
strand's allele. With two unknown strands, the mean is over ordered pairs of these alleles. The
mean is of the likelihoods themselves, and it is then multiplied by a fixed **escape** probability
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

A chain that no reference path passes through is an **off-reference chain**, and its sites have no
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

`--hp-prior`, when not 0, replaces $F$ at **homopolymer sites**. At a homopolymer site, the
reference allele and another candidate allele differ only in the length of one homopolymer run, by 1
to 49 copies of its base. The run must also be long: in the longer of the two alleles it has at
least `--hp-prior-run` bases, or it reaches an end of the allele. A run that reaches an end of the
allele can continue into the neighbouring site, so its full length is not known, and it counts as
long.

#### Forward–backward windows

vg runs the forward–backward algorithm over overlapping windows of a linkage chain. Each window
**keeps** a fixed number of consecutive sites, and the kept sites of successive windows cover the
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

A site can contain child chains, and their sites can contain chains in turn. **Nested calling**
genotypes the sites of each child chain and writes them in records of their own. Genotyping the
child chains of a site is called **descent**. A site's **generation** counts the descents that reach
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
in **symbolic** form: the walk with each crossing of a child chain replaced by one symbol for that
chain. A called allele whose symbolic form equals the reference allele's is written as the
reference allele at this site, and the records of the child chain's sites report the difference.

##### Block records

A called allele can also differ from the reference allele in several places, separated by nodes
that both share. `--atomize-blocks`, on by default with nested calling, then writes one record per
difference; `--no-atomize-blocks` turns it off. vg aligns the symbolic form of each called allele
to that of the reference allele, minimising edit distance. Each maximal stretch of nodes and chain
symbols that the alignment does not match is a **block**, and becomes a record. A block record's
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
- When the reference paths include a gRef fragment. `vg paths --compute-gref` adds a **gRef cover**
  (graph-reference cover) to a graph: a copy of each reference path, under a name starting
  `gref_`, and **gRef fragments**. A gRef fragment is a path named `gref_<reference>_<N>_alt` that
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

- In the **sweep**, vg goes through the reads once and computes every site's likelihoods and direct
  call. A child site gets a **provisional ploidy**: the number of the parent's direct-call alleles
  that cross it. Under the linkage model, a child that none of them crosses is genotyped too, at
  the parent's ploidy, since the model may move the parent onto an allele that crosses it; under
  direct calling it is not genotyped. When the child has more than one candidate
  allele, vg also computes its likelihoods and direct call at the other ploidy, because the
  barrier may choose either; a child with one candidate allele keeps its provisional ploidy unless
  it is dropped. Each site is **staged**: vg keeps what the site's records will be built from, and
  writes nothing.
- In the **barrier**, vg settles and phases the sites one generation at a time, parents before
  children. The linkage model settles generation 0, and vg phases it from the panel (see
  [From the panel](#from-the-panel)). Each generation-1 site then takes its ploidy from its
  parent's settled genotype, and vg uses the likelihoods and direct call that the sweep computed
  at that ploidy. A site at ploidy 0 is dropped with everything nested in it. The model then
  settles generation 1, each linkage chain starting from the panel haplotypes that the parent's
  strands copy. Each later generation follows in the same way.

After the last barrier pass, vg **renders** each staged site: it builds and writes the site's
records from its settled genotype and phase.

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

A diploid site is **phaseable** when its settled genotype holds two different candidate alleles, so
that its two possible phases differ. A phaseable site can still be homozygous in the VCF, when its
two alleles differ only inside a child chain.

A **nested haploid chain** is a linkage chain of ploidy-1 sites whose parent is diploid, or is a
site of another nested haploid chain. All its sites lie on one strand of the nearest diploid
ancestor, the strand that carries them, directly or through ploidy-1 parents (see
[Which child chains are genotyped](#which-child-chains-are-genotyped)). Its sites' `GT` is `a|.` on
strand 0 and `.|a` on strand 1.

### From the panel

vg phases each linkage chain from its **Viterbi path**: the most probable sequence of the linkage
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

Each record's `FORMAT/PS` names its **phase set**: the reference position of the first site of its
top-level linkage chain. A nested site takes its parent's phase set. Sites of one phase set are
phased relative to one another.

At some phaseable sites, neither strand of the Viterbi path copies a panel haplotype that carries
either called allele. Nothing then orders the pair. vg still writes it with `|`, in an arbitrary
order, and the site stays in its phase set. Read phasing, when on, can still order it.

Phasing is on wherever the linkage model runs. `--phased` makes vg call fail when the linkage model
does not run (see [Direct call or linkage](#direct-call-or-linkage)), and `--no-phased` turns
phasing off. Read phasing, re-genotyping and `--anchors-hom-split` start from that phase, so an
explicit `--read-phasing`, `--regenotype` or `--anchors-hom-split` with `--no-phased` is an error,
and a preset's are turned off. Where the linkage model runs, `--no-phased` also turns nested calling
off (an explicit `--nested` is then an error), and variation inside nested sites is reported in the
enclosing site's alleles. Nested calling needs phasing there because the barrier takes each nested
site's strand, and the panel haplotypes its linkage chain starts from, from its parent's phase.

### From the reads

A read that spans two phaseable sites shows directly whether their alleles lie on the same strand.
`--read-phasing` uses such reads to re-decide the phase of each phaseable site, nested and
off-reference sites included, within the phase sets the panel gave. It leaves out the chains that a
[block record](#block-records) spells out. Read phasing never changes a genotype. It is off by
default and on under `--preset ont`.

Read phasing takes one phase set at a time, and its phaseable sites in order of position (see
[Transitions](#transitions)). A read is identified by its name, so the two mates of a pair are one
read. Where both mates are used at one site, they come from one molecule and often read the same
bases, so the site keeps only the mate with the higher confidence (defined below); between equal
confidences, the one with the larger $q_{rs}$, then the larger $c_{rs}$. Each site then holds each
read at most once, in its reliability, its links and its coherence alike.

#### What each read says

Take a read $r$ of a phaseable site $s$, and let $a_0$ and $a_1$ be the site's alleles on strand 0
and strand 1 in its current phase. For $k \in \lbrace 0, 1 \rbrace$, let
$x_{rsk} = (1 - e_r) v_{a_k} p_{r a_k}$, where $e_r$ is the read's
[mismapping probability](#mismapping-probability) and $p_{r a_k}$ its
[relative likelihood](#read-term) under $a_k$ at $s$. The **allele-length weights**
$v_{a_0}, v_{a_1}$ are the [mixture weights](#mixture-weights) with $\lvert a_k \rvert$ in place
of $U_i(G)$. $\lvert a \rvert$ is the length of all of allele $a$'s node visits, boundary nodes
included, while $U_i(G)$ counts only the nodes that one allele visits and the other does not. A
read with $x_{rs0} + x_{rs1} = 0$ is not used at $s$. For each read used, we keep:

- $q_{rs} = x_{rs0} / (x_{rs0} + x_{rs1})$, the probability that the read carries $a_0$, given
  that it came from one of the two strands;
- $c_{rs} = (x_{rs0} + x_{rs1}) / (x_{rs0} + x_{rs1} + e_r)$, the probability that it did come from
  one of them, rather than being mismapped;
- its **confidence**,
  $-10 \log_{10}\left(1 - \max(x_{rs0}, x_{rs1}) / (x_{rs0} + x_{rs1} + e_r)\right)$, the
  phred-scaled probability that its better allele is wrong.

#### Links between sites

Two sites $s$ and $t$ are **linked** by the reads they share. For a shared read $r$, let

$$
m_r = q_{rs} q_{rt} + (1 - q_{rs})(1 - q_{rt})
$$

This is the probability that the read's alleles at the two sites lie on one strand, given the
sites' current phases. The read is assumed to report this relation truly with probability
$\gamma_r = c_{rs} c_{rt}$, and to be a coin flip otherwise. The link is

$$
\mathrm{link}(s, t) = \sum_{r} \log_{10} \frac{\gamma_r m_r + (1 - \gamma_r)/2}{\gamma_r (1 - m_r) + (1 - \gamma_r)/2}
$$

A positive link favours the two sites' current phases. A link's **size** is its absolute value.
`--phase-cap`, when not 0, limits the size of each link.

#### Reliable sites

A site's **reliability** is the mean confidence of its reads, and the site is **reliable** if this
is at least `--phase-min-q`. Consider a read with $e_r = \epsilon_{\min}$, the floor on $e_r$, at a
site whose two alleles have equal length. If the read fits one allele perfectly and the other not at
all, its confidence is
$-10 \log_{10}\left(\epsilon_{\min} / (\epsilon_{\min} + (1 - \epsilon_{\min})/2)\right)$, the
**heterozygous score ceiling**. A read that favours the longer of two unequal alleles can have a
higher confidence. vg rejects a `--phase-min-q` above the ceiling, since sites whose alleles are of
similar length could not reach it.

#### Deciding the phases

To **flip** a site is to reverse its phase. For each site, read phasing decides whether to flip it
against the phase the panel gave. It phases the reliable sites first, each relative to the one
before it, and then phases every other site on its own. A wrong link between reliable sites flips
every site after it, while a wrong decision for another site flips only that site. Each phase set
goes through four stages:

1. **Phase chain.** The reliable sites of the phase set, in order of position, form its **phase
   chain**, and each is linked to the next. The phase chain breaks wherever the size of a link is
   below `--phase-break` $\log_{10}$ units, and the breaks divide it into **pieces**. The first site
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
   **agrees** at $s$ if $q_{rs}$ lies on the winning strand's side of $1/2$. Every $q$ here is taken
   in the current phase. A site's **coherence** is the fraction of its reads that agree, counting
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
   and a tie keeps the panel's phase.

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

A read's **strand log-odds** at site $s$, $\Lambda_{rs}$, measures how strongly its other sites
place it on strand 0. Each other phaseable site $t$ at which $r$ is used adds one term: the
natural-log odds that the read came from strand 0, with $q_{rt}$ taken in the phase read phasing
chose. A site adds one term per read, from the mate it keeps (see
[From the reads](#from-the-reads)).

$$
\Lambda_{rs} = \sum_{t \neq s} \ln \frac{c_{rt} q_{rt} + (1 - c_{rt})/2}{c_{rt}(1 - q_{rt}) + (1 - c_{rt})/2}
$$

Leaving out $s$ keeps a site from confirming its own genotype. Each phase set labels its strands
independently, so a read's strand is usable only at the sites of the phase set it was used in. Its
$\Lambda_{rs}$ is 0 at a site of another phase set, and at every site if it was used in more than
one phase set. At a site that has no phase set, any read used in one phase set keeps its
$\Lambda_{rs}$.

#### Tempering

The terms of $\Lambda_{rs}$ are not independent evidence. A wrong phase at some of the read's
sites, or an error that the read makes the same way at several sites, enters several terms at once.
So $\Lambda_{rs}$ overstates how sure the strand is. The **tempered** strand log-odds is

$$
y_{rs} = \mathrm{logit}\left(C \sigma(\tau \Lambda_{rs}) + (1 - C)/2\right)
$$

where $\tau$ is the **temper** (`--regeno-temper`) and $C$ is `--regeno-ceiling`. With $C = 1$,
$y_{rs} = \tau \Lambda_{rs}$. A smaller $C$ keeps the probability of either strand between
$(1 - C)/2$ and $(1 + C)/2$.

Unless `--regeno-temper` is given, $\tau$ is fitted once, from the phase that read phasing gave
before the first correction. Each read at each phaseable site $s$ gives an **observation** when
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
on the genotypes it settles. A correction, the barrier and read phasing together make a **round**.
The new genotypes can change the phase, and with it $\Lambda$, so another round can follow. Each
round corrects the likelihoods from the sweep, not those of the previous round. `--regeno-passes`
caps the number of times the barrier runs, the first run included. The rounds stop sooner when the
correction changes no site's direct call, when the barrier changes no settled genotype, or when the
settled genotypes return to an earlier state. Every round's correction is settled by the barrier,
including the round that stops, so the genotypes are settled from the likelihoods that `GL` reports.
With `--regeno-passes 1` the correction is computed and reported but not applied. `--regeno-ledger`
writes one line for each site whose direct call the last correction changed.

## Output

### VCF fields

| Field | Meaning |
|---|---|
| `GT` | the settled genotype, phased where phasing ran |
| `GL` | $\log_{10} \mathcal{L}(G)$ for every genotype of the record's alleles, in VCF order |
| `GQ` | the difference between the log-likelihoods of the direct call and the **runner-up**, the genotype with the second-highest $\mathcal{L}(G)$, in phred units. It is multiplied by the **explained share**, the fraction of reads whose best allele is in the direct call (a tied read split as for `AD`), and by the `--depth-quality` factor where that applies |
| `GQI` | the same difference, with neither factor |
| `GQN` | the same difference divided by the **achievable gap** (below), held at 1 or less, and multiplied by the explained share |
| `GP` | one value: the natural log of the posterior probability of the direct call, computed from $\mathcal{L}(G)$ with a uniform prior over genotypes. (In the VCF specification, `GP` is a phred-scaled value per genotype.) |
| `QUAL` | phred-scaled posterior probability, under the same uniform prior, of the genotype whose alleles are all the reference allele; 0 when `GT` is all reference |
| `DP` | number of reads of the site |
| `AD` | for each allele in the record, the number of reads whose best allele it is, rounded; a read tied between alleles counts a fraction to each |
| `DR` | $N_{\mathrm{eff}} / \mu_G$: the effective read count over the count that the written genotype predicts |
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

`GQN` removes both effects by dividing by the **achievable gap**. This is the difference that the
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

A record is **moved** when the linkage model settles a genotype other than the direct call. Its
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
them as the per-site `GQ` is, for that genotype, but with the sweep's `--depth-quality` factor. Each
round starts again from the sweep's likelihoods and `GQ`, so `GQ` follows the last round's
correction, and is the sweep's where that correction left the best genotype alone. `GP`, `GQI`,
`GQN` and `lowconf` keep their values from before the correction, so they describe the direct
call made from the uncorrected likelihoods. `DR` does not depend on the likelihoods, and describes
the written genotype. A moved record takes the `GQ`, `GQN` and `lowconf`
described above, whether or not re-genotyping ran.

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
that haplotype reaches back to the first. Where neither does, the stretch between them is a **gap in
the walk**. The gap is filled with the reference, on a `ref` row, if the reference is a panel
haplotype that crosses it. The linkage model can have a strand copy a panel haplotype at a site that
the haplotype does not pass through, so a segment's haplotype need not cover the whole segment. Such
a segment is replaced by a `ref` row where the reference crosses it.

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

A **pin** is a point between two adjacent bases of the graph, with no sequence of its own. Each site
has two. The **start pin** lies just after the site's start boundary node, where an allele's walk
enters the interior. The **end pin** lies just before its end boundary node, where the walk leaves.

A site's reads are divided among its **slots**. A slot's number names a field of `GT`: slot 0 the
first allele, to the left of the `|`, and slot 1 the second. A phaseable site has one slot for each
allele. Any other diploid site has a single slot, 0, unless `--anchors-hom-split` divides its reads
between slots 0 and 1, which then carry the same allele. A haploid site also has a single slot, 0,
except at a site of a nested haploid chain whose `GT` is `.|a`, where it is 1. An **anchor** is one
pin together with the reads of one slot that cross it.

#### Placing reads in slots

For each distinct called allele $a$, let $x_{ra} = (1 - e_r) v_a p_{ra}$, where $v_a$ are the
allele-length weights of [From the reads](#from-the-reads), and $v_a = 1$ at a site that is not
phaseable. Each read goes to the slot whose allele has the largest $x_{ra}$, except where the read
phase decides (below). A read's **anchor confidence** is
$-10 \log_{10}\left(1 - x_{ra} / (\sum_b x_{rb} + e_r)\right)$, where $a$ is the allele of the slot
it is placed in. The sum runs over distinct alleles, so the two slots of a split site count their
allele once. A site's **anchor reliability** is the mean anchor confidence of the reads written for
it, after the read filters below, each read counted once: paired mates share a name, and count at
the higher of their confidences.

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
  `--split-min-q`. A split site leaves out a read whose strand is not usable there, because the
  read was used in another phase set or in more than one: its strand says nothing about this
  phase set's. Any other read whose $y_{rs}$ is 0 is placed by a coin flip derived from its name,
  so it takes the same slot at every site.

#### Rows of the anchor file

The file is tab-separated. Each anchor is written as an `A` row, followed by one `R` row for each
of its reads. Header lines start with `#`, and a reader should skip keys it does not recognise.
Among them, the `#read` lines list the read names, which `R` rows refer to by index. The `#H` lines
name the columns of each kind of row, and the `#note` lines describe them.

| `A` column | Meaning |
|---|---|
| `node` | the graph's ID of the pin's boundary node: the first node of `snarl` for the start pin, the second for the end pin. Under `-N` or `-O`, `snarl` names the translated segments, while `node` stays the graph's ID |
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
| Relative likelihood | `--gap-open`, `--gap-extend`, `--insertion-nats`, `--optimal-pairing`, `--no-optimal-pairing` |
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

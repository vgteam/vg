# Direct genotyping

**Direct genotyping** is the part of `vg call --read-likelihood` that genotypes each site on its
own, from the site's reads. For each genotype $G$ that the site could have, it computes the
genotype's **site likelihood** $\mathcal{L}(G)$, a number that measures how well $G$ explains the
reads. The genotype with the highest site likelihood is the site's **direct call**. No prior over
genotypes is applied, so the direct call is the maximum-likelihood genotype (ties are broken as
[Direct genotype](read-likelihood-genotyping.md#direct-genotype) says).

Direct genotyping hands the rest of the caller the site likelihood of every genotype, the direct
call, and two kinds of value for each of the site's reads: its mismapping probability, and one
relative likelihood for each candidate allele. Both are defined below. The
[linkage model](read-likelihood-linkage-model.md) takes the site likelihoods as its evidence.
[Read phasing](read-likelihood-read-phasing.md#read-phasing),
[re-genotyping](read-likelihood-read-phasing.md#re-genotyping-from-the-phase) and the
[anchor file](read-likelihood-genotyping.md#assembly-anchors---anchors-out) take the per-read
values.

The caller as a whole is described in
[read-likelihood-genotyping.md](read-likelihood-genotyping.md), which also says how a site's
settled genotype follows from its direct call. Its [Vocabulary](read-likelihood-genotyping.md#vocabulary)
defines the terms used here, such as site, allele, ploidy, genotype, read placement and the reads of
a site.

## Model of how reads arise

We treat the reads of a site as generated from the sample's genotype in five steps.

1. **Choose the coverage.** The genotype, and the depth of sequencing near the site, set how many
   reads the site is expected to receive. Longer alleles are expected to receive more. This number
   is the genotype's **expected read count**. The site's **effective read count** is drawn from a
   distribution around the expected read count. It counts each read as its probability of being
   correctly mapped rather than **mismapped**, put in the wrong place by the mapper, so every read
   counts as less than one.
2. **Choose the reads and their mapping qualities.** The site receives some number of reads, each
   with a mapping quality (MAPQ), which sets the probability that the read is mismapped and so how
   much it adds to the effective read count. The number of reads and their MAPQs are drawn from a
   distribution that, given the effective read count, is independent of the genotype.
3. **Decide whether each read is mismapped.** Independently for each read, with the probability
   that its MAPQ sets, the read is mismapped. Either it came from elsewhere in the genome, or it
   came from this site but its placement gives it the wrong start, walk or edits within the site.
   Otherwise the read is correctly mapped.
4. **Choose the haplotype that each correctly mapped read came from.** Each correctly mapped read
   came from one of the sample's haplotypes at the site. The probability that it came from a given
   haplotype is that haplotype's **mixture weight**.
5. **Copy the read from that haplotype's allele.** The read's bases inside the site are a copy of
   part of the allele, with sequencing errors. The error model is the one behind vg's alignment
   scores: base substitutions have probabilities set by base quality, and insertions and deletions
   have fixed penalties rather than probabilities (see [Relative likelihood](#relative-likelihood)).
   A mismapped read's bases come from the place in the genome where the read originated.

Three of the steps leave evidence about the genotype in the reads. The effective read
count reflects the alleles' lengths: a homozygous deletion shows mainly as reads that are missing.
The chance of mismapping makes a read with a doubtful placement count for less. Most of the
evidence comes from the copying step, because a read fits the allele it was copied from better
than it fits the others.

The model rests on five assumptions. Each is listed with the data that break it, and with the
error that follows.

- **Reads are independent given the genotype.** Reads from one haplotype share its allele, and the
  genotype accounts for that. The assumption fails for paired mates, and for reads that share a
  systematic error. Their evidence is then counted more than once, and the likelihoods grow
  over-confident with depth.
- **Whether a read is mismapped depends only on its MAPQ.** The mismapping probability is computed
  from MAPQ, held between a floor and a ceiling (see
  [Mismapping probability](#mismapping-probability)). Because a mismapped read is also taken to fit
  as well as its best allele (see [Read term](#read-term)), a read that fits every allele badly is
  no more likely to be treated as mismapped than one that fits well. The assumption fails for a read
  from another copy of a repeat that the mapper placed here with a high MAPQ, which then counts as
  evidence against the alleles it fits worst.
- **MAPQs say nothing about the genotype beyond the effective read count.** This is step 2. The
  assumption fails where one allele duplicates sequence elsewhere in the genome, so that its reads
  get lower MAPQs than the other allele's. Those reads then count for less, and the evidence for
  that allele is understated.
- **The effective read count depends only on the genotype and on the depth near the site.** The
  expected read count uses the depth measured near the site, and the effective read count's
  distribution is wider than a Poisson distribution (see [Depth term](#depth-term)). The
  assumption fails where mappability, base composition or a change in copy number moves the depth
  at one site away from that of its neighbourhood. The depth term then favours genotypes whose
  allele lengths fit the wrong depth.
- **Each haplotype passes through the site once.** An allele can visit a node more than once: a
  walk that loops inside the site, between its boundary nodes, is one allele. The assumption fails
  at a duplication, where a haplotype's path passes through the site twice. The model still gives
  that haplotype one allele, and explains the reads of both copies with it, so the site has more
  reads than the genotype predicts, and the depth term favours genotypes with longer alleles.

## Notation

| Symbol | Meaning |
|---|---|
| $A$ | the site's candidate alleles: the alleles vg considers at the site (see [Candidate alleles](read-likelihood-genotyping.md#candidate-alleles)) |
| $a, b$ | alleles in $A$ |
| $P$ | the ploidy at the site |
| $G = \lbrace g_1, \dots, g_P \rbrace$ | a genotype; $g_i$ is the allele of its $i$-th haplotype, with the haplotypes listed in any order. A homozygous genotype lists the same allele $P$ times |
| $R$ | the reads of the site |
| $r$ | a single read in $R$ |
| $e_r$ | the mismapping probability of read $r$ |
| $e_R$ | the mismapping probabilities of all the reads in $R$ |
| $\epsilon_{\min}, \epsilon_{\max}$ | the floor and ceiling on $e_r$ |
| $\Pr(R \mid G)$ | the probability of the site's reads under $G$ |
| $\Pr(r \mid G, e_r)$ | the probability of read $r$'s bases under $G$, given its mismapping probability |
| $\Pr(r \mid a)$ | the probability of read $r$'s bases if it was correctly mapped and copied from allele $a$ |
| $\Pr(r \mid \text{mismapped})$ | the probability of read $r$'s bases if it was mismapped |
| $\mathcal{L}^\ast(G)$ | a model likelihood of $G$: a function of $G$ proportional to $\Pr(R \mid G)$ |
| $\mathcal{L}(G)$ | the site likelihood of $G$, which approximates a model likelihood |
| $p_{ra}$ | the relative likelihood of read $r$ under allele $a$ |
| $w_i(G)$ | the mixture weight of $G$'s $i$-th haplotype |
| $N_{\mathrm{eff}}$ | the effective read count: the sum of $1 - e_r$ over the site's reads |
| $\mu_G$ | the expected read count: the number of correctly mapped reads that $G$ predicts |
| $\mathrm{ExtPois}(\lambda, \beta)$ | the extended Poisson distribution, located at about $\lambda$, with its width set by $\beta$ |
| $f_{\mathrm{Pois}}$ | the Poisson probability, extended to counts that are not whole numbers |
| $f_{\mathrm{ExtPois}}$ | the density of $\mathrm{ExtPois}(\lambda, \beta)$ |
| $Z(\lambda, \beta)$, $\hat Z(\lambda, \beta)$ | the normaliser of $f_{\mathrm{ExtPois}}$, and the approximation of it that vg computes |
| $\hat f_{\mathrm{ExtPois}}$ | $f_{\mathrm{ExtPois}}$ with $Z$ replaced by $\hat Z$: the density that vg computes |
| $\beta$ | the width parameter, `--depth-term` |
| $T_a$ | the length in bases of allele $a$, excluding the site's two boundary nodes; a node the allele visits twice counts twice |
| $U_i(G)$ | at ploidy 2, the length in bases of the nodes that $g_i$ visits and the other allele of $G$ does not, each node counted once |
| $\bar L$ | the mean length of the reads that begin in the site's rate window (see [Depth term inputs](#depth-term-inputs)) |
| $\kappa$ | the read-start rate: the expected number of correctly mapped reads that begin at each base of one haplotype near the site |

## Likelihood formula

Under the model, the probability of the site's reads given $G$ has one factor for each part of the
model:

$$
\Pr(R \mid G) = \left[ \prod_{r \in R} \Pr(r \mid G, e_r) \right] \times \Pr(e_R \mid N_{\mathrm{eff}}) \times \Pr(N_{\mathrm{eff}} \mid G)
$$

The product over reads comes from steps 3 to 5, which happen independently for each read once its
MAPQ is chosen. $\Pr(e_R \mid N_{\mathrm{eff}})$ is the probability of the reads' particular MAPQs
given the effective read count, and comes from step 2. $\Pr(N_{\mathrm{eff}} \mid G)$ comes from
step 1.

We call a function of $G$ a **model likelihood** of $G$ if it is proportional to $\Pr(R \mid G)$
when the reads are held fixed. All model likelihoods give the same ratio between any two genotypes,
and so the same most likely genotype. $\Pr(e_R \mid N_{\mathrm{eff}})$ is the same for every
genotype, so

$$
\mathcal{L}^\ast(G) = \left[ \prod_{r \in R} \Pr(r \mid G, e_r) \right] \times \Pr(N_{\mathrm{eff}} \mid G)
$$

is a model likelihood. We compute the site likelihood as

$$
\mathcal{L}(G) = \left[ \prod_{r \in R} \left( (1 - e_r) \sum_{i=1}^{P} w_i(G) p_{r g_i} + e_r \right) \right] \times \hat f_{\mathrm{ExtPois}}\left(N_{\mathrm{eff}} ; \mu_G, \beta\right)
$$

The product over reads is the **read term**. It approximates $\prod_{r \in R} \Pr(r \mid G, e_r)$
divided by a value that depends on the reads but not on $G$ (see [Read term](#read-term)). The last
factor is the **depth term**. It approximates $\Pr(N_{\mathrm{eff}} \mid G)$, which is the density
of an extended Poisson distribution $\mathrm{ExtPois}$ (see [Depth term](#depth-term)). So

$$
\mathcal{L}(G) \approx \mathcal{L}^\ast(G) \Big/ \prod_{r \in R} \max_{b \in A} \Pr(r \mid b)
$$

The divisor is the same for every genotype, so the site likelihood approximates a model likelihood.
vg works with $\ln \mathcal{L}(G)$, the sum of the logarithms of the factors.

### Read term

Take one read $r$ of the site. By steps 3 to 5 of the model, it arose in one of $P + 1$ mutually
exclusive ways: it was correctly mapped and came from $G$'s $i$-th haplotype, for one $i$ from 1
to $P$, or it was mismapped. Its probability is the sum over these cases:

$$
\Pr(r \mid G, e_r) = \sum_{i=1}^{P} (1 - e_r) w_i(G) \Pr(r \mid g_i) + e_r \Pr(r \mid \text{mismapped})
$$

$\Pr(r \mid a)$ is the probability of the read's bases inside the site if the read was copied from
allele $a$. We compute it, up to a factor that depends only on the read, from an alignment of the
read to the allele (see [Relative likelihood](#relative-likelihood)).
[Mixture weights](#mixture-weights) and [Mismapping probability](#mismapping-probability) give
$w_i(G)$ and $e_r$.

A mismapped read's bases come from a sequence that the model leaves out: elsewhere in the genome,
or this site with a start, walk or edits other than its placement's. We approximate the
probability of its bases by its probability under its best allele at the site:

$$
\Pr(r \mid \text{mismapped}) \approx \max_{b \in A} \Pr(r \mid b)
$$

The mapper placed the read here because it resembles this site, so the sequence that the read
originated from probably explains it about as well as the site's best allele does. This choice
also makes the mismapped term $e_r$ after the division below, which bounds each read's
factor.

Dividing one read's probability by a value that does not depend on $G$ divides $\mathcal{L}^\ast(G)$
by that value, and so leaves a model likelihood. We divide each read's probability by
$\max_{b \in A} \Pr(r \mid b)$, which leaves

$$
\frac{\Pr(r \mid G, e_r)}{\max_{b \in A} \Pr(r \mid b)} \approx (1 - e_r) \sum_{i=1}^{P} w_i(G) p_{r g_i} + e_r
$$

The right-hand side is the read's **factor** in the read term. In it, $p_{ra}$ is the read's
**relative likelihood** under allele $a$:

$$
p_{ra} = \frac{\Pr(r \mid a)}{\max_{b \in A} \Pr(r \mid b)}
$$

A relative likelihood is 1 for the read's best allele, and smaller for alleles that explain the
read less well. The mixture weights of $G$, summed over $G$'s haplotypes, are 1, so
$\sum_{i} w_i(G) p_{r g_i}$ is the expected value of the read's relative likelihood over which of
$G$'s haplotypes it came from. The read's factor therefore lies between $e_r$ and 1.

Dividing by the best allele's probability does three things. We compute $\Pr(r \mid a)$ only up to a
factor that depends on the read, and that factor cancels in $p_{ra}$, so it need not be computed.
Each read's factor lies between $e_r$ and 1, so its logarithm is finite and does not underflow, for
a read of any length. And a read's relative likelihoods can be read directly: 1 for its best allele,
and near 0 for an allele that fits it far worse. `--dump-likelihoods` writes them.

### Mixture weights

The mixture weight $w_i(G)$ is the probability, before we look at the read's bases, that a
correctly mapped read of the site came from $G$'s $i$-th haplotype. Summed over $G$'s haplotypes,
the weights are 1. We approximate them crudely, from the lengths of the sequence that tells $G$'s
alleles apart.

The weights change the factors of only some reads. If a read fits $G$'s alleles equally well, so
that $p_{r g_1} = p_{r g_2}$, its factor is the same whatever the weights, because they sum to 1.
Call a read **informative** for $G$ if it fits one of $G$'s alleles better than another. Only the
factors of informative reads depend on the weights, so we choose each haplotype's weight to be its
expected share of the informative reads. The weights count only sequence that tells the alleles
apart, where the expected read count (see [Depth term](#depth-term)) counts all of it.

We approximate that share in two steps, at ploidy 2. First, we take a read from $g_i$ to be
informative only where it overlaps nodes that $g_i$ visits and the other allele of $G$ does not.
Where $g_i$ has no such nodes, as the reference allele has none against an insertion, a read from
$g_i$ is informative only where it spans the junction at which the other allele's extra sequence
would be. Second, we treat $g_i$'s unique nodes as one stretch, of length $U_i(G)$, which is 0 at a
junction. An informative read can then start at $U_i(G) + \bar L - 1$ positions on the haplotype,
and the weights are proportional to that number:

$$
w_i(G) = \frac{U_i(G) + \bar L - 1}{\sum_{j=1}^{P} \left(U_j(G) + \bar L - 1\right)}
$$

A homozygous genotype has $U_i(G) = 0$ for both haplotypes, which therefore get equal weights, as
do two alleles whose unique sequence is equally long, such as the two alleles of a SNP.
At ploidy 1, $w_1(G) = 1$.

The approximation errs in two cases. First, an allele's unique nodes can form several stretches,
separated by nodes that both alleles visit. When the stretches are further apart than a read is
long, each stretch adds $\bar L - 1$ start positions at its edge, and the formula adds $\bar L - 1$
once in all. Second, $U_i(G)$ counts each node once, in either orientation, so two alleles that
visit the same nodes get equal weights. At an inversion, whose alleles visit the same nodes in
opposite orientations, equal weights are right by symmetry. At a repeat whose two alleles go round a
loop different numbers of times, they are wrong. When the loop is longer than a read, every read
from the allele with fewer turns fits the other allele too, so all the informative reads come from
the allele with more turns.

`--flat-mixture` sets $w_i(G) = 1/P$ instead, so that the effect of the weighting can be measured.
It also flattens the allele-length weights that read phasing uses (see
[Read phasing](read-likelihood-read-phasing.md#read-phasing)).

### Mismapping probability

MAPQ is the mapper's estimate, on the phred scale, of the probability that the read came from
elsewhere in the genome, the first of the two ways of being mismapped in step 3. vg's mappers
compute it from the scores of the read's distinct placements in the graph, so it measures whether
the read belongs at another locus, not whether it is aligned correctly through this site. The
probability that the read is mismapped is the sum of the probabilities of the two ways. We
approximate the second by a constant, the floor $\epsilon_{\min}$ (`--mismap-min`), and the sum by
the larger of its two terms, and we hold the result below a ceiling $\epsilon_{\max}$
(`--mismap-max`), to give $e_r$:

$$
e_r = \min\left(\max\left(10^{-\mathrm{MAPQ}_r / 10}, \epsilon_{\min}\right), \epsilon_{\max}\right)
$$

The larger of two probabilities is at least half their sum, so before the ceiling this
approximation is at least half the probability it stands for.

Every read is mismapped with probability at least $\epsilon_{\min}$, so the model can explain any
read that fits $G$'s alleles badly as mismapped. One read can therefore change the ratio of two
genotypes' likelihoods by at most a factor of $1 / \epsilon_{\min}$, and the higher the floor, the
less any single read's fit counts. (The read also adds $1 - e_r$ to $N_{\mathrm{eff}}$, which moves
the depth term.)

The ceiling applies to reads whose MAPQ gives a mismapping probability above it, such as MAPQ 0.
Many mappers give MAPQ 0 to a read with several equally good placements. Its unclamped $e_r$ would
be 1, which would make the read's factor 1 under every genotype. The ceiling decides how much such a
read still counts. Such reads are used unless `--read-min-mapq` excludes them.

`--no-mismap-term` sets every $e_r$ to $\epsilon_{\min}$, whatever the read's MAPQ, as if every
read were well mapped, so that the contribution of the mismapping case can be measured. The
effective read count and $\kappa$ are computed from $e_r$ too (see
[Depth term inputs](#depth-term-inputs)), so under this option every read counts as
$1 - \epsilon_{\min}$ in them.

### Depth term

The depth term represents step 1 of the model, $\Pr(N_{\mathrm{eff}} \mid G)$: how well the
effective read count fits $G$. $N_{\mathrm{eff}}$ is usually not a whole number, so it cannot
follow a Poisson distribution. We define a continuous distribution for it instead, the **extended
Poisson distribution** $\mathrm{ExtPois}(\lambda, \beta)$. It has two parameters: $\lambda > 0$,
which locates it, and $\beta > 0$, which sets its width. (`--depth-term 0` leaves the depth term
out, so that its factor is 1.)

It is built from the Poisson probability, written for a real $n \geq 0$ as

$$
f_{\mathrm{Pois}}(n ; \lambda) = \frac{\lambda^{n} \exp(-\lambda)}{\Gamma(n + 1)}
$$

where the gamma function $\Gamma$ extends the factorial to real numbers: $\Gamma(n + 1) = n!$ for
a whole number $n$. For a whole number $n$, $f_{\mathrm{Pois}}(n ; \lambda)$ is the Poisson
probability of $n$ events when $\lambda$ are expected. The density of
$\mathrm{ExtPois}(\lambda, \beta)$ at $n \geq 0$ is

$$
f_{\mathrm{ExtPois}}(n ; \lambda, \beta) = \frac{f_{\mathrm{Pois}}(n ; \lambda)^{\beta}}{Z(\lambda, \beta)}, \qquad Z(\lambda, \beta) = \int_0^\infty f_{\mathrm{Pois}}(x ; \lambda)^{\beta} \mathrm{d}x
$$

where the normaliser $Z(\lambda, \beta)$ makes the density integrate to 1. At $\beta = 1$ the
distribution is a continuous analogue of the Poisson distribution. For large $\lambda$ it is close
to a normal distribution with mean $\lambda$ and variance $\lambda / \beta$, where a Poisson
distribution has variance $\lambda$, so a $\beta$ below 1 widens it. We widen it because read
depth varies between sites more than a Poisson count would, for reasons the model leaves out.

Step 1 of the model draws $N_{\mathrm{eff}}$ from $\mathrm{ExtPois}(\mu_G, \beta)$, with $\beta$
set by `--depth-term`. So

$$
\Pr(N_{\mathrm{eff}} \mid G) = f_{\mathrm{ExtPois}}(N_{\mathrm{eff}} ; \mu_G, \beta)
$$

and the depth term is this likelihood with $Z$ replaced by an approximation $\hat Z$:

$$
\hat f_{\mathrm{ExtPois}}(n ; \lambda, \beta) = \frac{f_{\mathrm{Pois}}(n ; \lambda)^{\beta}}{\hat Z(\lambda, \beta)}
$$

We compute $\hat Z$ from the normal approximation: for large $\lambda$,
$f_{\mathrm{Pois}}(n ; \lambda)$ is close to a normal density in $n$ with mean and variance
$\lambda$, so

$$
\ln Z(\lambda, \beta) \approx \ln \hat Z(\lambda, \beta) = \frac{1 - \beta}{2} \ln (2 \pi \lambda) - \frac{1}{2} \ln \beta
$$

Against numerical integration at $\beta = 0.1$, $0.5$ and $1$, $\ln \hat Z$ is within 0.13 of
$\ln Z$ for every $\lambda \geq 2$, and within 0.02 for $\lambda \geq 30$. Below $\lambda = 2$,
$\ln \hat Z$ falls without bound as $\lambda$ shrinks while $\ln Z$ does not, so there we use its
value at $\lambda = 2$. At $\beta = 1$, $\ln \hat Z = 0$. At $\beta < 1$, $Z$ grows with $\lambda$,
about as $\lambda^{(1 - \beta)/2}$, so the normaliser counts against the genotype with the larger
expected read count, by about $\frac{1 - \beta}{2} \ln (\mu_G / \mu_{G'})$ in $\ln \mathcal{L}$
between genotypes $G$ and $G'$. That is 0 between genotypes whose alleles have equal lengths, such
as at a SNP. At $\beta = 0.1$, when one genotype expects twice as many reads as the other, it is
0.31.

The expected read count is

$$
\mu_G = \kappa \sum_{i=1}^{P} \left(T_{g_i} + \bar L - 1\right)
$$

$T_{g_i} + \bar L - 1$ is the number of positions at which a read of length $\bar L$ can start on
the haplotype and still include some of the interior of $g_i$. When $g_i$ has no interior
($T_{g_i} = 0$), it is the number of positions from which a read reaches across the junction
between the two boundary nodes.

## Computing the terms

The read term needs, for each read of the site, its relative likelihoods and its mismapping
probability. The depth term needs the effective read count, and the read-start rate and mean read
length that set the expected read count. Those two come from the reads that begin near the site
(the mean read length falls back to the site's own reads if none do).

### Read input

Reads with MAPQ below `--read-min-mapq`, secondary alignments, and alignments with no placement
are discarded as they are read in. The reads come from one of three sources:

- `--gam` or `--gaf-reads`: a GAM or [GAF](static/GAF.md) file, loaded into memory.
- `--gam` with `--gam-index`: a GAM file sorted and indexed by `vg gamsort --index`. Reads are
  fetched from it as they are needed, one range of consecutive node IDs at a time; node IDs in a
  typical pangenome graph increase roughly along the genome.
- `--gaf-base`: a [GAF-base](https://github.com/jltsiren/gbz-base) database of alignments,
  queried one range of node IDs at a time by the `gbz-base` program. `--gaf-base-binary` gives the
  path to that program. It reads the graph from the GBZ-base database given with `--gbz-base`, or
  else from the input graph.

`--read-window` sets the size of those ranges, in node IDs, for the two indexed sources. It changes
which reads are fetched together and the order in which they arrive. The results are the same for
every value of it. A floating-point sum depends on its order, so a site's reads are put in a fixed
order, by read name, then by where the alignment begins, then by their values, before anything is
summed over them. The reads that begin in a rate window (see
[Depth term inputs](#depth-term-inputs)) are tallied as counts per MAPQ, which do not depend on
order.

A read whose alignment crosses the site against the direction of the alleles is reverse-complemented
before it is scored. We decide the direction by a vote over the read's visits to the nodes that the
alleles visit, leaving out any node that two alleles visit in opposite orientations. The read is
reversed if more of those visits are in the opposite orientation to the alleles' than in the same
one. A tie, including a read that visits none of those nodes, leaves it as it is. The boundary nodes
decide the vote at an inversion, whose inverted nodes the alleles visit in both orientations.

### Relative likelihood

We compute a read's relative likelihoods by aligning the read to each candidate allele, scoring
each alignment, and comparing each score with the best one. The alignment is of node visits, not
of bases. The read's placement and the allele are both sequences of node visits, and two visits are
the **same visit** when they go to the same node in the same orientation. We align the read's
visits inside the site, boundary nodes included, to the allele's visits in much the way that
dynamic programming aligns two DNA sequences with affine gap costs, with node visits in place of
bases.

Such an alignment is a **pairing**. Each read visit is either paired with one allele visit or left
unpaired, and the pairs come in the same order in both sequences. A pair is one of two kinds:

- **A read visit paired with the same allele visit.** It scores the read's own edits in that node,
  as the mapper aligned them.
- **A substitution.** It pairs a read visit that the allele does not make with an allele visit that
  the read does not make. It compares the read's bases in its node with the allele node's sequence,
  base by base from the first base of each, over the shorter of the two lengths, and adds a gap
  for the difference in length. This is an approximation: it compares a read that starts partway
  into its node, or two nodes that differ by an internal indel, out of register.

So a visit that the read and the allele share is either paired with the same visit or left
unpaired. Visits left unpaired are gaps, measured in bases:

- A run of consecutive read visits left unpaired is an **insertion**, scored as one gap as long as
  their bases. A visit that the read and the allele share may be left unpaired, which can explain
  a read with poor edits in that node better. Before the read's first pair of same visits, each
  unpaired read visit is a gap of its own. A read visit with no bases, where the read deletes its
  whole node, adds nothing when left unpaired and does not break a run.
- A run of allele visits left unpaired after the read's first pair of same visits and before its
  last pair is a **deletion**, scored as one gap as long as their bases. Unpaired allele visits
  before the read's first pair of same visits, or after its last pair, lie outside the read and
  score 0, so that a read is not penalised for being short.

Every read base inside the site is therefore scored under every allele, and all alleles are scored
over the same read bases.

Scores use vg's quality-adjusted alignment scoring, in which a mismatch at a low-quality base costs
less, with gap scores set by `--gap-open` and `--gap-extend`. A read without base qualities is
scored without the quality adjustment. Scores are higher for better fits. A pairing determines a
base-level alignment of the read's bases in the site to the allele's sequence, and the pairing's
score $s_{ra}$ is the score of that alignment, except that each pair and each gap is scored on its
own, so a gap at the edge of one is not joined to a gap in the next.

#### Optimal pairing

With `--optimal-pairing`, we use optimal pairing: the pairing with the highest score $s_{ra}$,
found by dynamic programming over the read's visits against the allele's, as in affine-gap
alignment. Only the pairing of visits is searched; bases are not aligned again. When the product
of the numbers of read visits and allele visits exceeds a fixed limit, the dynamic programming is
banded (see [Fixed constants](read-likelihood-genotyping.md#fixed-constants)), so at the largest
sites the pairing found may not be the highest-scoring one.

#### Greedy pairing

Without `--optimal-pairing`, we use greedy pairing, which approximates optimal pairing in one pass
along the read's visits. It pairs each read visit with the next matching allele visit: the next
occurrence of the same visit after the last pair. Allele visits skipped over between two such pairs
are a deletion. When there is no matching visit, it looks for a simple insertion: if the allele's
next unpaired visit is one that the read makes later, the read visit is left unpaired. Otherwise it
pairs the read visit with the allele's next unpaired visit as a substitution, if neither of the two
visits is made elsewhere by the other sequence, and leaves the read visit unpaired if one is. Once
the allele's visits are used up, the remaining read visits are left unpaired. Each pair is final
once made.

Greedy pairing considers the same pairs as optimal pairing and scores the pairing it finds by the
same rules, including those for runs of unpaired read visits, so the two differ only in which
pairing they find.

#### From scores to relative likelihoods

We convert the pairing's score to a log-likelihood score in nats:

$$
\ell_{ra} = \alpha s_{ra} + \iota I_{ra}
$$

$\alpha$ is the scale that converts score units to nats (vg calls it the scorer's log base). vg
computes it from the substitution matrix so that $\alpha$ times a substitution score is the natural
log of a likelihood ratio: the probability of the aligned bases under the alignment model, over
their probability as unrelated random sequence. This is how vg interprets alignment scores as
probabilities elsewhere. The scorers with and without the quality adjustment have different log
bases. $I_{ra}$ is the number of insertions in the pairing in which the read has bases the allele
lacks: those inside the mapper's edits, each unpaired read visit with bases, and each substitution
whose read node is the longer. Each unpaired read visit counts once here, even where a run of them is scored as
one gap, so $I_{ra}$ depends on how many nodes an inserted sequence spans. $\iota$ is
`--insertion-nats`. A positive $\iota$ makes an insertion, where the read has bases the allele
lacks, cost less than a deletion of the same length, where the allele has bases the read lacks.
Optimal pairing chooses the pairing by $s_{ra}$ alone, and then adds $\iota$ for its insertions.

These scores stand in for the error model of the copying step. Substitution scores are calibrated
log-likelihood ratios; gap scores and $\iota$ are penalties, not normalised probabilities. Given the
alignment that the pairing describes, the probability of the read's bases in the site is taken to be
$\exp(\ell_{ra})$ times $B_r$, their probability as unrelated random sequence. $B_r$ is the same for
every allele, because every allele is scored over the same read bases. We take the probability of
one alignment, the best pairing found, in place of the sum over all alignments of the read to the
allele, so

$$
\Pr(r \mid a) \approx B_r \exp\left(\ell_{ra}\right)
$$

We score the read against every candidate allele, and then divide by the best. Substituting this
approximation into the definition of $p_{ra}$, $B_r$ cancels:

$$
p_{ra} = \frac{\Pr(r \mid a)}{\max_{b \in A} \Pr(r \mid b)} \approx \frac{B_r \exp\left(\ell_{ra}\right)}{B_r \exp\left(\max_{b \in A} \ell_{rb}\right)} = \exp\left(\ell_{ra} - \max_{b \in A} \ell_{rb}\right)
$$

### Depth term inputs

The depth term compares the effective read count $N_{\mathrm{eff}}$ with the expected read count
$\mu_G = \kappa \sum_{i} (T_{g_i} + \bar L - 1)$, under a distribution whose width $\beta$ sets.
Its inputs are:

1. **The read-start rate $\kappa$**, the expected number of correctly mapped reads that begin at
   each base of one haplotype near the site. It is measured over a **rate window** on the reference
   path. The reference path is cut into buckets of a fixed length (see
   [Fixed constants](read-likelihood-genotyping.md#fixed-constants)). A site belongs to the bucket
   that holds the reference position of its start boundary node. If that node has none, vg uses the
   end boundary node, and failing that, the nearest enclosing site. The site's rate window is its
   bucket and one bucket on each side. The window's **read rate** is the number of reads whose
   alignment begins on a reference node in the window, each counted as $1 - e_r$, divided by the
   reference length of the window. Only reference nodes count, in both the number and the length, so
   the rate is independent of how many non-reference nodes the window holds. Variation that the
   sample carries in the window still changes it: a deletion leaves reference bases on which no
   reads start. The counts are computed once per bucket and shared, and each site divides the read
   rate by the ploidy of its region (set by `--ploidy`, `--ploidy-regex` or `--ploidy-bed`) to give
   $\kappa$, a rate per haplotype. That is the site's own ploidy, except at a nested site that only
   some of its parent's alleles cross: the site is genotyped at a lower ploidy, but the window's
   reads come from every haplotype. When no read begins in the window, $\kappa = 0$, and the depth
   term is left out (its factor is 1). A site with no reference position anywhere among its
   enclosing sites, as in a graph without reference path positions, uses a fixed block of
   consecutive node IDs in place of the rate window, the block that contains the site's lowest node
   ID, with the summed length of its nodes as the length.
2. **The mean read length $\bar L$**, the unweighted mean sequence length of the reads that begin in
   $\kappa$'s rate window, or of the site's own reads if none does. The mixture weights use the same
   $\bar L$.
3. **The effective read count** $N_{\mathrm{eff}} = \sum_{r \in R} (1 - e_r)$, over the site's
   reads.
4. **The width parameter $\beta$**, `--depth-term`.

The expected read count $\mu_G$ is then computed for each genotype from $\kappa$, $\bar L$ and the
lengths of the genotype's alleles.

`--depth-count-raw` counts each read as 1 in place of $1 - e_r$, in both $N_{\mathrm{eff}}$ and
$\kappa$, so that the contribution of the mismapping probabilities to the depth term can be
measured. $N_{\mathrm{eff}}$ is then the number of the site's reads, a whole number.

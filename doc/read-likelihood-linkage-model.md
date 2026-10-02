# Linkage model

The **linkage model** is the second stage of `vg call --read-likelihood`'s genotyping. It chooses
the genotypes of neighbouring sites together, using the panel haplotypes. Alleles at nearby sites
are inherited together, so the sample's alleles at neighbouring sites tend to form combinations
that panel haplotypes also carry. The model uses this by treating each of the sample's strands as
a mosaic of panel haplotypes, as the
[Li–Stephens model](https://doi.org/10.1093/genetics/165.4.2213) does.

At each site the model uses the site likelihood $\mathcal{L}(G)$ of every genotype $G$, which
measures how well $G$ explains the site's reads; the allele that each panel haplotype carries there;
and the site's position. Direct genotyping computes the site likelihoods from each site's reads
alone, and the genotype with the highest is the site's **direct call**. The model runs along
**linkage chains**: sequences of sites in order of position, each holding the top-level sites of a
contig between changes of ploidy, or the sites of one child chain at one ploidy (and, at ploidy 1,
on one strand of the parent). Each time it runs, in a
[linkage pass](read-likelihood-genotyping.md#passes-and-rounds), it gives each site of a linkage
chain a **chosen genotype**, and it phases the chosen genotypes. The genotype it chooses in the last
pass is the site's **settled genotype**, the genotype that vg reports.

The caller as a whole is described in
[read-likelihood-genotyping.md](read-likelihood-genotyping.md), whose
[Vocabulary](read-likelihood-genotyping.md#vocabulary) defines terms used here, such as site,
strand, panel and phase. The site likelihood is derived in
[read-likelihood-direct-genotyping.md](read-likelihood-direct-genotyping.md).

## Differences from PanGenie

The model follows PanGenie ([Ebler et al. 2022](https://doi.org/10.1038/s41588-022-01043-w)). In
both, a hidden Markov model runs along a sequence of sites. The hidden state at a site is an
ordered pair of panel haplotypes, one for each of the sample's strands, and it implies a genotype:
the alleles that the two panel haplotypes carry there. Between sites, each strand can switch to
copying another panel haplotype. The forward–backward algorithm gives each state's posterior
probability, and a genotype's posterior is collected from the states that imply it. The Viterbi
algorithm gives the phase. This model differs from PanGenie's as follows.

- **Emissions.** PanGenie computes a state's emission from the counts, in the reads, of k-mers
  unique to the site. Here the emission is the site likelihood of the genotype that the state
  implies, which direct genotyping computes from the reads' alignments to the graph (see
  [States and emissions](#states-and-emissions)).
- **Wildcard and escape.** PanGenie can call only alleles that some panel haplotype carries. This
  model adds a **wildcard** haplotype, which can carry any allele, so that it can call an allele
  that no panel haplotype carries at the site. A strand that copies the wildcard, or a panel
  haplotype that does not pass through the site, has an unknown allele there. Each such strand
  multiplies the state's emission by a fixed **escape** penalty (see
  [Wildcard haplotype](#wildcard-haplotype)).
- **Switch probability.** PanGenie, as Li and Stephens did, redraws a strand's haplotype with
  probability $1 - \exp(-4 N_e r d / K)$, from the distance $d$, a recombination rate $r$, an
  effective population size $N_e$ and the panel size $K$. This model keeps the form
  $1 - \exp(-d/D)$, but sets the scale $D$ directly, with `--linkage-scale`, so that it does not
  change with $K$. It adds a floor, and raises the result to a power, the **linkage weight**
  (`--linkage-weight`). A redraw picks uniformly from the panel and the wildcard (see
  [Transitions](#transitions)).
- **Allele-frequency prior.** In both models, a genotype that many pairs of panel haplotypes carry
  is implied by many states, and so collects more posterior probability. PanGenie keeps this prior
  as the states imply it. This model raises it to a power, `--linkage-prior` (or `--hp-prior` at
  long homopolymers), which can strengthen or weaken it (see
  [Allele-frequency prior](#allele-frequency-prior)).
- **Windows.** PanGenie runs forward–backward along each chromosome as a whole. This model decodes
  a long linkage chain in overlapping windows of a fixed number of sites (see
  [Forward–backward windows](#forwardbackward-windows)).
- **Ploidy 1.** PanGenie's states are always pairs. Where the sample has one strand, the states
  here are single haplotypes, a panel haplotype or the wildcard: on a contig or region of ploidy 1,
  and at a nested site that only one strand passes through.
- **Nested sites.** PanGenie treats each top-level bubble as one site, with the variation nested in
  it folded into its alleles. This model also runs on sites nested inside other sites. It decodes
  parents before children, and starts each child chain from the panel haplotypes that the parent's
  strands copy (see [Linkage chains](#linkage-chains)).
- **Phasing.** PanGenie's Viterbi path is unrestricted, so the genotypes it implies need not be the
  posterior ones. Here the Viterbi path is restricted to states compatible with each site's
  chosen genotype, so phasing orders each genotype's alleles and keeps the genotype (see
  [Phasing from the panel](#phasing-from-the-panel)).

## Notation

| Symbol | Meaning |
|---|---|
| $K$ | the number of panel haplotypes |
| $\ast$ | the wildcard haplotype |
| $h$ | a panel haplotype or the wildcard |
| $(h_0, h_1)$ | a state at ploidy 2: strand 0 copies $h_0$, and strand 1 copies $h_1$ |
| $a(h)$ | the allele that panel haplotype $h$ carries at the site, when it passes through the site |
| $A_{\mathrm{c}}$ | the site's compact allele set |
| $\mathcal{L}(G)$ | the site likelihood of genotype $G$, from direct genotyping |
| $E(h_0, h_1)$, $E(h)$ | the emission of a state at ploidy 2, and at ploidy 1 |
| $\epsilon_{\mathrm{esc}}$ | the escape penalty |
| $d$ | the distance in bases between two consecutive sites |
| $\rho(d)$ | the probability that a strand redraws its haplotype between two sites $d$ bases apart (a redraw can pick the same haplotype) |
| $\rho_{\min}$ | the floor on the bracketed term of $\rho(d)$ |
| $D$ | the distance at which $1 - e^{-d/D}$ reaches $1 - 1/e$, `--linkage-scale` |
| $\omega$ | the linkage weight, `--linkage-weight` |
| $F$ | the exponent on the allele-frequency prior, `--linkage-prior`, or, when not 0, `--hp-prior` at a homopolymer site |
| $c_G$ | the number of states with no unknown allele that imply genotype $G$ |
| $n_a$ | the number of panel haplotypes that carry allele $a$ |

## States and emissions

At ploidy 2 the hidden state at a site is an ordered pair $(h_0, h_1)$: strand 0 copies $h_0$
there, and strand 1 copies $h_1$. Each of $h_0$ and $h_1$ is a panel haplotype or the wildcard, and
the two can be the same. At ploidy 1 the state is a single $h$.

A panel haplotype carries at most one allele at a site. It can take the walks of more than one
candidate allele there, when it passes through the site twice or is stored as several paths. vg
then takes it to carry the one that comes last in the site's list of candidate alleles (see
[Candidate alleles](read-likelihood-genotyping.md#candidate-alleles)), and drops the others.

At each site, the model works over the site's **compact allele set** $A_{\mathrm{c}}$, a subset of
its candidate alleles: the alleles of the site's direct call, and every allele that some panel
haplotype carries there. The model can choose for a site only a genotype made of these alleles. A
site whose compact allele set is larger than a fixed limit (see
[Fixed constants](read-likelihood-genotyping.md#fixed-constants)) is left out of the model. It
keeps its direct call, written unphased, and the model links the sites on either side of it
directly.

A state whose strands both have known alleles implies a genotype, made of those alleles. The
state's **emission**, the evidence that the reads give for it, is that genotype's site likelihood:

$$
E(h_0, h_1) = \mathcal{L}\left(\lbrace a(h_0), a(h_1) \rbrace\right)
$$

and at ploidy 1, $E(h) = \mathcal{L}(\lbrace a(h) \rbrace)$.

## Wildcard haplotype

The model adds a **wildcard** haplotype $\ast$ to the panel, so that it can call an allele that no
panel haplotype carries. The wildcard can carry any allele of the compact allele set. Its allele
is averaged out independently at each site, so a strand that stays on the wildcard links nothing
between sites.

A strand has an **unknown allele** at a site when it copies the wildcard, or a panel haplotype that
does not pass through the site. The emission of a state with one unknown strand is the mean of
$\mathcal{L}(G)$ over the alleles of the compact allele set, each taken in turn as the unknown
strand's allele, times a fixed **escape** penalty $\epsilon_{\mathrm{esc}}$. When strand 1's
allele is unknown,

$$
E(h_0, h_1) = \frac{\epsilon_{\mathrm{esc}}}{\lvert A_{\mathrm{c}} \rvert} \sum_{b \in A_{\mathrm{c}}} \mathcal{L}\left(\lbrace a(h_0), b \rbrace\right)
$$

and likewise when strand 0's is. With two unknown strands, the mean is over ordered pairs of these
alleles, and the escape penalty is applied once for each strand:

$$
E(h_0, h_1) = \frac{\epsilon_{\mathrm{esc}}^2}{\lvert A_{\mathrm{c}} \rvert^2} \sum_{b \in A_{\mathrm{c}}} \sum_{b' \in A_{\mathrm{c}}} \mathcal{L}\left(\lbrace b, b' \rbrace\right)
$$

At ploidy 1, the emission of a strand with an unknown allele is

$$
E(h) = \frac{\epsilon_{\mathrm{esc}}}{\lvert A_{\mathrm{c}} \rvert} \sum_{b \in A_{\mathrm{c}}} \mathcal{L}\left(\lbrace b \rbrace\right)
$$

## Transitions

Between consecutive sites $d$ bases apart, each strand independently keeps copying the same
haplotype or, with probability $\rho(d)$, draws one uniformly from the $K$ panel haplotypes and the
wildcard, possibly the same one. A strand that copies $h$ at one site therefore copies $h'$ at the
next with probability

$$
(1 - \rho(d)) [h' = h] + \frac{\rho(d)}{K + 1}
$$

where $[h' = h]$ is 1 when $h'$ is $h$ and 0 otherwise. The switch probability is

$$
\rho(d) = \left[ \rho_{\min} + (1 - \rho_{\min}) \left(1 - e^{-d/D}\right) \right]^{\omega}
$$

$d$ is the difference between two sites' positions (see
[Vocabulary](read-likelihood-genotyping.md#graph)), so a site's own length counts in its distance
to the next. $D$ is `--linkage-scale`. $\omega$ is `--linkage-weight`, an exponent on the switch
probability: a larger $\omega$ makes switches rarer and linkage stronger. `--linkage-weight 0`
turns the model off: every site keeps its direct call and is not phased from the panel. (The
formula at $\omega = 0$ would instead give $\rho = 1$.) $\rho_{\min}$ is a small fixed floor on the
bracketed term, so that $\rho(d) \geq \rho_{\min}^{\omega} > 0$ and every switch has a positive
probability.

The sites of an
[off-reference chain](read-likelihood-genotyping.md#which-child-chains-are-genotyped) have no
reference position of their own. vg takes an allele of the parent's chosen genotype that crosses
such a site, the one first in the parent's compact allele set where both do, and places the site
at the parent's position plus the site's offset along that allele. Two sites of one off-reference
chain are then as far apart as they are along that allele, and the chain's first site is as far
from the parent as its offset along that allele. No distance is known between such a site and any other
site with a reference position, and there $\rho = 1$.

## Genotype posterior

The forward–backward algorithm gives each state's posterior probability at each site. A genotype's
posterior collects the probability of the states that imply it:

- A state with no unknown strand implies one genotype.
- A state with one unknown strand, the other carrying allele $a$, shares its probability among the
  genotypes $\lbrace a, b \rbrace$, for each allele $b$ of the compact allele set, in proportion to
  their site likelihoods.
- A state with two unknown strands shares its probability among all genotypes of the compact
  allele set, each counted once, in proportion to their site likelihoods.

At ploidy 1 a state is a single haplotype, and a state with an unknown strand shares its
probability among the alleles in the same way.

## Allele-frequency prior

A genotype that many pairs of panel haplotypes carry is implied by many states, so it collects more
probability. The rescaling below is applied to each site's posteriors after forward–backward, so
it changes that site's chosen genotype but not its neighbours' posteriors or the Viterbi path. The
exponent $F$, `--linkage-prior`, sets the prior's strength:

- The probability from states with both alleles known is multiplied by $c_G^{F-1}$, where $c_G$ is
  the number of such states that imply $G$.
- The probability from states with one unknown strand is multiplied by $n_a^{F-1}$, where $n_a$ is
  the number of panel haplotypes, the wildcard aside, that carry $a$.
- States with two unknown strands are not rescaled.

At ploidy 1, the probability of the states that carry allele $a$ is multiplied by $n_a^{F-1}$, and
states with an unknown allele are not rescaled. The posteriors are then normalised to sum to 1.
$F = 1$ keeps the prior that the states imply, $F = 0$ removes it, and $F > 1$ strengthens it. The
model chooses the genotype with the highest posterior.

## Homopolymer sites

Sequencing errors in a long homopolymer run tend to recur in many reads at the same site, and the
site likelihood counts each such read as independent evidence. A larger exponent $F$ at a site like
this gives the panel's allele frequencies more weight against these reads.

`--hp-prior`, when not 0, replaces $F$ at **homopolymer sites**. At a homopolymer site, the
reference allele and some other candidate allele differ only in the length of one homopolymer run,
by at most a fixed number of copies of its base (see
[Fixed constants](read-likelihood-genotyping.md#fixed-constants)). The run must also be long: in
the longer of the two alleles it has at least `--hp-prior-run` bases, or it reaches an end of the
allele's sequence between its boundary nodes. A run that reaches an end can continue into a
boundary node, so its full length is unknown, and it counts as long.

## Forward–backward windows

vg runs the forward–backward algorithm over overlapping windows of a linkage chain. Each window
**keeps** a fixed number of consecutive sites, and the kept sites of successive windows cover the
chain without overlapping. A window also decodes a fixed margin of sites on each side of its kept
sites, fewer at the ends of the chain, and discards their posteriors. The windows are decoded
independently of one another. Every kept posterior therefore has the margin's number of sites on
each side, except near the ends of the linkage chain, and approximates the posterior from decoding
the whole chain at once. The first window of a chain below the top level starts as
[Linkage chains](#linkage-chains) describes, and every later window starts from a uniform
distribution. A linkage chain at ploidy 2 that starts from its parent's panel haplotypes is
decoded as one window, whatever its length.

## Linkage chains

The model runs on each linkage chain separately, one level of nested sites at a time, parents
before children (see [The linkage pass](read-likelihood-genotyping.md#the-linkage-pass)). A chain
below the top level is decoded after its parent's chain has been decoded and phased. The model then
starts the chain from the state that the parent's Viterbi path chose at the parent site (see
[Phasing from the panel](#phasing-from-the-panel)); read phasing, which comes later, does not change
it. A chain at its parent's ploidy holds the parent as its first site, fixed at the parent's chosen
genotype, with all its probability on that state. A ploidy-1 chain under a diploid parent does not
hold the parent. It starts from the panel haplotype $h$ that the parent's strand carrying the chain
copies, moved by one transition over the distance $d$ from the parent to the chain's first site:
probability $1 - \rho(d) + \rho(d)/(K + 1)$ on $h$ and $\rho(d)/(K + 1)$ on each other state, so
that the first site's reads can overrule the parent's haplotype. The chain starts from a uniform
distribution instead when the parent was not phased, or, at ploidy 1, when that strand copies the
wildcard or a panel haplotype that does not pass through the chain's first site.

## Phasing from the panel

vg phases each linkage chain from its **Viterbi path**: the most probable sequence of the model's
states along the chain, found window by window (below), among the states **compatible** with each
site's chosen genotype: a state with both alleles known that implies it, a state with one unknown
strand whose other strand carries one of its alleles, and a state with two unknown strands. Phasing
orders each chosen genotype's alleles and keeps the genotype; where the path's state has no known
allele, the order of a heterozygous site's alleles is arbitrary (see
[Phasing](read-likelihood-genotyping.md#phasing)). For the path, each compatible state's emission is
the chosen genotype's $\mathcal{L}(G)$, times $\epsilon_{\mathrm{esc}}$ for each unknown strand.
Every compatible state with both alleles known therefore has the same emission, so the path is
chosen by the transitions, with each unknown strand penalised by $\epsilon_{\mathrm{esc}}$. Read
phasing, when on, re-decides the phases afterwards (see
[Read phasing](read-likelihood-read-phasing.md#read-phasing)).

A site decoded at an earlier level, such as a parent decoded with its child chain, is held to the
state it was phased with, where that state is compatible with its chosen genotype, and
otherwise only to its genotype. On every linkage chain, the path is decoded over overlapping
windows of the sizes given in [Forward–backward windows](#forwardbackward-windows). Each window
after the first is held, at the previous window's last kept site, to the state that the previous
window chose there.

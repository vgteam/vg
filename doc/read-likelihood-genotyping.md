# Read-likelihood genotyping (`vg call --read-likelihood`)

`vg call --read-likelihood` genotypes one sample from its reads aligned to a pangenome graph. At
each variant site it computes, for every candidate genotype, the probability of the sample's reads
there, and calls the genotype with the highest probability. When the graph also stores a panel of
haplotypes, the per-site probabilities are instead combined along each chromosome by a hidden
Markov model over the panel, which also phases the calls. This page describes the model and names
the options that control each part. Default values are listed by `vg call --help`; the few
constants that have no option are listed under [Fixed constants](#fixed-constants), with where
they are defined.

A minimal run, on a GBZ graph (a vg graph file that also stores haplotypes) and reads aligned by
`vg giraffe` in GAM format (vg's alignment format):

```
vg call graph.gbz --read-likelihood --gam reads.gam > calls.vcf
```

`vg call` finds the graph's sites itself, or reads them from a file made by `vg snarls` (`-r`).

- [Vocabulary](#vocabulary)
- [The model](#the-model)
- [Notation](#notation)
- [The likelihood](#the-likelihood)
- [How each term is computed](#how-each-term-is-computed)
- [Calling a genotype](#calling-a-genotype)
- [The linkage model](#the-linkage-model)
- [Nested sites](#nested-sites)
- [Phasing](#phasing)
- [Output](#output)
- [Options](#options)

## Vocabulary

### The graph

- **Node.** A node holds a DNA sequence and can be traversed in either *orientation*: forward,
  reading its sequence, or reverse, reading its reverse complement. A *walk* is a sequence of
  oriented node visits along the graph's edges.
- **Site.** A snarl of the graph: a subgraph separated from the rest of the graph by two
  *boundary nodes*, a start and an end. The boundary nodes belong to the site, and its other nodes,
  including those of any sites nested in it, are its *interior*. A site is read from its start
  boundary node to its end boundary node.
- **Chain.** A series of sites joined end to end, each site's end boundary node being the next
  site's start. Sites nest: a site can contain chains of smaller sites, its *child chains*.
- **Allele.** A walk through a site from its start boundary node to its end boundary node, also
  called a *traversal*.
- **Reference.** The *reference paths* are the graph paths on which `vg call` reports positions:
  all paths of REFERENCE or GENERIC sense, or the paths chosen with `-p`, `-P` or `-S`. A site's
  *reference allele* is the walk a reference path takes through it. In the VCF the reference allele
  is allele 0, and the other alleles are numbered from 1.

### The sample and its reads

- **Haplotype.** One copy of the genome. The sample has $P$ haplotypes, where $P$ is its *ploidy*:
  1 or 2, and it can differ between regions (see [Ploidy](#ploidy)). At each site, each of the
  sample's haplotypes carries one allele.
- **Panel.** The haplotypes stored in the graph's GBWT, the haplotype index inside a GBZ file (or a
  separate GBWT file given with `-g`). A panel haplotype can be stored as several paths that each
  cover part of a chromosome, so it does not necessarily pass through every site.
- **Genotype.** The multiset of the $P$ alleles that the sample's haplotypes carry at a site, such
  as $\lbrace 0, 1 \rbrace$. Which haplotype carries which allele is the genotype's *phase*, which
  is decided separately. `vg call` genotypes one sample per run.
- **Read placement.** Reads are used as a mapper, such as `vg giraffe`, aligned them to the graph.
  A read's *placement* is
  the walk its alignment takes, the *edits* inside each node (the runs of matching bases,
  mismatches, insertions and deletions against the node's sequence), and its mapping quality
  (MAPQ).
- **Reads of a site.** The reads whose placement visits an interior node of the site, or moves
  between two different nodes of the site. The second case includes a read that crosses directly
  from the start boundary node to the end boundary node, which supports an allele that deletes the
  site's interior. A read that stays inside one boundary node fits every allele equally and is not
  used.

## The model

We treat the reads at a site as generated from the sample's genotype $G$ by this process:

1. **How many reads.** The number of reads that come from the site is Poisson distributed, with a
   mean that grows with the length of $G$'s alleles.
2. **Where each read came from.** Independently for each read: with some probability the read
   belongs to another locus and was placed here by mistake (it is *mismapped*), and then its bases
   say nothing about $G$. Otherwise it came from one of $G$'s haplotypes, chosen with a probability
   called that haplotype's *mixture weight*. Haplotypes that carry the same allele get equal
   weights; otherwise a haplotype's weight grows with the amount of sequence that only its allele
   contains.
3. **What the read says.** The part of the read inside the site is a copy of part of the allele
   that its haplotype carries, with sequencing errors. The probability of the read given that
   allele comes from its alignment score against the allele.

Some parts of the computation are heuristics rather than consequences of this process:

- The probability of a read given an allele is vg's alignment score converted to natural-log
  units. vg's scores are log-odds of the read's bases against a random-sequence background, so the
  conversion gives the log-likelihood only up to a term that depends on the read and not on the
  allele. Each read's likelihoods are therefore used relative to its best allele. The gap
  penalties do not come from an error model.
- The read's score against each allele is derived from the mapper's alignment of the read, as
  described under [Read log-likelihood](#read-log-likelihood), not from aligning the read to each
  allele again.
- A mismapped read is given the likelihood of its best allele at the site, which is the same under
  every genotype. So a read that fits every allele badly is no more likely to be mismapped than one
  that fits well.
- The probability that a read is mismapped is its MAPQ converted to a probability, held between a
  floor and a ceiling. The floor also stands for a read placed at the right locus whose alignment
  through this particular site is wrong, which MAPQ does not measure.
- The read count in the Poisson part counts each read by its probability of not being mismapped,
  so it is not a whole number.
- The logarithm of the Poisson part is multiplied by a weight below 1, because read depth varies
  between sites for reasons the model leaves out, such as mappability and base composition.
- Reads are treated as independent, although paired mates, and reads that share a systematic
  error, are not.

## Notation

| Symbol | Meaning |
|---|---|
| $A$ | the site's candidate alleles (see [Candidate alleles](#candidate-alleles)); $a, b \in A$ |
| $P$ | the ploidy at the site |
| $G = \lbrace g_1, \dots, g_P \rbrace$ | a genotype; $g_1, \dots, g_P$ list its alleles in any order, and nothing below depends on that order |
| $R$ | the reads of the site; $r \in R$ |
| $\ell_{ra}$ | log-likelihood score of read $r$ given allele $a$, in nats |
| $p_{ra}$ | relative likelihood of read $r$ given allele $a$, in $[0, 1]$ |
| $e_r$ | probability that read $r$ is mismapped |
| $\epsilon_{\min}, \epsilon_{\max}$ | floor and ceiling on $e_r$ |
| $w_i(G)$ | mixture weight of the haplotype carrying $g_i$ |
| $U_i(G)$ | length of the sequence in $g_i$ that the other allele of $G$ lacks |
| $T_a$ | length of allele $a$, excluding the site's two boundary nodes |
| $\bar L$ | mean read length near the site |
| $\kappa$ | read starts per base per haplotype near the site |
| $N_{\mathrm{eff}}$ | $\sum_{r \in R} (1 - e_r)$, the expected number of reads that came from the site |
| $\mu_G$ | expected number of reads from the site under $G$ |
| $\beta$ | weight of the read-count term |

## The likelihood

$$
\ln \mathcal{L}(G) = \sum_{r \in R} \ln \left[ (1 - e_r) \sum_{i=1}^{P} w_i(G) p_{r g_i} + e_r \right] + \beta \ln \mathrm{Pois}\left(N_{\mathrm{eff}} ; \mu_G\right)
$$

The sum over reads is the *read term*. For each read it averages the read's relative likelihood
over the haplotype of $G$ that produced it, and adds the case that the read is mismapped, in which
it has relative likelihood 1, that of its best allele (see
[Relative likelihood](#relative-likelihood)). A homozygous genotype lists the same allele $P$
times. The second part is the *depth term*, where

$$
\mu_G = \kappa \sum_{i=1}^{P} \max\left(T_{g_i} + \bar L - 1, 1\right)
$$

$$
\ln \mathrm{Pois}(n ; \mu) = n \ln \mu - \mu - \ln \Gamma(n + 1)
$$

$T_{g_i} + \bar L - 1$ is the number of positions at which a read overlapping the interior of
$g_i$ can start, and the floor of 1 keeps $\mu_G$ above zero. The gamma function extends the
Poisson distribution to a count $n$ that is not a whole number.

$\mathcal{L}(G)$ is a likelihood only up to a constant factor that is the same for every genotype
of the site, because each read's relative likelihoods leave out a factor that depends on the read
alone. It can be used to compare genotypes at one site, not to compare sites.

## How each term is computed

The sections below cover the terms in the order in which they are computed.

### The reads

The reads of a site are chosen as described under [Vocabulary](#vocabulary). Reads with MAPQ below
`--read-min-mapq` are discarded as they are read in. The reads come from one of three sources:

- `--gam` or `--gaf-reads`: a GAM or GAF file, loaded into memory.
- `--gam` with `--gam-index`: a GAM file sorted and indexed by `vg gamsort -i`, from which reads
  are fetched as they are needed, one window of consecutive node IDs at a time.
- `--gaf-base`: a GAF-Base database of alignments, queried one window at a time by running the
  `gbz-base` program (`--gaf-base-binary`) against the graph given with `--gbz-base`, or else
  against the input graph.

`--read-window` sets the window size, in node IDs, for the two indexed sources. It changes which
reads are fetched together, and so the order in which a site sees its reads, but not which reads a
site uses.

A read whose alignment crosses the site against the direction of the alleles is
reverse-complemented before it is scored. The direction is decided by a vote: the read is reversed
if more of its visits to the site's nodes are in the opposite orientation to the alleles' visits
than in the same orientation.

### Read log-likelihood

A read's placement says which nodes it visits and, inside each node, which of its bases match. An
allele is also a sequence of node visits. We pair the read's visits to the site's nodes (boundary
nodes included) with the allele's visits, in order, and add up the scores of the pairs and of the
visits left unpaired. Here a *visit* means a node in one orientation. The rules are:

| Read visit paired with allele visit | Score |
|---|---|
| the same visit | the score of the read's own edits in that node |
| a visit the allele never makes, paired with a visit the read never makes | the bases of the two nodes compared one by one from their first base, over the shorter length, plus a gap for the difference in length |
| any other pair | not allowed |

- A run of consecutive read visits left unpaired is an *insertion*, scored as one gap as long as
  their bases. A visit that the read and the allele share can be left unpaired in this way.
- A run of allele visits skipped between two pairs is a *deletion*, scored as one gap as long as
  their bases.
- Allele visits before the read's first pair or after its last lie outside the read and score
  nothing. Before the first pair, each unpaired read visit is scored as a gap of its own.

So every read base inside the site is scored under every allele, and all alleles are scored over
the same read bases, the read's *scoring window*.

Scores use vg's alignment scoring, in which a better fit scores higher. Matches and mismatches
score from a substitution matrix adjusted by base quality, in which a mismatch at a low-quality
base costs less; a read without base qualities uses the plain matrix. A gap of length $k$ scores
$-(o + (k - 1) x)$, where $o$ is `--gap-open` and $x$ is `--gap-extend`. The mapper's edits say
which bases match, mismatch, or are inserted or deleted, and their scores are computed here with
these settings. The total score $s_{ra}$ of the pairing is converted to nats:

$$
\ell_{ra} = \lambda s_{ra} + \iota n_{ra}
$$

$\lambda$ is the scorer's *log base*, the constant that makes vg's substitution scores
natural-log likelihood ratios. $n_{ra}$ is the number of gaps in which the read has bases the
allele lacks, counting those inside the mapper's edits, and $\iota$ is `--insertion-nats`. A
positive $\iota$ raises $\ell_{ra}$ for each such gap, so that extra read bases count against an
allele less than missing ones.

### The two walks

Choosing the pairing is itself a search, done in one of two ways. We call them the greedy walk and
the optimal walk, after the way they step through the two lists of visits.

- The default *greedy walk* makes one pass along the read's visits. A read visit is paired with
  the next occurrence of the same visit in the allele, if there is one, and once a first pair has
  been made, any allele visits passed over become a deletion. Otherwise, if the allele's next
  unpaired visit occurs later in the read, the read visit is left unpaired as an insertion.
  Otherwise the two visits are paired and their bases compared, whether or not the rules allow
  that pair. Once the allele's visits are used up, the remaining read visits are insertions. Each
  unpaired read visit is scored as a gap of its own.
- With `--realign`, the *optimal walk* finds the highest-scoring pairing that the rules allow, by
  dynamic programming over the read's visits against the allele's. When the product of the two
  numbers of visits exceeds a fixed limit, it only considers allele visits within a fixed distance
  of a diagonal that follows the visits the read and the allele share.

So the optimal walk applies the rules exactly at all but the largest sites, and the greedy walk
approximately. Inside a node that the read and the allele share, both walks score the mapper's
edits; neither aligns bases again.

### Relative likelihood

$\ell_{ra}$ is known only up to a term that differs from read to read, for the reason given under
[The model](#the-model). We remove it by using

$$
p_{ra} = \exp\left(\ell_{ra} - \max_{b \in A} \ell_{rb}\right)
$$

which is 1 for the read's best allele at the site and smaller for alleles that explain the read
worse. The values for one site form a matrix with a row per read and a column per allele, and each
row's maximum is 1.

### Mismapping probability

$$
e_r = \min\left(\max\left(10^{-\mathrm{MAPQ}_r / 10}, \epsilon_{\min}\right), \epsilon_{\max}\right)
$$

$\epsilon_{\min}$ is `--mismap-min` and $\epsilon_{\max}$ is `--mismap-max`. Each read's
contribution to the read term lies between $\ln e_r$ and 0, so one read changes the difference
between two genotypes' log-likelihoods by at most $-\ln \epsilon_{\min}$. The floor therefore sets
how strongly a single read can count against an allele. The ceiling applies to reads with MAPQ 0
or close to it, whose unclamped $e_r$ is near 1, and keeps them contributing a little rather than
nothing. Such reads are used unless `--read-min-mapq` excludes them. `--no-mismap-term` sets every
$e_r$ to $\epsilon_{\min}$, in $N_{\mathrm{eff}}$ and $\kappa$ as well as in the read term.

### Mixture weights

A read's term depends on the weights only when the read fits the alleles of $G$ differently: if
$p_{r g_1} = p_{r g_2}$, the average over haplotypes is the same whatever the weights. A read can
fit $g_i$ better than the other allele only if it overlaps sequence that $g_i$ has and the other
allele lacks. Treating that sequence as one stretch of length $U_i(G)$, such a read can start at
$U_i(G) + \bar L - 1$ positions, so each haplotype is weighted by that number:

$$
w_i(G) = \frac{\max\left(U_i(G) + \bar L - 1, 1\right)}{\sum_{j=1}^{P} \max\left(U_j(G) + \bar L - 1, 1\right)}
$$

$U_i(G)$ is the total length of the nodes that $g_i$ visits and the other allele of $G$ does not,
each node counted once. It is 0 for both haplotypes of a homozygous genotype, which therefore get
equal weights, as do two alleles whose unique sequence is equally long, such as two alleles that
differ by one base. At ploidy 1, $w_1(G) = 1$. `--flat-mixture` uses $w_i(G) = 1/P$ instead, so
that the effect of the weighting can be measured.

### Depth term

The read term takes the reads as given. The depth term asks whether their number fits $G$: a
homozygous deletion, for example, shows mainly as reads that are missing.

- $T_a + \bar L - 1$ is the number of positions at which a read overlapping the interior of allele
  $a$ can start. An allele with no interior, which deletes the whole interior, still gets
  $\bar L - 1$ positions, those of the reads that span its junction.
- $\kappa$ is measured over a *rate window*: the block of a fixed number of consecutive node IDs
  that contains the site's lowest node ID. Graphs built for vg number their nodes roughly in order
  along the genome, so this block is a stretch of genome around the site. $\kappa$ is the number
  of reads whose alignment begins in the window, each counted as $1 - e_r$, divided by the total
  length of the window's nodes and by the site's ploidy. It is computed once per window and shared
  by the sites in it.
- $\bar L$ is the mean length of the reads that begin in the rate window, or of the site's own
  reads if none does. The mixture weights use the same $\bar L$.
- $N_{\mathrm{eff}}$ counts each read of the site as $1 - e_r$. With `--depth-count-raw`, each read
  counts as 1, in both $N_{\mathrm{eff}}$ and $\kappa$.
- $\beta$ is `--depth-term`; 0 turns the depth term off.

The ratio $N_{\mathrm{eff}} / \mu_G$ at the called genotype is reported as `DR`, whether or not the
depth term is on.

## Calling a genotype

### Candidate alleles

- When the graph is a GBZ with at least two panel haplotypes, the candidate alleles are by default
  the distinct walks that panel haplotypes take through the site (*haplotype enumeration*). `-z`
  asks for this explicitly, and `-g` gives the panel as a separate GBWT file. Haplotype
  enumeration only offers alleles that some panel haplotype carries.
- Otherwise, or with `--enumerate-support`, the candidates are the walks with the most read
  support, found by Yen's k-shortest-paths algorithm over the node and edge coverage in a
  `vg pack` file given with `-k` (*support enumeration*). A fixed maximum number of walks is kept,
  which is larger with `-T`.
- The reference allele is always a candidate.

`--max-snarl-edges` skips a site with more edges than its limit, and the sites nested inside it are
then genotyped on their own; if it has none, nothing is called there. The limit exists because
Yen's search is slow on very large sites, so by default it is lifted when the candidates come from
haplotype enumeration, which does no search.

### Direct call or linkage

Every genotype of $P$ candidate alleles is scored; there is no pruning. The likelihoods are then
used in one of two ways:

- **Directly**, under support enumeration or when the linkage model is off (`--linkage-weight 0`):
  the call is the genotype with the highest $\mathcal{L}(G)$. An exact tie goes to the homozygous
  reference genotype if it is one of the tied genotypes, and otherwise to the tied genotype that
  comes first in VCF genotype order.
- **As emissions of the linkage model**, a hidden Markov model over the panel, under haplotype
  enumeration: the call is the genotype with the highest posterior probability under the model
  described in the next section.

A site with no reads gets no genotype, and the sites nested in it are then genotyped on their own.
Either way, the genotype chosen is the site's *settled* genotype, from which its *record*, its line
in the VCF, is written.

## The linkage model

Alleles at nearby sites are inherited together, so the sample's alleles at neighbouring sites tend
to form combinations that panel haplotypes also carry. We use this by modelling each of the
sample's haplotypes as a mosaic of panel haplotypes (the Li–Stephens model), with
$\mathcal{L}(G)$ as the emission probability, following PanGenie (Ebler et al. 2022). We call this
model the *linkage model*.

### Linkage chains, strands and states

The model runs along a *linkage chain*: a sequence of sites in reference order. At the top level,
a linkage chain holds the sites of one contig, split wherever the ploidy changes. Below the top
level, it holds the sites of one child chain (see [Nested sites](#nested-sites)). We call each of
the sample's haplotypes a *strand*, numbered 0 and 1; the word refers to a haplotype of the sample,
not to a strand of DNA. At ploidy 2 the hidden state at a site is an ordered pair $(h_1, h_2)$
naming the panel haplotype that each strand copies there. The pair implies a genotype, made of the
alleles that $h_1$ and $h_2$ carry at the site, and the state's emission is that genotype's
$\mathcal{L}(G)$. At ploidy 1 the state is a single panel haplotype.

### The wildcard

The panel is extended by a *wildcard* haplotype, which can carry any candidate allele, so that
genotypes that no panel pair carries can still be called. A strand on the wildcard, or on a panel
haplotype that does not pass through the site, has an unknown allele. The emission of such a state
is $\mathcal{L}(G)$ averaged uniformly over the candidate alleles the unknown strand could carry (an
average of likelihoods, not of their logarithms), multiplied by a fixed *escape* probability for each
unknown strand. The wildcard carries only candidate alleles, so it cannot add alleles to $A$.

### Transitions

Between consecutive sites $d$ bases apart, each strand independently keeps copying the same
haplotype or, with probability $\rho(d)$, picks a new one uniformly from the panel and the
wildcard:

$$
\rho(d) = \left[ \rho_{\min} + (1 - \rho_{\min}) \left(1 - e^{-d/D}\right) \right]^{\omega}
$$

$d$ is the difference between the two sites' reference positions; below the top level it is
measured along the parent site's settled allele, and where no distance is known, $\rho = 1$. $D$ is
`--linkage-scale`, the distance over which linkage decays. $\omega$ is `--linkage-weight`, an
exponent on the switch probability: a larger $\omega$ makes switches rarer and linkage stronger,
and $\omega = 0$ turns the linkage model off (see [Direct call or linkage](#direct-call-or-linkage)).
$\rho_{\min}$ is a small fixed floor, so that a switch is never impossible.

### Posterior and allele-frequency prior

The forward–backward algorithm gives each state's posterior probability at each site. A
genotype's posterior collects the probability of the states that imply it:

- A state whose strands both copy panel haplotypes that pass through the site implies one
  genotype.
- A state with one unknown strand, the other carrying allele $k$, shares its probability among the
  genotypes $\lbrace k, b \rbrace$, in proportion to their likelihoods.
- A state with two unknown strands shares its probability among all genotypes, in proportion to
  their likelihoods.

Summing over states favours genotypes that many panel pairs spell, which acts as a prior from the
panel's allele frequencies. The exponent $F$, `--linkage-prior`, sets its strength. The probability
collected from states of the first kind is multiplied by $c_G^{F-1}$, where $c_G$ is the number of
ordered panel pairs that spell $G$. That from states of the second kind is multiplied by
$n_k^{F-1}$, where $n_k$ is the number of panel haplotypes that carry $k$. States of the third kind
are not rescaled. The posteriors are then normalised to sum to 1. $F = 1$ keeps the prior that the
states imply, $F = 0$ removes it, and $F > 1$ strengthens it. At ploidy 1 the states are single
haplotypes, and the probability of those that carry allele $a$ is multiplied by $n_a^{F-1}$. The
called genotype is the one with the highest posterior.

`--hp-prior`, when not 0, replaces $F$ at sites where the reference allele and another candidate
differ only in the length of one homopolymer run, and the run is at least `--hp-prior-run` bases
long. Sequencing errors in long homopolymer runs tend to recur in many reads at the same site, which
the read term counts as independent evidence; a larger exponent there keeps the panel's share of
the decision.

### Windows

The forward–backward algorithm runs over overlapping windows of a fixed number of sites, so that
windows can be processed in parallel, and the results in a fixed margin at each window's edges are
discarded. This approximates a single pass over the whole linkage chain.

## Nested sites

*Nested calling* genotypes each child chain in records of its own. It is on by default with
`--read-likelihood`, and `--nested` turns it on with the other genotypers. `--no-nested` turns it
off: each site is then genotyped against its full walks, variation inside nested sites is reported
in the enclosing site's alleles, and a nested site is genotyped on its own only when its parent
could not be genotyped.

### Generations

Sites form a tree: a site can contain chains of smaller sites, which can contain chains in turn. A
site's *generation* is its depth in this tree. Sites not inside any other site are generation 0,
the sites in chains directly inside them are generation 1, and so on.

### Symbolic alleles

Two alleles of a site that take the same route except inside a child chain differ only in that
chain. To report such a difference once, we compare each called allele with the reference allele
in *symbolic* form: the walk with each pass through a child chain replaced by one symbol for that
chain. A called allele whose symbolic form equals the reference allele's is written as the
reference allele at this site (*symbolic collapsing*), and the records of the child chain's sites
report the difference. Genotyping the child chains of a site is called *descent*.

### A child chain's ploidy

A child chain is present only on the parent's haplotypes whose alleles pass through it, so its
ploidy is the number of the parent's settled alleles that cross it. An allele that crosses it twice
counts once.

- **0**: the sample has no copy of the chain, and no record is written for it or for anything
  nested inside it.
- **1**: the chain is genotyped at ploidy 1. When the genotypes are phased, its `GT` is written
  `a|.` or `.|a`, and the position of the allele says which of the parent's strands carries the
  chain. Unphased, it is written as a single allele.
- **2**: the chain is genotyped at ploidy 2.

### Off-reference chains

A chain that no reference path passes through (an *off-reference chain*) has no position on the
reference. It is genotyped in two cases:

- With `--anchors-out`, unless `--no-off-ref-nesting` is given, because its sites have anchors
  even though they have no VCF records. Genotyping these chains can change the phase written for
  other records.
- When the reference paths include a gRef fragment. `vg paths --compute-gref` adds a *gRef cover*
  to a graph: a copy of each reference path under a name starting `gref_`, and *fragments*, paths
  named `gref_<reference>_<N>_alt` that run through sequence the reference does not cover. A chain
  that a selected fragment passes through gets records, with the fragment as their contig.

The environment variable `VG_CALL_NO_REF_NESTED` also turns this genotyping on, for testing.

### The sweep and the barrier

Under the linkage model, a site's settled genotype is decided only after the model has seen its
whole linkage chain, and a child chain's ploidy depends on its parent's settled genotype. The reads
are therefore read once, in a pass we call the *sweep*, in which every site is scored: its
$\mathcal{L}(G)$ is computed for every genotype (a child chain's at both ploidy 1 and ploidy 2),
and its record is held back rather than written. After the sweep, the *barrier* settles genotypes
one generation at a time. The linkage model decides generation 0. Each generation-1 chain then
takes its ploidy from its settled parent (a chain at ploidy 0 is dropped with everything under
it), and the model decides it, starting from the parent's settled state, the pair of panel
haplotypes that the parent's strands copy (see [From the panel](#from-the-panel)). The same
happens for each later generation. After the barrier, each record is *rendered*, that is built,
from its settled genotype. Under direct calling a genotype is settled as soon as it is scored, and
the barrier only renders the records. [Re-genotyping](#re-genotyping-from-the-phase) runs the
barrier again.

### Several variants in one site

A called allele can differ from the reference allele in several places, separated by matching
sequence. `--atomize-blocks`, on by default with nested calling, writes one record per difference;
`--no-atomize-blocks` turns it off. We align the symbolic form of each called allele to that of the
reference allele by edit distance, and each maximal run of unmatched steps, a *block*, becomes a
record. A block record's `GT` gives the allele that each strand carries over that block. The blocks
come from symbolic forms, so block emission needs nested calling.

A site's block records share the site's ID (the VCF `ID` column) and the site's `AD`, `GL`, `GQ`,
`GQI`, `GQN`, `GP`, `QUAL`, `DP`, `DR` and `BL`, because the likelihood is computed for the whole
site. Summing or averaging these fields over a site's records therefore counts the site's evidence
more than once. `INFO/SB` gives each record's index and the number of records the site wrote.

## Phasing

### From the panel

The linkage model also finds the most probable sequence of states through each linkage chain (the
Viterbi path), over the same windows, each window continuing from the state the previous window
chose. At every site the path is restricted to states that imply the settled genotype, where a
state with an unknown strand qualifies if that strand could carry the allele it needs. So phasing
orders a genotype without changing it. The order is the phase: strand 0's allele is written to the
left of the `|` in `GT`. `FORMAT/PS` names the *phase set*: the reference position of the first
site of the top-level linkage chain, which nested sites share with their parent.

Phasing is on wherever the linkage model runs, and `--phased` makes a run fail when the linkage
model cannot run. `--no-phased` turns phasing off. Where the linkage model runs, it also turns
nested calling off, because a child chain's ploidy and strand come from its parent's phased pair;
variation inside nested sites is then reported in the enclosing site's alleles.

### From the reads

A read that spans two heterozygous sites shows directly whether their alleles lie on the same
strand. `--read-phasing` uses such reads to re-decide the order of each heterozygous site's pair,
within the phase sets the panel gave. It works on the diploid heterozygous sites of each phase set,
nested ones included, in reference order, and it changes no genotype. Reads are identified by name,
so paired mates count as one read.

#### What each read says

For a read $r$ at a heterozygous site $s$ with settled pair $(a_1, a_2)$, let
$x_{rsk} = (1 - e_r) v_k p_{r a_k}$ for $k \in \lbrace 1, 2 \rbrace$. The *allele-length weights*
$v_1, v_2$ are computed like the mixture weights, but from each allele's full length (all its
nodes, boundary nodes included) in place of $U$, so that each allele's weight depends on its own
length alone. A read with $x_{rs1} + x_{rs2} = 0$ is not used at $s$. For each read we keep:

- $q_{rs} = x_{rs1} / (x_{rs1} + x_{rs2})$, the probability that the read carries $a_1$, given
  that it came from one of the two strands;
- $c_{rs} = (x_{rs1} + x_{rs2}) / (x_{rs1} + x_{rs2} + e_r)$, the probability that it did come from
  one of them;
- its *confidence*, $-10 \log_{10}\left(1 - \max(x_{rs1}, x_{rs2}) / (x_{rs1} + x_{rs2} + e_r)\right)$,
  the phred-scaled probability that its better allele is wrong (the anchor file calls it the read's
  *score*). At a heterozygous site it can be at most
  $-10 \log_{10}\left(\epsilon_{\min} / (\epsilon_{\min} + (1 - \epsilon_{\min})/2)\right)$,
  the *heterozygous score ceiling*.

A site is *reliable* if the mean confidence of its reads is at least `--phase-min-q`. A
`--phase-min-q` above the heterozygous score ceiling is rejected, since no site could be reliable.

#### Links between sites

Two sites $s$ and $t$ are *linked* by the reads they share. For a shared read, let
$m_r = q_{rs} q_{rt} + (1 - q_{rs})(1 - q_{rt})$ be the probability that its alleles at the two
sites lie on one strand in the current orders, and assume that the read reports the truth with
probability $\gamma_r = c_{rs} c_{rt}$ and is a coin flip otherwise. The link is

$$
\mathrm{link}(s, t) = \sum_{r} \log_{10} \frac{\gamma_r m_r + (1 - \gamma_r)/2}{\gamma_r (1 - m_r) + (1 - \gamma_r)/2}
$$

A positive link agrees with the two sites' current orders. `--phase-cap`, when not 0, limits the
size of each link.

#### Deciding the orders

The orders are decided in four steps:

1. **Chain.** The reliable sites of a phase set, in order, form the *phase chain*, and each is
   linked to the next. The chain breaks wherever the size of a link is below `--phase-break`
   $\log_{10}$ units. Within each unbroken piece, a site is flipped relative to the site before it
   when their link is negative, and the first site of the piece keeps the panel's order.
2. **Relink.** Each break is decided from the links between the last `--phase-relink` chain sites
   before it and the first `--phase-relink` after it, taken relative to the two pieces' current
   orders. If their sum says that the pieces disagree, every site of the later piece is flipped. If
   no read links them, the panel's order stands.
3. **Coherence.** For each read that spans two or more chain sites, and each such site $s$, we ask
   which strand the read's other chain sites put it on, and whether its allele at $s$ agrees. The
   strand is the one with the higher sum over the other sites $t$ of
   $\log_{10}(c_{rt} q_{rt} + (1 - c_{rt})/2)$ for strand 0, or of
   $\log_{10}(c_{rt}(1 - q_{rt}) + (1 - c_{rt})/2)$ for strand 1, with each $q$ taken in the
   current order. A site's *coherence* is the fraction of its reads that agree. Sites with enough
   such reads and coherence below `--phase-coherence` are removed from the chain, and steps 1 and 2
   are repeated on the sites left, up to `--phase-coh-rounds` times. `--phase-coherence 0` skips
   this step.
4. **Hang.** Every other heterozygous site of the phase set is oriented from its links to its
   nearest chain sites, $\lfloor H/2 \rfloor + 1$ on each side for $H$ = `--phase-hang`, each link
   counted with its sign and size. A vote of weight `--phase-prior`, in the same $\log_{10}$ units,
   is added for keeping the panel's order relative to the nearest chain site.

When a parent site's order changes, the strand of its nested haploid chains changes with it.

### Re-genotyping from the phase

The mixture weights are the same for every read, because one site on its own cannot tell which
strand a read came from. Once sites are phased, a read that spans other heterozygous sites can.
`--regenotype`, which needs `--read-phasing`, uses this to correct each site's likelihoods. The
quantities $q_{rt}$, $c_{rt}$ and $v_k$ are those of [From the reads](#from-the-reads), and
$\sigma$ is the logistic function.

#### Strand log-odds

A read's *strand log-odds* at site $s$ is the natural-log odds that it comes from strand 0, summed
over the other phased heterozygous sites $t$ that it overlaps in the same phase set, with each
$q_{rt}$ taken in the order that read phasing settled:

$$
\Lambda_{rs} = \sum_{t \neq s} \ln \frac{c_{rt} q_{rt} + (1 - c_{rt})/2}{c_{rt}(1 - q_{rt}) + (1 - c_{rt})/2}
$$

Leaving out $s$ keeps a site from confirming its own genotype. A read seen in more than one phase
set has no usable strand, and its $\Lambda_{rs}$ is 0.

#### Tempering

Reads are not independent, so $\Lambda_{rs}$ overstates how sure the strand is. The *tempered*
strand log-odds is

$$
y_{rs} = \mathrm{logit}\left(C \sigma(\tau \Lambda_{rs}) + (1 - C)/2\right)
$$

where $\tau$ is the *temper* (`--regeno-temper`) and $C$ is `--regeno-ceiling`. With $C = 1$,
$y_{rs} = \tau \Lambda_{rs}$, and a smaller $C$ keeps the probability of either strand further
from 1.

Unless it is given, $\tau$ is fitted once, from the phased sites of the first round. Each read at
each site $s$ with $\Lambda_{rs} \neq 0$ is an observation: does the sign of $\Lambda_{rs}$ agree
with the strand that the read's own allele at $s$ points to? The observations are grouped into bins
of equal size by $\vert \Lambda_{rs} \vert$, and $\tau$ is chosen from a fixed grid to minimise the
squared difference, weighted by bin size, between $C \sigma(\tau \vert \Lambda \vert) + (1 - C)/2$
at each bin's mean and the bin's rate of agreement. With too few observations, $\tau$ is 0 and
re-genotyping changes nothing.

#### The correction

At a heterozygous genotype with alleles $a$ and $b$, each read gets its own weights

$$
\pi_{ra} = \frac{v_a e^{y_{rs}}}{v_a e^{y_{rs}} + v_b}, \qquad \pi_{rb} = 1 - \pi_{ra}
$$

with $a$ on strand 0. The correction added to $\ln \mathcal{L}(G)$ is

$$
\sum_{r \in R} \left[ \ln\left((1 - e_r)(\pi_{ra} p_{ra} + \pi_{rb} p_{rb}) + e_r\right) - \ln\left((1 - e_r)(v_a p_{ra} + v_b p_{rb}) + e_r\right) \right]
$$

The other assignment of the two alleles to the strands is scored too, and the larger correction is
kept. Homozygous genotypes are unchanged and have no assignment to choose, so this choice can only
favour heterozygous genotypes.

#### Nested haploid chains

A nested haploid chain has one strand, so there is nothing to reweight. There, a read whose
tempered strand log-odds point to the parent's other strand is made less informative: each of its
relative likelihoods $p$ becomes $\eta p + 1 - \eta$, where $\eta = \min(1, e^{y})$ and $y$ is its
tempered strand log-odds towards the chain's strand. $\eta$ is 1 when $\tau = 0$ or when the read
points to the chain's strand, so those reads are unchanged. `--no-regeno-haploid` turns this off.

#### Rounds

The corrected likelihoods go back to the linkage model, the barrier runs again, and read phasing
runs again on the new genotypes. New genotypes change the phase and so $\Lambda$, so the cycle can
repeat. Each round corrects the likelihoods from the sweep, not those of the previous round.
`--regeno-passes` caps the total number of barrier passes, the first included, and the cycle stops
sooner if no settled genotype changes or the genotypes return to an earlier state. With
`--regeno-passes 1` the correction is computed and reported but not applied. `--regeno-ledger`
writes one line for each site whose best genotype the correction changes.

## Output

### VCF fields

| Field | Meaning |
|---|---|
| `GT` | the settled genotype, phased where phasing ran |
| `GL` | $\log_{10} \mathcal{L}(G)$ for every genotype of the record's alleles, in VCF order |
| `GQ` | the difference between the log-likelihoods of the best and second-best genotypes, in phred units, multiplied by the *explained share*: the fraction of reads whose best allele is a called allele |
| `GQI` | the same difference, not multiplied by the explained share |
| `GQN` | the same difference as a fraction of the largest difference the site could give, multiplied by the explained share |
| `GP` | natural log of the posterior probability, from $\mathcal{L}(G)$ under a uniform prior over genotypes, of the genotype with the highest $\mathcal{L}(G)$ |
| `QUAL` | phred-scaled posterior probability, under the same uniform prior, of the homozygous reference genotype; 0 for a homozygous reference call |
| `DP` | number of reads of the site |
| `AD` | for each allele in the record, the number of reads whose best allele it is, rounded; a read tied between alleles counts a fraction to each |
| `DR` | $N_{\mathrm{eff}} / \mu_G$ at the called genotype |
| `BL` | mean over reads of $\max_a \ell_{ra}$, each read's best log-likelihood score at the site |
| `FORMAT/PS` | the phase set |
| `INFO/SB` | the record's index and the number of records its site wrote under `--atomize-blocks` |
| `FILTER=noreads` | the site had no reads |
| `FILTER=lowconf` | `GQN` is below `--min-confidence` |

A read whose best allele is not in the call usually fits the called genotype and its runner-up
about equally, so the difference between them hardly reflects it. Multiplying by the explained
share lowers `GQ` when the call leaves reads unexplained, which makes `GQ` a score for ranking
calls rather than a posterior. `GL`, `GQ` and `GP` are also over-confident at high depth, because
reads are treated as independent. `AD` need not sum to `DP`: every candidate allele was scored,
but only the alleles written in the record have an entry.

`GQ` depends on depth, since the likelihood difference is a sum over reads. It also depends on
ploidy: at ploidy 1 the runner-up is a different allele, which most reads can tell apart from the
call, while a diploid heterozygote's runner-up differs from it on only one strand, which fewer
reads can. `GQN` removes both effects. Its denominator, the *achievable gap*, is the difference the
read term alone would give between the called genotype and the runner-up if every read came from
one of the called genotype's haplotypes, in proportion to its mixture weights, fitted that
haplotype's allele with relative likelihood 1 and every other allele with 0, and had
$e_r = \epsilon_{\min}$. `GQN` is `.` when there is no difference to normalise (no reads, or a
single possible genotype), and otherwise lies in $[0, 1]$, except on the records described next.

When the linkage model settles a genotype other than the one with the highest $\mathcal{L}(G)$,
`GQ` becomes $-10 \log_{10}(1 - \text{posterior})$ times the explained share, capped at `GQI`. The
posterior includes the panel's allele-frequency prior, and the cap keeps the reported confidence
from exceeding what the reads alone support. `GQN` becomes the called genotype's margin over the
best other genotype in `GL`, as a fraction of the achievable gap and multiplied by the explained
share. It lies in $[-1, 1]$ and is negative when the linkage model called against the reads.
`lowconf` is then decided from this `GQN`. Records whose genotype the linkage model did not change
keep the per-site `GQ`.

When re-genotyping is applied (`--regeno-passes` above 1), `GL` holds the corrected likelihoods,
and at a site whose best genotype the correction changed, `GQ` is recomputed from them as the phred difference between the two best
genotypes. `GQN` is not changed.

`--no-share-quality` writes the unmultiplied value as `GQ`. `--depth-quality` multiplies `GQ` by
$e^{-D_q \vert \ln \mathrm{DR} \vert}$, where $D_q$ is its value, at records whose called alleles
change length by at least a fixed amount. `--min-confidence` marks records with `FILTER=lowconf`
and does not remove them. None of these three changes a genotype.

With `-A`, `--top-down` or `--bottom-up`, and whenever off-reference chains are genotyped, records
carry vg's nesting INFO tags. `INFO/CH` counts the levels of non-reference sequence that
separate a record from the linear reference: a gRef fragment that hangs off the reference is level
1, a fragment that hangs off a level-1 fragment is level 2, and so on. `INFO/LV`, `INFO/PS`,
`INFO/RC`, `INFO/RS` and `INFO/RD` are described in the VCF header.

### Mosaic (`--mosaic-out`)

The mosaic file describes each of the sample's strands as a walk through the graph, stating which
panel haplotype it copies where. It is the phasing in another form, so it needs the linkage model,
turns phasing on, and cannot be combined with `--no-phased`.

The file is tab-separated. Header lines start with `#` and are identified by their first field:

- `#mosaic-version`, the format version;
- `#graph`, the input graph;
- `#sample`, the sample name (`-s`);
- `#reference`, one line for each reference path that positions refer to;
- `#gref-fragments`, the number of gRef fragments among the reference paths, when there are any;
- `#decoding`, how the strands were chosen, `constrained-viterbi`: the Viterbi path restricted to
  the called genotypes, as in [From the panel](#from-the-panel);
- `#patch`, `#nested` and `#unexplained`, which record the choices of the three options described
  below;
- `#haplotype`, one line per panel haplotype, giving its index and name;
- `#note`, text describing the columns;
- `#H`, the column names.

A reader should skip header keys it does not recognise.

A data line starts with `H` and describes a *run*: a maximal stretch of consecutive sites, on one
strand, over which the strand copies one panel haplotype. A strand's runs join into one walk, but
the walk has to start again where a run cannot continue the previous one: where a gap is left
unfilled, or where the direction of travel reverses, as at an inversion. Each such walk is a
*fragment*.

| Column | Meaning |
|---|---|
| `contig` | reference contig of the run's sites |
| `strand` | `0` or `1`; strand 0 carries the allele to the left of the `\|` in `GT` |
| `fragment` | with `contig` and `strand`, identifies one walk |
| `ref_start`, `ref_end` | reference positions of the run's first and last sites, in the coordinates of the `#reference` paths; for orientation only, since the nodes define the run |
| `start_node`, `end_node` | oriented node IDs (node ID times 2, plus 1 if reverse) where the run starts and ends |
| `hap_index` | the panel haplotype's row in the `#haplotype` table; `ref` on a row filled with the reference, where no panel haplotype can be followed; `*` on a row where the strand is on the wildcard |
| `haplotype` | the panel haplotype's name, as its sample name and haplotype number joined by `#`; for a `ref` row, the reference's name; `*` on a wildcard row |
| `sites` | number of called sites in the run, or `.` for a `ref` row that fills a gap between two runs |
| `gbwt_offset` | with `start_node`, a position in the graph's GBWT from which the haplotype can be followed to `end_node`; `.` if there is none |

Consecutive rows of a fragment meet at a shared node, so a fragment expands to one walk in the
graph. A run never crosses the end of one of the paths that store its panel haplotype, so a
haplotype stored as several paths gives several rows. `gbwt_offset` is valid only for the graph
named in `#graph`, and `hap_index` only within one file; `haplotype` is the name to compare across
files. A haploid contig or region has rows for strand 0 only.

A strand's sites include the nested sites on its walk, so a change of panel haplotype at a nested
site starts a new row. Three options change how the runs are formed, and the header records each
choice:

- Between two runs, the panel haplotype of the first is carried on to the second where it can
  be. Where it cannot, the gap is filled with the reference, on a `ref` row, and so is a run whose
  own haplotype the graph does not carry across it. `--no-mosaic-patch-gaps` leaves such gaps
  unfilled, starting a new fragment after each.
- `--no-mosaic-nested` leaves nested sites out of the runs, so a strand follows its enclosing
  site's haplotype through them.
- Where a strand is on the wildcard, the panel cannot name a haplotype for it. By default those
  sites are left out of the runs, so the neighbouring haplotype is carried through them.
  `--mosaic-break-unexplained` writes a `*` row for them instead.

### Assembly anchors (`--anchors-out`)

The anchor file records, for each genotyped site, which reads support which of the sample's
strands, for use in pangenome-guided assembly. An *anchor* is a point between two adjacent bases
of the graph, with no sequence of its own, at which reads are tied to the site. Each site has two:
the *start
anchor*, just after its start boundary node, where the walk enters the interior, and the *end
anchor*, just before its end boundary node, where it leaves.

A site's reads are divided among its *slots*, one for each distinct called allele. At a
heterozygous site, slot 0 holds the reads of the allele to the left of the `|` in `GT` and slot 1
those of the allele to the right. A homozygous site has a single slot, as does a nested haploid
chain, whose slot is the strand that carries it. Each read goes to the slot of the called allele
with the larger $x_{rk} = (1 - e_r) v_k p_{r a_k}$, with the allele-length weights $v_k$ of
[From the reads](#from-the-reads). Its *confidence* in that slot is
$-10 \log_{10}\left(1 - \max_k x_{rk} / (\sum_k x_{rk} + e_r)\right)$, as in read phasing, and a
site's *reliability* is the mean confidence of the reads written for it.

Each anchor row (`A`) gives the boundary node, the site's ID (the VCF `ID` column, also written for
off-reference sites, which have no VCF record), the slot, the allele as an index into the site's
candidate alleles (not the VCF allele number), and three site-level values: `GQN`, the explained
share, and the reliability. The read rows (`R`) that follow give each read's index in the header's
table of read names, the direction in which it crosses the anchor (0 in the site's direction, 1
against it), the 0-based offset in the read as sequenced of the base next to the anchor, and the
read's confidence. The file's header describes every column.

Options select sites and reads:

- `--anchors-het-only` writes anchors only at heterozygous sites, and `--anchors-leaf-only` only at
  sites with no child chains.
- `--anchors-min-gqn` skips sites whose `GQN` is below its value.
- `--anchors-min-q` skips reads whose confidence is below its value.
- `--anchors-keep-off-call` keeps reads whose best candidate allele was not called; by default they
  are left out.
- `--anchors-reads` skips anchors with fewer reads than its value.
- `--anchors-end-new` writes a slot's end anchor only if it has at least that many reads that the
  slot's start anchor lacks.

With `--read-phasing`, the tempered strand log-odds of [Re-genotyping](#re-genotyping-from-the-phase)
are also available, with $\tau$ fitted as described there when `--regenotype` is off:

- At a heterozygous site, a read's slot is chosen from its $x_{rk}$ with $v_k$ replaced by the
  per-read weights $\pi$ of [The correction](#the-correction). `--no-anchors-phase-hets` uses
  the $x_{rk}$ alone, and `--anchors-strict-hets` uses only the sign of the strand log-odds.
- `--anchors-hom-split` divides a homozygous site's reads between two slots, carrying the same
  allele, by the sign of their tempered strand log-odds. A site is split only when at least
  `--split-min-side` reads on each side have a tempered strand log-odds of size at least
  `--split-min-q`.

### Ploidy

A site's ploidy comes from `-d`. `-R` overrides it per contig with a comma-separated list of
`REGEX:PLOIDY` rules, the first rule whose regular expression matches a reference path's name
applying to that path. `--ploidy-bed` overrides both per region, with a BED file of
`CHROM START END PLOIDY` lines: `CHROM` is spelled as in the output VCF, intervals are 0-based and
half-open and must not overlap, and a site takes the ploidy of the interval containing its start.
Every ploidy must be 1 or 2. A nested chain takes its ploidy from its parent, as described
[above](#a-child-chains-ploidy). A linkage chain ends wherever the ploidy changes. At ploidy 1,
`GT` is a single allele.

## Options

`vg call --help` gives each option's default. `--preset ont` sets several of the options below to
values suited to Oxford Nanopore reads, and `vg call --help` lists which. An option given
explicitly overrides the preset. Where an option has a `--no-` form, the two set the same thing and
the one given last wins. The read-likelihood options are rejected without `--read-likelihood`, and
the options that modify `--anchors-out`, `--mosaic-out` and `--regenotype` are rejected when those
are not in use.

| Part | Options |
|---|---|
| Reads | `--gam`, `--gaf-reads`, `--gam-index`, `--gaf-base`, `--gbz-base`, `--gaf-base-binary`, `--read-window`, `--read-min-mapq` |
| Candidate alleles | `--enumerate-support`, `-k`, `-g`, `-z`, `--max-snarl-edges` |
| Read log-likelihood | `--gap-open`, `--gap-extend`, `--insertion-nats`, `--realign`, `--no-realign` |
| Mismapping | `--mismap-min`, `--mismap-max`, `--no-mismap-term` |
| Mixture weights | `--flat-mixture` |
| Depth term | `--depth-term`, `--depth-count-raw` |
| Linkage | `--linkage-weight`, `--linkage-scale`, `--linkage-prior`, `--hp-prior`, `--hp-prior-run` |
| Nested sites | `--nested`, `--no-nested`, `--atomize-blocks`, `--no-atomize-blocks`, `--no-off-ref-nesting` |
| Ploidy | `-d`, `-R`, `--ploidy-bed` |
| Phasing | `--phased`, `--no-phased`, `--read-phasing`, `--no-read-phasing`, `--phase-min-q`, `--phase-coherence`, `--phase-coh-rounds`, `--phase-break`, `--phase-relink`, `--phase-hang`, `--phase-prior`, `--phase-cap` |
| Re-genotyping | `--regenotype`, `--no-regenotype`, `--regeno-temper`, `--regeno-ceiling`, `--regeno-passes`, `--regeno-haploid`, `--no-regeno-haploid`, `--regeno-ledger` |
| Quality fields | `--no-share-quality`, `--depth-quality`, `--min-confidence` |
| Mosaic | `--mosaic-out`, `--mosaic-patch-gaps`, `--no-mosaic-patch-gaps`, `--no-mosaic-nested`, `--mosaic-break-unexplained` |
| Anchors | `--anchors-out`, `--anchors-reads`, `--anchors-min-gqn`, `--anchors-min-q`, `--anchors-het-only`, `--anchors-leaf-only`, `--anchors-keep-off-call`, `--anchors-end-new`, `--anchors-phase-hets`, `--no-anchors-phase-hets`, `--anchors-strict-hets`, `--anchors-hom-split`, `--split-min-q`, `--split-min-side` |
| Debugging and evaluation | `--dump-likelihoods` (writes each read's $e_r$ and $p_{ra}$ at every site as TSV), `--regeno-shuffle` (randomises the sign of each $\Lambda_{rs}$ before re-genotyping, as a control for how much the phase contributes); `--flat-mixture` and `--anchors-strict-hets`, listed above, also serve to measure the parts they replace |
| Presets | `--preset` |

### Fixed constants

These constants have no option. Each is named as it appears in the source.

| Constant | Defined in | What it sets |
|---|---|---|
| `LinkageModel::Params::escape` | `src/linkage_model.hpp` | the escape probability for a strand with an unknown allele |
| `LinkageModel::Params::rho_min` | `src/linkage_model.hpp` | the floor $\rho_{\min}$ on the switch probability |
| `LinkageModel::Params::window`, `margin` | `src/linkage_model.hpp` | sites per forward–backward window, and sites discarded at each edge |
| `RATE_WINDOW`, in `local_read_stats` | `src/allele_likelihood.cpp` | node IDs per rate window, for $\kappa$ and $\bar L$ |
| `banded` and `band`, in `score_read_against_allele` | `src/allele_likelihood.cpp` | when the optimal walk is restricted, and the half-width of the restriction in allele visits |
| `max_yens_traversals` | `src/subcommand/call_main.cpp` | the most candidate alleles support enumeration keeps, with and without `-T` |
| the `min_length` argument of `set_depth_quality` | `src/read_likelihood_caller.hpp` | the change in allele length at which `--depth-quality` applies |
| the minimum counted reads in `read_phase_flips` | `src/read_phasing.cpp` | the reads a site needs before low coherence can remove it from the phase chain |
| `RegenotypeParams::fit_bins`, `fit_min_per_bin`, and the grid in `fit_calibration` | `src/regenotype.hpp`, `src/regenotype.cpp` | the bins, the fewest observations per bin, and the candidate values used to fit $\tau$ |

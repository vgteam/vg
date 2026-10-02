# Read phasing and re-genotyping

`vg call --read-likelihood` phases each linkage chain from the panel, as the linkage model's Viterbi
path (see [Phasing from the panel](read-likelihood-linkage-model.md#phasing-from-the-panel)). Two
optional steps then use the reads' evidence again, with statistics of their own rather than the
linkage model's:

- **Read phasing** (`--read-phasing`) re-decides the phases from reads that span several sites.
- **Re-genotyping** (`--regenotype`) uses the phase to correct each site's likelihoods, so that
  the next linkage pass chooses the genotypes from the corrected likelihoods.

Both work from the evidence that direct genotyping keeps for each read of a site: its mismapping
probability $e_r$ and its relative likelihood $p_{ra}$ under each candidate allele $a$ (see
[Direct genotype](read-likelihood-genotyping.md#direct-genotype)). Read phasing ends every round,
and re-genotyping starts every round after the first (see
[Passes and rounds](read-likelihood-genotyping.md#passes-and-rounds)). The caller as a whole, and
the terms used here, such as phaseable site, phase set, strand and nested haploid chain, are
described in [read-likelihood-genotyping.md](read-likelihood-genotyping.md); how the options turn
these steps on is under [Phasing](read-likelihood-genotyping.md#phasing) there.

## Read phasing

A read that spans two phaseable sites shows directly whether their alleles lie on the same strand.
Read phasing uses such reads to re-decide the phase of each phaseable site, within the phase sets
the panel gave. It changes phases and keeps every genotype.

Read phasing takes one phase set at a time, and its phaseable sites in order of position, with the
sites of an off-reference chain placed as the linkage model places them (see
[Transitions](read-likelihood-linkage-model.md#transitions)).

### What each read says

Take a read $r$ of a phaseable site $s$, and let $a_0$ and $a_1$ be the site's alleles on strand 0
and strand 1 in its current phase. For $k \in \lbrace 0, 1 \rbrace$, let $x_{rsk} = (1 - e_r)
v_{a_k} p_{r a_k}$, where $e_r$ is the read's mismapping probability and $p_{r a_k}$ its relative
likelihood under $a_k$ at $s$ (see [Direct
genotype](read-likelihood-genotyping.md#direct-genotype)). The **allele-length weights** are

$$
v_{a_k} = \frac{\lvert a_k \rvert + \bar L - 1}{\left(\lvert a_0 \rvert + \bar L - 1\right) +
\left(\lvert a_1 \rvert + \bar L - 1\right)}
$$

where $\lvert a \rvert$ is the total sequence length of allele $a$'s walk, boundary nodes included,
and $\bar L$ is the mean read length of
[Depth term inputs](read-likelihood-direct-genotyping.md#depth-term-inputs). A read overlapping
$a$ can start at $\lvert a \rvert + \bar L - 1$ positions, so the weights estimate the share of the
site's reads that come from each strand. `--flat-mixture` sets both to $1/2$. A read with
$x_{rs0} + x_{rs1} = 0$ is not used at $s$. For each read used, vg keeps:

- $q_{rs} = x_{rs0} / (x_{rs0} + x_{rs1})$, the probability that the read carries $a_0$, given
  that it came from one of the two strands;
- $c_{rs} = (x_{rs0} + x_{rs1}) / (x_{rs0} + x_{rs1} + e_r)$, the probability that it did come from
  one of them, rather than being mismapped, with a mismapped read taken to fit with relative
  likelihood 1, as in the read term;
- its **confidence**,
  $-10 \log_{10}\left(1 - \max(x_{rs0}, x_{rs1}) / (x_{rs0} + x_{rs1} + e_r)\right)$, the
  phred-scaled probability that its better allele is wrong.

A read is identified by its name, so the two mates of a pair are one read. Where both mates are
used at one site, they come from one molecule and often read the same bases, so the site keeps
only the mate with the higher confidence (on a tie, the larger $q_{rs}$, then the larger $c_{rs}$),
and everything below counts each read once.

### Links between sites

Two sites $s$ and $t$ are **linked** by the reads they share. For a shared read $r$, let

$$
m_r = q_{rs} q_{rt} + (1 - q_{rs})(1 - q_{rt})
$$

This is the probability that the read carries the same strand's allele at both sites, under the
sites' current phases. The read is assumed to report this relation truly with probability
$\gamma_r = c_{rs} c_{rt}$, and to be a coin flip otherwise. The link is

$$
\mathrm{link}(s, t) = \sum_{r} \log_{10} \frac{\gamma_r m_r + (1 - \gamma_r)/2}{\gamma_r (1 - m_r) + (1 - \gamma_r)/2}
$$

A positive link favours the two sites' current phases. A link's **size** is its absolute value.
`--phase-cap`, when not 0, limits the size of each link.

### Reliable sites

A site's **reliability** is the mean confidence of its reads, and the site is **reliable** if this
is at least `--phase-min-q`. Consider a read with $e_r = \epsilon_{\min}$, the floor on $e_r$ set
by `--mismap-min`, at a site whose two alleles have equal length. If the read fits one allele
perfectly and the other not at all, its confidence is
$-10 \log_{10}\left(\epsilon_{\min} / (\epsilon_{\min} + (1 - \epsilon_{\min})/2)\right)$, the
**heterozygous score ceiling**. A read that favours the longer of two unequal alleles can have a
higher confidence. Under read phasing, vg rejects a `--phase-min-q` above the ceiling, since sites
whose alleles are of similar length could not reach it.

### Deciding the phases

To **flip** a site is to reverse its phase. For each site, read phasing decides whether to flip it
against the phase the panel gave. It phases the reliable sites first, each relative to the one
before it, and then phases every other site on its own. A wrong link between reliable sites flips
every site after it, while a wrong decision for another site flips only that site. Each phase set
goes through four stages:

1. **Phase chain.** The reliable sites of the phase set, in order of position, form its **phase
   chain**, and each is linked to the next. The phase chain breaks wherever the size of a link is
   below `--phase-break` $\log_{10}$ units, and the breaks divide it into **pieces**. The first site
   of each piece keeps the panel's phase. Each later site of the piece is flipped relative to the
   site before it when their link, taken in the panel's phases, is negative.
2. **Relink.** At each break, read phasing decides from the links across it whether to flip the
   later piece. Take the last `--phase-relink` sites of the piece before the break and the first
   `--phase-relink` sites of the piece after it, or all of a piece's sites if it has fewer. Sum the
   links between every site on one side and every site on the other, each in the current phases. If
   the sum is negative, every site of the later piece is flipped. Otherwise, as when no read links
   the two sides, the later piece keeps the phases stage 1 gave it. Breaks are decided from left to
   right, so each piece is oriented against the piece before it as that piece now stands.
3. **Coherence.** This stage removes phase-chain sites whose reads disagree with the rest of the
   phase chain. Take a read that spans two or more phase-chain sites, and one of those sites, $s$.
   The read's other phase-chain sites $t$ vote for the strand it came from. Strand 0 scores $\sum_t
   \log_{10}(c_{rt} q_{rt} + (1 - c_{rt})/2)$, strand 1 scores $\sum_t \log_{10}(c_{rt}(1 - q_{rt})
   + (1 - c_{rt})/2)$, and the higher score wins. The read **agrees** at $s$ if $q_{rs}$ lies on the
   winning strand's side of $1/2$. Every $q$ here is taken in the current phase. A site's
   **coherence** is the fraction of its reads that agree, counting only reads that span another
   phase-chain site. A site with at least a fixed number of counted reads (see [Fixed
   constants](read-likelihood-genotyping.md#fixed-constants)) and a coherence below
   `--phase-coherence` is removed from the phase chain. If any site is removed, stages 1 and 2 run
   again on the sites left, starting again from the panel's phases, and then this stage runs again.
   This stops when the stage removes no site, or after it has removed sites `--phase-coh-rounds`
   times (once if that is 0). The stage runs only on a phase chain of at least 3 sites, and removes
   sites only when at least 2 would be left. `--phase-coherence 0` skips this stage.
4. **Hang.** Each phaseable site of the phase set that is not in the phase chain, because it is
   unreliable or stage 3 removed it, is then phased on its own. Its nearest phase-chain sites are
   used, up to $\lfloor H/2 \rfloor + 1$ on each side, so that $H$, `--phase-hang`, is shared
   between the two sides. The link to each of them votes for the phase it implies, with a weight
   equal to its size. A further vote of weight `--phase-prior`, in the same $\log_{10}$ units,
   keeps the site's phase relative to the nearest phase-chain site as the panel gave it: it
   favours flipping the site if that site was flipped, and keeping its phase if not. The site
   takes the phase with the larger total vote, and a tie keeps the panel's phase.

When read phasing flips a site, its alleles change strands, and so does every nested haploid chain
that takes its strand from them, directly or through another nested haploid chain. A diploid
nested site keeps its own phase, which read phasing decides like any other site's.

## Re-genotyping from the phase

In the site likelihood, every read of a site has the same [mixture
weights](read-likelihood-direct-genotyping.md#mixture-weights): the probabilities that it came from
each strand, before its bases at the site are seen. Once sites are phased, a read's alleles at the
other phaseable sites it spans show which strand it more likely came from. Re-genotyping uses this
to give each read its own weights, tilted towards the strand its other sites place it on, and so
corrects each site's likelihoods. The quantities $e_r$, $p_{ra}$, $q_{rt}$, $c_{rt}$ and the
allele-length weights $v$ are those of [Read phasing](#read-phasing). $\sigma$ is the logistic
function, and $\mathrm{logit}$ is its inverse.

### Strand log-odds

A read's **strand log-odds** at site $s$, $\Lambda_{rs}$, measures how strongly its other sites
place it on strand 0. Each other phaseable site $t$ at which $r$ is used adds one term: the
natural-log odds that the read came from strand 0, with $q_{rt}$ taken in the phase read phasing
chose, from the mate the site keeps.

$$
\Lambda_{rs} = \sum_{t \neq s} \ln \frac{c_{rt} q_{rt} + (1 - c_{rt})/2}{c_{rt}(1 - q_{rt}) + (1 - c_{rt})/2}
$$

Leaving out $s$ keeps a site from confirming its own genotype. Each phase set labels its strands
independently, so a read's strand is usable only at the sites of the phase set it was used in. Its
$\Lambda_{rs}$ is 0 at a site of another phase set, and at every site if it was used in more than
one phase set. At a site that the linkage model did not phase, and so has no phase set, any read
used in one phase set keeps its $\Lambda_{rs}$.

### Tempering

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
$\Lambda_{rs} \neq 0$ and $q_{rs} \neq 1/2$. Each of the two points to a strand: $\Lambda_{rs}$ to
strand 0 when it is positive, $q_{rs}$ to strand 0 when it is above $1/2$, and each to strand 1
otherwise. The observation records whether they point to the same strand. The observations are
sorted by $\vert \Lambda_{rs} \vert$ and grouped into bins. Every bin but the last holds the larger
of a fixed number and a fixed fraction of the observations, and the last holds the remainder (see
[Fixed constants](read-likelihood-genotyping.md#fixed-constants)). For each bin, the tempered
probability of the strand, $C \sigma(\tau \vert \Lambda \vert) + (1 - C)/2$ at the bin's mean $\vert
\Lambda \vert$, predicts its rate of agreement. $\tau$ is chosen from a fixed grid, with $C$ held at
its given value, to minimise the squared difference between the predicted and observed rates, with
each bin weighted by its size. The fit treats the strand that $q_{rs}$ points to as true, so errors
in $q_{rs}$ lower the fitted $\tau$. With too few observations, $\tau$ is 0 and re-genotyping leaves
the likelihoods as they were.

### Likelihood correction

Take a genotype of two different alleles, $a$ on strand 0 and $b$ on strand 1, whose allele-length
weights are $v_a$ and $v_b$. Each read gets its own weights

$$
\pi_{ra} = \frac{v_a e^{y_{rs}}}{v_a e^{y_{rs}} + v_b}, \qquad \pi_{rb} = 1 - \pi_{ra}
$$

The correction added to $\ln \mathcal{L}(G)$ sums over the site's reads $R$:

$$
\sum_{r \in R} \left[ \ln\left((1 - e_r)(\pi_{ra} p_{ra} + \pi_{rb} p_{rb}) + e_r\right) - \ln\left((1 - e_r)(v_a p_{ra} + v_b p_{rb}) + e_r\right) \right]
$$

A genotype other than the chosen one has no phase, so each heterozygous genotype, the chosen one
included, is scored under both assignments of its alleles to the strands, and the larger correction
is kept. Homozygous genotypes are unchanged and have no assignment to choose, so this choice can
only favour heterozygous genotypes.

Both terms of the correction use the allele-length weights $v$, while the site likelihood's read
term uses the mixture weights. The two agree when the alleles have the same length and neither
visits a node twice; otherwise the corrected likelihood only approximates the one with per-read
weights.

### Nested haploid chains

At a site of a nested haploid chain, each genotype is a single allele on one strand, so there are no
weights to tilt. Instead, a read is made less informative when the phase places it on the other
strand. Let $y$ be the read's tempered strand log-odds towards the chain's strand, and
$\eta = \min(1, e^{y})$. Each of the read's relative likelihoods $p$ becomes $\eta p + 1 - \eta$.
$\eta$ is 1 when $\tau = 0$ or when the read points to the chain's strand, so those reads are
unchanged. `--no-regeno-haploid` turns this off.

# Read-likelihood genotyping (`vg call --read-likelihood`)

`vg call` has two genotypers. The default one aggregates read *support* — how many reads cover
each allele — and judges a genotype against those counts. `--read-likelihood` instead asks, for
every read and every candidate allele, *how well does this read fit this allele*, and genotypes
from an explicit `P(reads | genotype)`.

This page specifies that model completely: every term in the objective, where each number comes
from, every parameter that changes it, and what the reported qualities do and do not mean.

- [The objective](#the-objective)
- [The read term](#the-read-term)
- [The depth term](#the-depth-term)
- [Where the fit scores come from](#where-the-fit-scores-come-from)
- [Which alleles are considered](#which-alleles-are-considered)
- [The linkage layer](#the-linkage-layer)
- [Genotype quality and the VCF fields](#genotype-quality-and-the-vcf-fields)
- [Parameter reference](#parameter-reference)
- [What this model does not give you](#what-this-model-does-not-give-you)

## The objective

A genotype `G` is a multiset of allele indices of size *ploidy* (`{0,0}`, `{0,1}`, `{1,2}`, …).
Every genotype of the site's ploidy over its alleles is enumerated and scored — exhaustively, with
no pruning. The score is:

```
ln P(reads | G) =  Σ_r  ln [ (1 − e_r) · Σ_{h ∈ G} w_h · rel(r,h)  +  e_r ]      (read term)
                 +  w_d · ln Poisson( N_eff ; λ_G )                              (depth term)
```

Both terms are on by default. A description of only the first describes a configuration nobody
runs.

The sum over `r` runs over reads overlapping the site. The inner sum over `h` runs over the
haplotypes of the genotype — twice over the same allele for a homozygote.

**Which genotype is then called depends on one more thing.** With `--linkage-weight 0` it is
simply the argmax of the score above. At the default `--linkage-weight 2` the scores become the
emissions of an HMM across neighbouring sites, and the call is the argmax of *that* posterior —
see [the linkage layer](#the-linkage-layer). The expression above is unchanged either way; it is
what the linkage layer consumes, not something the linkage layer replaces.

## The read term

### `rel(r,h)` — the read's fit to an allele

`rel(r,h) ∈ [0,1]` is read `r`'s likelihood under allele `h`, **divided by that read's own best
score at this site**. Every read's best-fitting allele therefore scores exactly 1, and the others
are expressed relative to it.

`rel(r,h) = 0` means allele `h` cannot place the read at all. That is strong evidence *against*
`h`, not missing data, and the model treats it that way.

Why normalise per read rather than per allele. The scoring produces `ln P(read | allele)` only up
to an unknown per-read constant. That constant cancels when comparing genotypes — but only while
nothing else enters the per-read expression, and the mismapping term does. Dividing each row by
its own maximum puts every component on one dimensionless scale, so the constant cancels by
construction. It also removes the numerical hazards: "no valid placement" is `0`, not `−∞`.

Do **not** read `rel` as a distribution over alleles. Normalising the other way — making a read's
alleles sum to 1 — gives a posterior under a uniform prior, which double counts the prior and
destroys the absolute-fit information the mismapping term needs.

**Consequence worth internalising:** because each row is divided by its own maximum, anything
that improves a read's fit to *all* of a site's alleles equally divides straight out. Only
*allele-differential* changes can move a genotype.

### `e_r` — the mismapping probability

```
e_r = clamp( phred_to_prob(MAPQ_r), --mismap-min, --mismap-max )      default clamp [0.02, 0.7]
```

`e_r` is the probability that read `r` did not come from this locus at all. Such a read explains
every genotype equally well, which is why it enters as an additive background rather than as a
weight.

Because the bracket is `(1 − e_r)·(…) + e_r` and the mixture is in `[0,1]`, the bracket lies in
`[e_r, 1]`. Its log is therefore bounded below by `ln(e_r)`: **no single read can penalise a
genotype without limit.** That bound is the entire purpose of the floor.

This is why `--mismap-min` is not really "the minimum mismapping rate". At the original floor of
1e-8 a single read that matched nothing contributed −18.4 nats, enough to veto a genotype on its
own. MAPQ measures confidence that a read is in the right *place*; it says nothing about whether
the read's path through this particular site is right. A locally misaligned read is still MAPQ 60.
The floor is what caps that read's veto, so `--mismap-min` reads as *"P(this read's evidence here
is unreliable, for any reason)"*.

`--no-mismap-term` does not set `e_r = 0`. It sets `e_r = --mismap-min` for every read, removing
the MAPQ *dependence* while keeping the floor — because removing the bound entirely is not a
configuration you want.

### `w_h` — the mixture weights

`w_h` is haplotype `h`'s share of the reads, summing to 1 over `G`. By default:

```
w_h ∝ U_h + R − 1
```

where `U_h` is the length of sequence unique to `h` among the genotype's members, and `R` is the
mean read length. That quantity is the number of read start positions that can produce a read
capable of *distinguishing* `h` — so the weight is proportional to the evidence `h` should
generate.

`--flat-mixture` restores `w_h = 1/|G|`. That is wrong wherever the alleles differ in length: it
asserts both haplotypes produced half the reads, which over a heterozygous deletion is false, and
losing large heterozygous deletions is exactly what it caused. Where alleles are the same length
the two are identical — equal lengths give exactly 1/2 — so SNV calls are unchanged bit for bit.

Note `U_h` is *unique* content, not whole traversal length. Weighting by whole traversal length
was tried and is worse.

## The depth term

The read term has no opinion about reads that are **absent**. That is precisely the evidence a
homozygous deletion presents: not reads arguing for the deletion, but the missing reads that the
reference allele would have required. The depth term supplies it.

```
w_d · ln Poisson( N_eff ; λ_G )

N_eff = Σ_r (1 − e_r)                              expected reads genuinely from this locus
λ_G   = rate · Σ_{h ∈ G} max(L_h + R − 1, 1)       expected count under G
rate  = reads per haplotype per bp, from a local window
w_d   = --depth-term, default 0.1
```

`L_h` is allele `h`'s traversal length and `R` the mean read length, so `L_h + R − 1` is again the
number of start positions producing an overlapping read.

`rate` is measured over the read source's own fetch window on the same contig — thousands of node
IDs wide, and already fetched and cached to answer the site's own query, so it costs no extra I/O.
Crucially both sides of the comparison weight a read by `(1 − e_r)`, so the correction is
*relative*: a site whose mapping quality matches its neighbourhood's is unaffected, and only a
site that differs from its surroundings moves. `--depth-count-raw` counts whole reads on both
sides instead.

`w_d` is small (0.1) because the depth term is a far cruder statistic than the read term and is
dominated by it wherever reads are informative. `--depth-term 0` disables the term.

### `DR`, the depth ratio

The two quantities above also give the one diagnostic this model reports about depth. `DR` is
their ratio, evaluated at the **called** genotype:

```
DR = N_eff / λ_G
```

So `DR = 1.0` means the site saw exactly as many reads as the call implies. Because both sides
weight a read by `(1 − e_r)`, it inherits the same *relative* property as the term itself: a site
whose mapping quality matches its neighbourhood's sits at 1 regardless of how good that
neighbourhood is, and only a site out of step with its surroundings moves.

Reading it:

- **Well above 1** — more reads than the call can account for. Collapsed repeats, where reads from
  several copies of a sequence pile onto the one locus that represents them.
- **Near 0.5** — the call claims about twice the sequence actually covered. This is what a missed
  heterozygous deletion looks like: the genotype asserts two long alleles where one is short.

`DR` is reported **whether or not `--depth-term` is armed**, so the signal can be measured before
it is trusted, and it is what `--depth-quality` scales `GQ` by. It is omitted from a record
entirely — rather than emitted as a sentinel — where the called genotype predicts no reads at all
and the ratio is undefined.

## Where the fit scores come from

The raw `ln P(read | allele)` is read off the **alignment the mapper already produced**, not
recomputed by dynamic programming. For each allele the read's path is walked against the allele's
path: on shared nodes the read's own edits are scored, and where the paths diverge the read's
bases are scored against the allele's bases directly. A read that cannot be placed on an allele
at all scores `−∞` there, which normalises to `rel = 0`.

This is what makes the model fast enough to score every enumerated allele exhaustively with no
pruning and no read subsampling — the cost is roughly 10²–10³× below alignment.

It is also the model's main approximation. Reading an alignment off against a different allele is
not the same as aligning to it, and the two differ whenever the alleles differ by a shift. Both a
bounded-shift correction and full WFA realignment were implemented and measured; neither paid,
because the row normalisation divides out any improvement common to a site's alleles, and optimal
realignment changed which allele a read preferred in 40 cases out of 91,914.

## Which alleles are considered

The model scores whatever the traversal finder enumerates. Two sources:

- **Panel enumeration** (`GBWTTraversalFinder`) — the traversals the haplotype panel actually
  spells. **This is the default under `--read-likelihood`** on a GBZ carrying at least two
  haplotype-sense paths, and it needs no pack file. `-z` selects it explicitly and is redundant
  here; `-g FILE` selects it from a separate GBWT index rather than the GBZ's own.
- **Support enumeration** (`FlowTraversalFinder`) — traversals found by flow through node and edge
  weights from a `vg pack` file. Selected by `--enumerate-support`, which then requires `-k`.

Panel enumeration measured better on every small-variant class tested and on three of four
structural-variant sets, which is why it is the default here. It is *not* the default for the
support-based caller, where it loses SV F1.

The limit is structural and worth knowing: **panel enumeration can never propose an allele no
haplotype carries.** Its ceiling is the panel's content. On a panel that represents your sample's
variation well this is a large net win; on a thin or poorly matched panel it is a recall ceiling,
and `--enumerate-support` is the way out. A panel of fewer than two haplotypes falls back to
support enumeration automatically.

## The linkage layer

`--linkage-weight` adds a layer that is **not part of the expression above**. The per-site model
is unchanged; linkage re-decides genotypes afterwards, using `ln P(reads | G)` as its emission.

The motivation: each snarl is genotyped independently, so the emitted call set is a concatenation
of per-site argmaxes — a pair of haplotypes free to switch panel haplotype at every site at no
cost. Measured on chr20 against a 33-haplotype panel, the typical adjacent called pair is 1.8×
more likely than independence predicts, and pairs no single panel haplotype carries are nine times
commoner where the reads were undecided (2.8% at GQ < 10) than where they were confident (0.3% at
GQ ≥ 40). The information is real and concentrated exactly where the per-site likelihood is flat.

The structure follows PanGenie (Ebler et al. 2022): hidden states are **ordered pairs** of panel
haplotypes plus a wildcard, `(H+1)²` states; transitions are Li–Stephens over reference distance;
genotypes come from forward–backward posteriors summed over the states implying each genotype.
The emission is this caller's own likelihood rather than PanGenie's k-mer model.

Per-strand switch probability between sites `g` bp apart:

```
ρ(g) = [ ρ_min + (1 − ρ_min)(1 − exp(−g / --linkage-scale)) ] ^ --linkage-weight
```

Tempering the switch *probability* by an exponent, rather than scaling a log-transition, keeps it
a probability at every setting. `--linkage-weight 0` makes every transition uniform, the chain
memoryless, and the posterior collapses to the per-site likelihood — the caller is then bit-for-bit
unchanged.

A **wildcard haplotype** may take any allele at any site. States involving it are discounted by
`escape` (1e-2) in the emission, once per wildcard strand, and its allele is marginalised
*conditional on the reads* rather than spread uniformly — spreading it uniformly ignores the
emission, and at a site whose reads are decisive for an allele the panel does not carry, that
hands the win to a genotype the reads had ruled out.

The wildcard is not optional. A state implies a genotype, so without it any genotype no panel pair
can spell is unreachable — and the graph need not contain the sample being genotyped. Omitting it
makes the model suppress novel alleles, which presents as a precision gain while destroying recall.

### `--linkage-freq-prior` — the allele-frequency prior

This is the dominant parameter of the linkage layer, and the least self-explanatory, because it
does not add a prior. It controls **how much of a prior the state space was already imposing.**

**Where the prior comes from.** A genotype is not a state; it is what a state *implies*. The
posterior for a genotype is therefore the sum over every ordered haplotype pair spelling it — and
a common allele is carried by many panel haplotypes, so genotypes containing it are spelled by
many more pairs. Summing states silently weights genotypes by panel allele frequency.

For example, take a 32-haplotype panel where allele 1 is carried by 20 haplotypes and allele 0 by
12. The number of ordered pairs spelling each genotype is then:

```
(0,1) → 12·20 + 20·12 = 480       (1,1) → 20·20 = 400       (0,0) → 12·12 = 144
```

Before a single read is consulted, the state space already prefers `0/1` over `0/0` by 480:144.
That is a real allele-frequency prior, and it is not optional — it is a consequence of counting
states at all.

**What the parameter does.** The summed mass is divided back out by `multiplicity^(1−F)`. Writing
the sum over the pairs spelling a genotype as `multiplicity × mean_pair_mass`, the surviving mass
is exactly:

```
posterior mass  =  multiplicity^F  ×  mean_pair_mass      (mass from known-known states)
```

so `F` is the exponent on the panel count and nothing else. `mean_pair_mass` is the forward–backward
mass an average such state carries, which is where the reads and the transitions enter.

| `F` | Effect |
|---|---|
| `0` | Prior removed entirely. Each genotype counted once, however many pairs spell it. |
| `1` | Divides by `multiplicity^0 = 1` — a no-op. The state space's own count stands. |
| `> 1` | Amplified **beyond** what the state space implies. `F` is an exponent, and nothing about it stops at 1. |

The default is `5`. On the example above that takes the `0/1`:`0/0` prior ratio from 3.33 to
`3.33^5 ≈ 410`, or 6.0 nats. Carrying the same example up the axis shows why the parameter has an
optimum at all rather than improving indefinitely:

| `F` | prior ratio | nats |
|---|---|---|
| 1 | 3.3 | 1.2 |
| 5 (default) | 411 | 6.0 |
| 8 | 15,200 | 9.6 |
| 12 | 1,880,000 | 14.5 |

The prior grows geometrically while the evidence a site carries does not, so past some point it
stops informing the reads and starts overruling them. That the measured inversion (below) lands
between 8 and 12 is consistent with this arithmetic; the arithmetic alone does not predict where,
since it depends on how much the reads at a given site actually say.

**It leaks through two channels, and both are corrected.** Besides pairs of known haplotypes, a
state pairing a known haplotype with the wildcard also occurs once per haplotype carrying that
allele — so `carriers[k]`, the number of panel haplotypes carrying allele `k`, weights the
half-wildcard mass exactly as multiplicity weights the known-known mass. Both are divided by their
own count to the `1−F`. A test caught this: neutralising only the first left a site whose reads
were flat still tilted toward the common allele, which is precisely what `F = 0` is supposed to
prevent. States with *both* strands latent have no multiplicity to correct and are left alone.

**Measured behaviour.** On chr20 and chr6 against 34-haplotype panels, crossed against the
transition weight, it improves every variant class at every weight, monotonically, through 1 and
well past it, peaking near 5–8. At the joint optimum (`weight 2`, `freq_prior 5`) it accounts for
most of the layer's total gain — small-variant genotype F1 +0.0099 on chr20 and +0.0074 on chr6
against no linkage at all, against +0.0047 and +0.0036 for the transition model alone. Beyond 8 it
inverts: by 12 the prior overwhelms the reads, SNV F1 falls below the no-linkage baseline and
structural-variant recall collapses.

The effect is almost entirely small indels — deletions most, insertions next, SNVs flat to within
0.0002 across the whole axis. That is the mechanism showing through: SNVs are already settled by
the reads, so a prior has nothing to act on.

**Caveat, not retired.** Every one of those numbers is from a 34-haplotype panel. Multiplicity is a
far coarser statistic over three haplotypes than over thirty-three, and a large exponent over a
count that barely varies is not the same operation. On the two 4-haplotype graphs the difference
between `F = 0` and `F = 5` is under 0.0004 on every class, in both directions — safe there, and
not useful there.

One structural note for anyone reading the code: `LinkageModel::Params` defaults this to `0`, so a
model constructed directly is the plain HMM with no prior, which is what the unit tests compare
against. `vg call` defaults its flag to `5`.

Inference runs over windows of `2000` sites with a `250`-site margin discarded at each end, because
exact inference over a whole chain would serialise a caller that is otherwise parallel over snarls.
Linkage is spent by 10–30 kb, so this is near-exact.

The layer needs panel enumeration, and declines silently without it. Asking for it explicitly where
it cannot run is a hard error rather than a silent no-op. Its benefit is concentrated on panels of
tens of haplotypes and is close to inert on thin ones.

## Genotype quality and the VCF fields

| Field | Meaning |
|---|---|
| `GT` | The called genotype: the argmax of the objective, or of the linkage posterior when that layer is on. |
| `GL` | `log10 P(reads \| genotype)` for every genotype, in VCF order. |
| `GQ` | Phred gap between best and second-best genotype, **scaled by the fraction of reads the call explains**. |
| `GQI` | The same gap *without* that scaling. Equals `GQ` under `--no-share-quality`. |
| `GP` | Log-scaled posterior of the called genotype under a uniform prior. |
| `DP` | Reads overlapping the site. |
| `AD` | Reads whose best-fitting allele is this one. **Does not sum to `DP`** — see below. |
| `DR` | `N_eff / λ_G` at the called genotype — see [above](#dr-the-depth-ratio). `1.0` means the read count matches the call. Absent where the ratio is undefined. |
| `BL` | Mean over reads of the best raw alignment score any allele gave them. |

**Why `GQ` is scaled.** A read whose best-fitting allele lies *outside* the called genotype enters
every genotype's likelihood and therefore cancels out of the best-versus-second-best comparison.
The likelihood ratio simply cannot see it. Scaling by the explained share restores it. The
consequence is that `GQ` here is a *quality score for ranking*, not a calibrated posterior;
`GQI` is the pure likelihood ratio for anyone who wants it.

**Why `AD` does not sum to `DP`.** Two reasons: a read fitting several alleles equally splits its
vote, and — much more importantly — only alleles that reached the VCF record get a column, while
the genotyper scored every allele the site offered. Where many alleles were enumerated and few
emitted, most reads best-fit something absent from the record. That shortfall is itself
informative: it is how much of the evidence the emitted alleles fail to explain.

**`--depth-quality A`** scales `GQ` by `exp(−A · |ln DR|)` at records whose called alleles change
length by at least 50 bp, so a call whose read count is implausible for the sequence it claims
ranks lower. Off by default, and affects ranking only. See [`DR`, the depth
ratio](#dr-the-depth-ratio) for how to read the value it acts on.

## What this model does not give you

**`GL`, `GQ` and `GP` are not calibrated probabilities.** Reads are treated as independent. Mates
overlapping the same site are not independent, so the product accumulates confidence like `R`
rather than `√R`. Derived qualities are therefore over-confident, and increasingly so with depth.
Use them for ranking and filtering, not as posteriors.

**`rel` values are relative, not absolute.** They are valid for argmax, for a likelihood
difference such as `GQ`, and for a posterior over genotypes. They are not absolute probabilities
of a read.

**Only allele-differential evidence matters.** Anything improving a read's fit to every allele at
a site equally cancels in the row normalisation. This is a property of the design, not a defect,
but it explains why several plausible-sounding improvements measure as no-ops.

**Scope of the tuning.** Every default here was set on HG002 against HPRC graphs, chr20 and chr6,
at 4 and 34 haplotypes — four datasets, one sample. The defaults are defensible on that evidence
and it is not whole-genome, not multi-sample, and the panel-dependent parameters
(`--linkage-*`, and panel enumeration itself) are the ones most likely to behave differently on a
sample the panel represents poorly.

## Parameter reference

Every flag below requires `--read-likelihood`; `vg call` errors if one is passed without it.

### Reads in — exactly one source required

| Flag | Default | Effect |
|---|---|---|
| `--gam FILE` | — | Read alignments, held in memory. |
| `--gaf-reads FILE` | — | The same, as GAF. |
| `--gam-index FILE` | — | `.gai` for `--gam` (`vg gamsort -i`), so reads are fetched per site instead of all held in memory. |
| `--gaf-base FILE` | — | GAF-Base database, queried per site via `gbz-base`. |
| `--gbz-base FILE` | input graph | Graph to resolve `--gaf-base` queries against. |
| `--read-window N` | 4096 (`--gaf-base`), 256 (`--gam-index`) | Node-ID window per indexed fetch. |
| `--read-min-mapq N` | 0 | Drop reads below this MAPQ outright. |

### Allele enumeration

| Flag | Default | Effect |
|---|---|---|
| `--enumerate-support` | off | Enumerate from read support instead of the haplotype panel. Requires `-k/--pack`. |

### Model terms

| Flag | Default | Effect |
|---|---|---|
| `--depth-term W` | 0.1 | Weight `w_d` on the Poisson depth term. 0 disables it. |
| `--depth-count-raw` | off | Count whole reads rather than `1 − e_r` in `N_eff` and `DR`. |
| `--mismap-max P` | 0.7 | Upper clamp on `e_r`. Governs how far a low-MAPQ read is discounted; matters most on haplotype-rich graphs. |
| `--mismap-min P` | 0.02 | Lower clamp on `e_r`, and so the bound on any one read's veto, `ln(P)`. |
| `--no-mismap-term` | off | Fix `e_r = --mismap-min` for every read, removing the MAPQ dependence. |
| `--flat-mixture` | off | `w_h = 1/ploidy` instead of length-weighted. Restores the pre-correction model exactly. |

### Linkage between sites — needs panel enumeration

| Flag | Default | Effect |
|---|---|---|
| `--linkage-weight W` | 2 | Exponent tempering the switch probability. 0 is off and reproduces the per-site caller exactly. |
| `--linkage-scale N` | 10000 | Distance scale of the linkage decay, in bp. Nearly flat from 10–40 kb. |
| `--linkage-freq-prior F` | 5 | Exponent on the state space's implied allele-frequency prior. 0 removes it, 1 keeps it as presented, and it inverts past ~8. |

### Quality reporting — ranking only, never changes a genotype

| Flag | Default | Effect |
|---|---|---|
| `--no-share-quality` | off | Report `GQ` as the raw likelihood ratio, so `GQ == GQI`. |
| `--depth-quality A` | 0 (off) | Scale `GQ` by `exp(−A·|ln DR|)` at records whose called alleles change length by ≥ 50 bp. |

### Debugging

| Flag | Default | Effect |
|---|---|---|
| `--dump-likelihoods F` | — | Write the per-site reads × alleles matrix to `F` as TSV. |

## Worked invocation

The simplest form, on a GBZ carrying a haplotype panel — no pack file needed, panel enumeration
and the linkage layer both on by default:

```
vg call graph.gbz --read-likelihood --gam reads.gam -t 8 > calls.vcf
```

The like-for-like comparison against the support caller, holding enumeration constant so that
only the genotyping model differs:

```
vg call graph.gbz --read-likelihood --enumerate-support -k reads.pack --gam reads.gam > calls.vcf
```

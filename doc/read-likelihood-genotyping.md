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
- [Phasing and the mosaic genome](#phasing-and-the-mosaic-genome)
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

`--linkage-weight` adds a layer that is **not part of the expression above**. The per-site model is
unchanged; linkage decides the genotype, using `ln P(reads | G)` as its emission.

"Decides", not "re-decides afterwards": the per-site argmax is an input to that decision, not an
answer the layer edits. Sites are genotyped while their reads are resident and staged; the layer then
settles every generation, deepest last; and each record is built from the settled genotype. Nothing is
written before it is decided, so there is no afterwards to revise in.

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

The two strands are transitioned independently and the joint transition is their product, which is
what keeps a step at O(m²): the sum over all `(a′,b′)` predecessors collapses into one row sum, one
column sum and one grand total. A jump redraws a strand's panel haplotype **uniformly over all m**,
including the one it is already on, so ρ is the probability of forgetting and redrawing rather than of
changing haplotype — the chance of actually landing elsewhere is `ρ·(m−1)/m`.

**A known limitation of this form.** Because the exponent is applied to a per-step *probability*, ρ
does not compose across steps: the same span divided into more sites is stickier than the same span
taken in one — measured at 5.8× for a 1,348 bp span split seven ways at the default weight. The
model's effective switch rate therefore depends on how many sites are in the chain as well as on
genomic distance, so a change to which sites enter the layer re-tunes the linkage strength without
anyone touching a parameter. A form that composes is `1 − ρ = ((1 − ρ_min)·exp(−g/scale))^(1/weight)`;
adopting it would change what the weight's numeric value means and needs a refit.

**Distance below the top level.** The step between adjacent sites is the reference difference,
below the top level as at it. Measuring it along the parent's settled traversal instead was built
and tested across three arms spanning no-link to perfect-link, and came out unmeasurable — the two
alternatives disagree in sign between chr20 and chr6 — so the machinery was removed rather than kept
behind a switch. A site with no reference position transitions uniformly.

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

## Phasing and the mosaic genome

The linkage layer's genotype decoding is *marginal* — each site's posterior is summed over the
states implying each genotype, and the argmax taken independently. That yields genotypes, but a
sequence of per-site argmaxes need not be spellable by any single pair of haplotypes, so it is not
a phasing.

Phasing runs a second decoding of the same model: **max-product**, giving the single most probable
path of haplotype pairs through the chain. The order of that pair is the phase. It is on wherever the
linkage layer runs, so the invocation is just:

```
vg call graph.gbz --read-likelihood --gam reads.gam > calls.vcf
```

`--no-phased` turns it off. `--phased` still exists and asks for it explicitly, which changes only
what happens when it cannot be delivered: the default declines quietly, an explicit request errors.

**A caveat that belongs beside every phased number here.** Blocks are chromosome-length, but the
switch rate is 2.4% per adjacent heterozygous pair, so the orientation re-randomises every forty sites
or so and blockwise Hamming sits near 48%. A long block says which sites share one `PS`, not that the
phase is trustworthy across a chromosome. Nested haploid records are worse again -- they contribute a
41% switch rate against that 2.4% baseline -- so their strand assignment should be read as provisional.

`GT` becomes `0|1`, and `FORMAT/PS` names the phase block. (Note the existing `INFO/PS` under `-A`
is vg's parent-snarl pointer — a different field in a different namespace.)

**Constrained to the emitted calls.** At each site the states are restricted to those spelling the
genotype actually written out, so phasing re-orders a call without re-deciding it. Unconstrained
Viterbi maximises path probability and will contradict the VCF; that is not wanted, so it is not
what runs. Feasibility is guaranteed under panel enumeration — if allele *i* is carried by some
haplotype and *j* by another, the pair spelling *i/j* exists — and the wildcard covers the rest.
Measured, the constraint costs essentially nothing: 0.98× the switches on chr20, 0.96× on chr6.

**Phase comes from the panel, not from reads.** A read-based phaser links two heterozygous sites
only when a read or fragment spans both, so its blocks are read-length. Here linkage carries as far
as the transition model allows, so **a phase block is a whole chain** — chromosome-scale. That is a
much stronger claim than a read-based phaser makes, and switch error has to be read beside it,
because shorter blocks make switch error small for free.

Measured against the phased T2T-Q100 HG002 truth, with HG002 excluded from the graph:

| graph | panel | switch error | block N50 |
|---|---|---|---|
| chr20, 4 haplotypes | 4 | 3.43% | 66.2 Mb |
| chr20, 34 haplotypes | 34 | 2.30% | 66.2 Mb |
| chr6, 34 haplotypes | 34 | 1.74% | 172.1 Mb |

Phasing quality is panel quality: the same chromosome with a richer panel loses a third of its
switch errors. Expect a sample the panel represents poorly to do worse.

### `--mosaic-out FILE`

The same path, written as the genome it implies. The mosaic is piecewise — 1.75% of chr20's site
slots begin a new segment — so it is run-length encoded, one line per maximal run of one strand on
one panel haplotype:

```
#mosaic-version	5
#patch	reference
#nested	kept
#unexplained	connected
#graph	chr20_0_chr20.gbz
#sample	HG002
#reference	CHM13#0#chr20
#decoding	constrained-viterbi
#note	ref_start/ref_end are advisory, in the #reference coordinate system; ...
#note	segments are maximal runs on one panel haplotype; ...
#note	hap_index is internal to this run; ...
#haplotype	0	CHM13#0
#haplotype	9	recombination#30
	... one line per panel haplotype, 34 of them on this graph
#note	gbwt_node/gbwt_offset is the GBWT position of the haplotype at start_node: ...
#note	a segment never spans a GBWT fragment boundary, ...
#note	a fragment is a PATH: its rows join end to start, on the same oriented node, ...
#note	start_node and end_node are ORIENTED node ids, id * 2 + is_reverse ...
#H	contig	strand	fragment	ref_start	ref_end	start_node	end_node	hap_index	haplotype	sites	gbwt_offset
H	chr20	0	0	22	533	229637742	229638340	4	recombination#19	16	19
H	chr20	0	0	945	2268	229638340	229641408	9	recombination#30	14	9
```

Tab-separated. Lines beginning `#` are header or comment; every data line begins `H`.

**Parsing the header.** Key on the first field, not on line position. `#note` lines carry prose for a
human and no data — they are truncated above — and they appear both before and after the
`#haplotype` table rather than in one block, so a parser that counts lines to find the table will
break. Treat an unrecognised `#` key as a comment and skip it: downstream tools extend this header
(the whole-genome assembly in the companion eval repo replaces `#graph`/`#reference` with a
per-contig `#contig` table, because per-contig runs do not share one graph). The keys `vg call`
itself writes are exactly `#mosaic-version`, `#graph`, `#sample`, `#reference`, `#decoding`, `#note`,
`#haplotype` and `#H`.

**Format version.** `#mosaic-version` is `5`, and unlike 3 and 4 it is **not** a compatible
extension: `end_node` changed meaning and the column set changed, so a version 4 reader must refuse
a version 5 file.

**A strand is one walk.** `end_node` is the *next* segment's `start_node`, so consecutive segments of
one fragment meet at the same oriented node and concatenate -- counting each junction once -- into a
single path. That is the invariant the format exists for; version 4 met it at 3% of boundaries.
`(contig, strand, fragment)` is the path identity, and one path per strand is the normal case.

**A fragment is a path or it is nothing.** Two things break a walk, and both end the fragment rather
than being written into it:

- *A row that spans no graph.* A nested snarl is CONTAINED in its parent, not sequential with it, so
  a recombination at a nested site can open a segment whose start the right-extension then moves to
  the container's end -- past the row's own sites. chr20 had two, at 32,599,461 and 32,599,482,
  claiming two and four called sites across no graph at all, and they were the only rows the round
  trip could not expand. The ROW is dropped and its SITES are kept: removing the sites instead moved
  run boundaries elsewhere and took chr20 from 3 fragments to 7.
- *A row the carried direction cannot walk.* chr20 has one, an inversion at 32,709,971-32,717,434 --
  node 119177446 down to 119168165, backwards in node id, which the haplotype traverses REVERSE at
  both ends. The carry is authoritative precisely so a row cannot pick an orientation that breaks
  contiguity, but where the row would otherwise carry no position at all there is nothing left to
  break, so the other direction is tried and the row stands alone in its own fragment. It cannot
  join its neighbours -- that is what an inversion boundary is, and X+ followed by X- is not a walk.

The result on chr20 is 9 fragments, every junction met, no row without a position, and every
fragment expanding to an exact walk in the graph.

**Which haplotype covers the stretch between two segments' sites is arbitrary.** No called site lies
in it, so nothing distinguishes the earlier haplotype from the later, and a recombination anywhere
inside is equally consistent. Extending the earlier one is a convention: the crossover is
*bracketed* by that stretch, not located within it.

**Three options, each recorded in the header so a consumer never infers which way the file was
written:**

| header | default | the other setting |
|---|---|---|
| `#patch` | `reference` -- a gap no panel haplotype spans is filled with the reference, marked `ref` | `none` (`--no-mosaic-patch-gaps`) leaves the gap and breaks the fragment |
| `#nested` | `kept` -- a haplotype switch inside a nested chain is its own segment | `merged` (`--no-mosaic-nested`) folds it into the enclosing run: fewer switches, the parent's route through the child snarl |
| `#unexplained` | `connected` -- a stretch the panel cannot explain is crossed by carrying the flanking haplotype through | `broken` (`--mosaic-break-unexplained`) emits an unwalkable `*` row alone in its own fragment |

The last two are approximations, and deliberately so: each reconstructs the carried haplotype's
sequence across a few sites rather than the called alleles, in exchange for a contiguous path. On
chr20 that is 12,225 sites for `#nested merged` and 21 for `#unexplained connected`.

| column | meaning |
|---|---|
| `contig` | reference contig the segment's sites sit on |
| `strand` | `0` or `1`. Strand identity is consistent along a contig and matches the order of the phased `GT` in the VCF: `strand 0` carries the allele left of the `\|`. |
| `ref_start` | reference position of the **first site** in the run — not the start of the sequence the segment covers. In the `#reference` coordinate system, and **advisory**. |
| `ref_end` | reference position of the **last site** in the run |
| `start_node` | snarl start node of the first site. Authoritative. |
| `end_node` | snarl end node of the last site. Authoritative. |
| `hap_index` | row in the `#haplotype` table. Internal to this run — see below. |
| `haplotype` | `sample#phase`, or `*` where the panel does not explain the strand |
| `sites` | number of called sites in the run |
| `gbwt_node` | oriented GBWT node at `start_node`: `start_node * 2 + is_reverse`. `.` if unresolvable. |
| `gbwt_offset` | offset within that node's GBWT record, so `(start_node, gbwt_offset)` is a GBWT position outright -- `extract()` it and follow `LF()` to `end_node`, with no locate and no r-index |

`hap_index` takes four values. An **integer** names a panel haplotype by its `#haplotype` row.
**`ref`** marks a stretch the reference was substituted across: with `sites` = `.` it filled a gap
between segments, and with a site count it replaced a haplotype the graph does not carry across that
segment. **`*`** appears only under `#unexplained broken` and is the one row kind that is not
walkable. `nested_sites` and `max_depth` were version 4 diagnostics and are gone; they are derivable
from the VCF's per-record depths.

**Ordering.** Grouped by contig; within a contig all `strand 0` lines precede all `strand 1` lines;
within a strand, `ref_start` increases. Segments on one strand partition that contig's sites
exactly — every site belongs to exactly one segment, so the `sites` column sums to the site count.

**Reconstruction, and why the GBWT position is there.** `extract({gbwt_node, gbwt_offset})` gives
the haplotype's path from that point; follow `LF()` to `end_node`. No search is involved, and that
is the point of carrying the position at all. A node plus a haplotype *name* is not enough to walk
from: turning a name into a path in the GBWT at a given node is a `locate()`, which needs a
document-array sample or an r-index — neither of which a plain GBZ carries — and otherwise degrades
to scanning the component. This project's chr20 r-index is 86 MB against a 77 MB GBZ, so that index
costs more than the graph it indexes. Resolving the position once at write time costs the caller a
few thousand lookups and gives every reader an O(1) entry point instead.

Two consequences follow from the position being a *position*:

- **A segment never spans a GBWT fragment boundary.** One position can only walk one fragment, so
  the caller splits a run wherever the underlying GBWT path ends and another begins. A haplotype
  present as several fragments therefore yields several segments. On chr20 this splits 2 of 3,675
  segments — the panel's paths are near-contiguous, so the cost is negligible.
- **`gbwt_offset` is not portable across graphs.** It is a rank among the sequences visiting that
  node, and a graph containing more sequences gives the same offset to a different path. The
  positions are valid against the `#graph` named in the header and nothing else. `start_node` and
  `end_node`, by contrast, are plain node IDs and travel anywhere those nodes exist.

**`hap_index` is internal; `haplotype` is portable.** `hap_index` is GBWT metadata order over
distinct `(sample, phase)` pairs in *this* graph, so two runs on two chunks of the same genome will
number the same haplotypes differently. The `#haplotype` table binds the two for this file; anything
crossing files must key on the `sample#phase` name. Note also that a haplotype in this sense is a
`(sample, phase)` pair and deliberately collapses the GBWT fragments that make it up — which is
exactly why a segment can carry a haplotype name and still need a specific fragment's position.

**Reference coordinates are advisory, and named.** `#reference` says which path `ref_start`/
`ref_end` are measured against, and it has to: a graph can carry several references — the HPRC
graphs name both `CHM13` and `GRCh38` as reference samples — so a contig name alone does not fix a
coordinate system. The coordinates are a convenience for eyeballing and for range queries. Node IDs
are the anchors, and where node order and reference order disagree — inside the chr20 centromere,
for instance — the reference columns are the ones that mislead.

**What the format does not tell you, and cannot.** Consecutive segments do **not** abut. A run ends
at the last site the outgoing haplotype explains and the next begins at the first site the incoming
one does, so there is an interval between them — 35 bp in the first chr20 example — that belongs to
neither line. That is not an oversight in the encoding: the model observes only that a switch
happened *between two sites*, and has no evidence about where in the intervening sequence the
recombination fell. A consumer reconstructing sequence has to adopt a rule (extend the outgoing
haplotype to the next segment's `start_node` is the obvious one) and should know it is choosing,
not reading.

**`#graph` is a path, not an identity.** It records the graph argument as invoked, so it can be
relative, and it will not detect a rebuilt graph with different node IDs. Reading a mosaic against
the wrong GBZ produces a plausible wrong genome rather than an error. A checksum belongs here and
is not yet implemented.

chr20 at 34 haplotypes gives 7,877 segments over 228,548 sites in **650 KB**, 187 KB gzipped — 84
bytes per segment. The same two walks written out as explicit node lists measure ~40.6 MB — one
haplotype is 2,031,992 steps — so the mosaic is smaller by a factor of about 62.

A whole genome, measured rather than projected: **143,365 segments over 4,742,752 sites in 11.05 MB**,
3.46 MB gzipped, at 80.8 bytes per segment. 99.82% of segments carry a GBWT position and 390 are
fragment splits. Every diploid contig's two strands agree on their site total and the strand-0 total
equals the VCF's record count exactly, which is the property that makes the file a description of a
genome rather than a list of observations.

The trade is that the mosaic is written *by reference* and cannot be read without the GBZ it names,
which is why the header carries the graph. An explicit path list is self-contained and ~137x larger.

**A few segments carry no position.** Genome-wide 254 of 143,365 (0.18%), and on chr20 10 of 3,675,
have `.` in both GBWT columns:
the named haplotype does not visit the segment's start node in either orientation, which happens
where the panel explains a strand through a haplotype that enters the region slightly later. All 10
are in the centromere (26–31 Mb), the same region where node IDs stop tracking reference order.
Such a segment still names its haplotype and its node anchors; only the O(1) entry point is missing,
and a consumer needing it there must fall back to a search.

**Cost.** Resolving positions and splitting on fragment boundaries adds about 3% to a chr20 run
(150 s against 146 s). The obvious implementation — resolving the haplotype's position at every
site, or walking `LF()` across the segment to detect a fragment end — costs 2.3x and 1.2x
respectively. Instead the caller resolves only the two endpoints of a candidate segment and compares
them; if they disagree it binary-searches for the boundary and recurses. Splits are rare (2 on
chr20), so the search almost never runs.

### Haploid contigs

A haploid chain — chrY, or chrX outside the pseudoautosomal regions, in a male sample — is not a
degenerate diploid one. Its state is a single haplotype rather than a pair, so it gets its own
decoding: `H+1` states instead of `(H+1)²`, an ordinary Li–Stephens chain, and the textbook
stay-or-jump max-product step rather than the four-case reduction two coupled strands force.

Set it with `-d 1`, or per contig with `-R 'chrY:1'`. Then:

- `GT` is a single allele. Not `a|a`, which would claim a homozygous diploid call.
- The mosaic carries **strand 0 only**. A second strand would assert two copies where there is one.
- There is no phase to infer on one strand, so what the mosaic gives is purely the ancestry: which
  panel haplotype explains each stretch. On chrY that is the entire answer.

`-d` and `-R` set ploidy per contig, which cannot express a within-contig split — and for chrX that
matters, because the pseudoautosomal regions are diploid and the rest is not. Use `--ploidy-bed` to
give ploidy per *region* and call such a contig in one pass. Linkage and the mosaic still break at
each ploidy boundary, because a chain is a maximal run of one ploidy, but the break now falls where
the biology does instead of where a run was spliced.

`--phased` and `--mosaic-out` both need `--linkage-weight` above 0. `--mosaic-out` names a path, so it
is always an explicit request and always an error with the layer off. Phasing is on by default, so it
follows the layer: where there is no haplotype panel to decode against, the layer declines and phasing
declines with it, and the run comes out unphased rather than refusing to start. Asking for `--phased`
outright against a disabled layer is still an error.

### A haploid record on a diploid contig is not the same thing

Under nested calling -- on by default with `--read-likelihood` -- a site inside a chain that only one
parent allele crosses is genotyped at ploidy 1 --
the other haplotype has no sequence there to call, because the parent's other allele deletes the chain.
That is one strand of a diploid locus, not a haploid locus, and the two must not be written alike:

- A haploid *contig* gets a single allele, as above. There is no other strand.
- A haploid *nested site* gets `a|.` or `.|a`. The position of the allele in the pair is which strand
  of the parent it sits on, and `.` marks the strand that carries nothing.

Without the pair there is nowhere in the record to say which strand, and the assignment existed only in
the mosaic -- where no phasing tool would look for it, and where whatshap could not have read it anyway,
since it refuses a file of mixed ploidy. `.` rather than `*`: `*` would say "absent because something
deleted it", which is more precise, but it is an ALT allele and adding one changes the arity of `AD`,
`GL` and `GQI`, which are written long before the strand is known.

A nested site whose parent's final genotype contradicts its ploidy never acquires that ploidy in the
output, because the record does not exist yet when the contradiction is resolved. Every site is
genotyped during the read sweep and *staged*; the barrier then settles each generation's genotypes,
handing each chain the ploidy its settled parent implies; and only then is any record built. A chain
the settled genotype crosses on neither haplotype is dropped before it becomes a record at all, since
the sample has no copy of it for a call to sit on.

There are consequently no `nested_*` FILTERs, and there is nothing for them to flag. They were removed
along with the machinery that made them necessary: when a record was written before its genotype was
final, the only way to change it afterwards was to patch the line, and a patch can neither add an ALT
nor withdraw a line -- so a settled genotype naming an allele the line did not carry, or settling on
the reference, left a record that had to be marked rather than fixed. Deciding before rendering
removes the class instead of labelling it. On chr20 both populations are now zero by construction
rather than by inspection: no settled genotype is unrenderable, and no record carries a hom-ref
genotype.

The same reordering is why the genotype in the VCF is always the one the model settled on. There is no
second pass that revises it: the allele list, the symbolic-reference test that decides whether a line
is written at all, `QUAL`, and the arity of `AD`/`GL`/`GQI` are all derived from the settled genotype
when the record is built.

## Genotype quality and the VCF fields

| Field | Meaning |
|---|---|
| `GT` | The called genotype: the argmax of the objective, or of the linkage posterior when that layer is on. |
| `GL` | `log10 P(reads \| genotype)` for every genotype, in VCF order. |
| `GQ` | Phred gap between best and second-best genotype, **scaled by the fraction of reads the call explains**. |
| `GQI` | The same gap *without* that scaling. Equals `GQ` under `--no-share-quality`. |
| `GQN` | The same gap as a **fraction of what this site could have achieved**, in `[0,1]`. Unlike `GQ` it means the same thing at any depth and any ploidy — see below. |
| `GP` | **Natural**-log-scaled posterior of the called genotype under a uniform prior. Nats, not log10: unlike `GL`, exponentiate with `e`, so `-2.303` means `p = 0.1`. |
| `DP` | Reads overlapping the site. |
| `AD` | Reads whose best-fitting allele is this one. **Does not sum to `DP`** — see below. |
| `DR` | `N_eff / λ_G` at the called genotype — see [above](#dr-the-depth-ratio). `1.0` means the read count matches the call. Absent where the ratio is undefined. |
| `BL` | Mean over reads of the best raw alignment score any allele gave them. |

**Why `GQ` is scaled.** A read whose best-fitting allele lies *outside* the called genotype enters
every genotype's likelihood and therefore cancels out of the best-versus-second-best comparison.
The likelihood ratio simply cannot see it. Scaling by the explained share restores it. The
consequence is that `GQ` here is a *quality score for ranking*, not a calibrated posterior;
`GQI` is the pure likelihood ratio for anyone who wants it.

### `GQN`, and why `GQ` is not comparable between sites

`GQ` answers "how much better is the winner than the runner-up", in absolute log-likelihood. That
quantity is not comparable across sites, and it moves on two axes.

**Depth.** The gap is a sum over reads, so it grows with coverage. Paired on identical sites across
a coverage titration, the gap per read is close to constant (3.2 at 5x to 3.8 at 30x), so this axis
is nearly linear — 30 reads give about twice the gap of 15.

**Ploidy, which is the larger axis and the surprising one.** At ploidy 1 the runner-up genotype is a
different allele outright, so a read that fits the call gives the runner-up nothing but the mismap
floor. At ploidy 2 a heterozygote's runner-up differs on one strand only, so a read discriminates
about half as much — *and only half the reads discriminate at all*, since reads from the shared
allele fit both genotypes. Per read, at the default `--mismap-min 0.02`:

| called | runner-up | achievable gap per read |
|---|---|---|
| `{A}` haploid | `{B}` | `−ln(0.02)` = **3.91** nats |
| `{A,B}` diploid het | `{A,A}` | `½·(ln 0.51 − ln 1) + ½·(ln 0.51 − ln 0.02)` = **1.28** |
| `{A,A}` diploid hom | `{A,B}` | `−ln(0.51)` = **0.67** |

So a hemizygous call carries roughly six times the `GQ` of a diploid homozygote on identical
evidence. Measured on HG002, chrX hemizygous calls run a median `GQ` of 247 where chr7 diploid
homozygotes at the same depth run 46.

`GQN` divides the observed gap by the achievable one, giving *what fraction of the discrimination
this site could offer did the data actually deliver*. It carries the same explained-share discount
as `GQ`, for the same reason.

Two details that are easy to get wrong, both of which were got wrong first and fixed by measurement
over a coverage titration of chr20 (diploid, 5–30x) and chrX (haploid, 2.5–14.6x). The figure quoted
is the mean spread of observed precision at a claimed score across those arms, lower being better:

- **The denominator is built at the mismap floor, not at the reads' own `e_r`.** An ideal read is
  well-fitting *and* well-mapped. Using each read's own `e_r` folds the site's unreliability into
  both sides of the ratio, where it cancels: a window of MAPQ-0 reads offers `−ln(0.7)` = 0.36 per
  read, so a badly mapped site gets a denominator small enough to make a weak call look strong. That
  version scored **0.427** — worse than not normalising at all.
- **Both terms of the heterozygote's difference are kept.** Reads from the allele the runner-up also
  carries actively favour the runner-up, and dropping them overstates a het's achievable gap by
  about a quarter (1.62 against 1.28).

| | mean spread |
|---|---|
| raw `GQ` | 0.347 |
| `GQ / DP` | 0.494 |
| `GQN` with per-read `e_r` | 0.427 |
| `GQN` at the floor, no explained-share | 0.316 |
| **`GQN` as shipped** | **0.259** |

**`GQ/DP` is the tempting fix and it is worse than doing nothing.** It halves the spread on the
diploid series alone (0.101 to 0.050), which is why it has to be checked across ploidies: it
corrects the smaller axis, leaves the larger one, and compresses the score range so that what
remains does more damage.

`GQN` is emitted as `.` where there was no gap to normalise — no reads, or a site offering one
genotype. That is not the same as `0` and should not be filtered as though it were.

### Using `GQN` — practical guide

**The one rule to take away: do not threshold `GQ` on low-coverage data.** Requiring `GQ ≥ 10`
takes a 5x diploid contig's F1 from **0.888 to 0.669**, because at that depth roughly half the true
calls sit below it. `GQ` is an absolute likelihood ratio, so its scale moves with both depth and
ploidy. `GQN` is the same evidence expressed as a fraction of what the site could have delivered,
so one number means the same thing everywhere.

Filter with either the caller or `bcftools` — they are equivalent, and the flag only exists so you
do not have to post-process:

```bash
vg call ... --read-likelihood --min-confidence 0.05        # marks as FILTER=lowconf
bcftools view -i 'FMT/GQN>=0.05' calls.vcf.gz              # same threshold, after the fact
bcftools view -f PASS calls.vcf.gz                         # drop what --min-confidence marked
```

`--min-confidence` **marks, never drops**, so nothing is lost by using it and a downstream
`bcftools view -f PASS` completes the job when you want the records gone.

#### What a threshold buys, measured

At `GQN ≥ 0.05`, over a coverage titration of chr20 (diploid) and chrX non-PAR (haploid):

| input | precision | recall | F1 |
|---|---|---|---|
| 5x diploid | 0.9712 → **0.9759** | 0.8174 → 0.7993 | 0.8877 → 0.8788 |
| 30x diploid | 0.9828 → **0.9863** | 0.9230 → 0.9139 | 0.9520 → 0.9487 |
| 2.5x haploid | 0.8863 → **0.9256** | 0.8047 → 0.7936 | 0.8435 → **0.8545** |
| 14.6x haploid | 0.9389 → **0.9840** | 0.9280 → 0.9184 | 0.9334 → **0.9501** |

Precision rises on **every** arm for one to two points of recall. F1 improves on the haploid arms
and falls slightly on the diploid ones.

#### Choosing a threshold

**There is no good default, and that is a finding rather than an omission.** F1 weights precision
and recall equally; under that weighting a threshold helps haploid contigs and hurts diploid ones,
so no single value wins everywhere. Which error you would rather make is not something the caller
can know. Hence:

| your situation | suggestion |
|---|---|
| diploid autosomes, any coverage | **no threshold**. Low coverage already costs recall; a filter costs more. |
| haploid contigs (male chrX/chrY) | **`--min-confidence 0.05`**. Low coverage costs *precision* there, and this is worth +0.017 F1 at 14.6x. |
| precision matters more than recall | **`--min-confidence 0.05`**, any ploidy. Raises precision everywhere for 1–2% recall. |
| you were going to filter on `GQ` | use `GQN` instead — see the first paragraph. |

Values above ~0.15 buy little further precision and cost recall steadily; 0.05 is where the curve
turns on every arm measured.

#### Reading a value

`GQN` is *the fraction of the discrimination this site could offer that the data actually
delivered*. `0.8` means the reads did 80% of what a perfect pile-up would have done here; `0.02`
means the site had the potential to be decisive and the reads were not. It is bounded in `[0,1]`
by construction, so a value near 1 is as good as this site gets — not as good as some other site
gets, which is exactly the comparison `GQ` cannot support.

Two things it is **not**. It is not a calibrated probability: `GQN = 0.5` does not mean a 50%
chance of being right. And it is not comparable to `GQ` numerically — they are different scales,
and `GQI` remains the raw likelihood ratio for anyone who wants it.

**Why `AD` does not sum to `DP`.** Two reasons: a read fitting several alleles equally splits its
vote, and — much more importantly — only alleles that reached the VCF record get a column, while
the genotyper scored every allele the site offered. Where many alleles were enumerated and few
emitted, most reads best-fit something absent from the record. That shortfall is itself
informative: it is how much of the evidence the emitted alleles fail to explain.

### Ploidy by region

`-d` and `--ploidy-regex` set ploidy per *contig*, which cannot express the case that actually
arises. A male sample's chrX is haploid **except** in the pseudoautosomal regions, where X and Y
recombine and two copies are present. `--ploidy-bed FILE` takes a BED of `CHROM START END PLOIDY`
and overrides the contig ploidy wherever an interval covers a site:

```
chrX	0	2394370	2
chrX	153926003	154259566	2
```

`CHROM` matches the contig as the **output VCF** spells it (`chrX`, not `CHM13#0#chrX`). Intervals
are BED half-open and 0-based, so `[0, E)` covers VCF POS `1..E`. Anything no interval covers keeps
the ploidy from `-d`/`--ploidy-regex`. Overlapping intervals are rejected rather than resolved by
precedence — a BED saying two things about one base has no correct reading.

**Why getting ploidy right matters more than it sounds.** Under ploidy 2 a site whose reads split
evenly between two alleles is a heterozygote — an ordinary genotype, and at worst half wrong. Under
ploidy 1 there is no such genotype, so the model must pick one allele on almost no evidence and
about half the time picks the wrong one, with a likelihood gap near zero. Anything producing a
balanced pile-up on a haploid contig therefore becomes a confident-looking false positive rather
than a merely half-wrong call. Two diverged paralogous copies collapsed onto one locus do exactly
that: on HG002's chrX, two ~200 kb such loci produce **29% of the whole chromosome's false
positives**, and their false calls have a median minor-allele fraction of 0.375.

A ploidy BED is the fix where the region genuinely is diploid. Where it is not — where the extra
reads come from a paralog somewhere else entirely — no ploidy setting helps, and the reads' balance
across `AD` against the called alleles is the signal that says so.

**Haploid linkage reaches the VCF.** It did not always: a haploid record's `GT` is a bare allele,
and the code that patched genotypes after the linkage pass built the genotype it expected as
`"i/j"` regardless, so its guard rejected every haploid change. The layer ran, reported the changes
in the progress line, and discarded them -- chrY and non-pseudoautosomal chrX received no linkage
correction at all, and because phasing and the mosaic are built from the *post*-linkage genotypes,
the mosaic described genotypes the VCF did not contain. Measured after the fix, haploid linkage is
worth **+0.017 F1** at 2.5x and **+0.008** at 14.6x on chrX.

**Linkage and the mosaic break at every ploidy boundary.** A chain is a maximal run of one ploidy
on one contig, not a whole contig. That is not a limitation to work around: the transition model
moves probability between adjacent sites through *pairs* of panel haplotypes, and across a ploidy
change there is no correspondence to carry, exactly as there is none between two contigs. Each side
gets its own phase set.

**What this replaces.** Before it, calling a male chrX meant running the contig twice at different
ploidies and splicing the VCFs on the PAR boundaries. Against that two-pass output, one pass with
the BED above reproduces all **97,068** non-PAR sites genotype for genotype, and switches ploidy
exactly at the boundaries with nothing leaking either way.

The splice was also quietly wrong. `bcftools view -t "^chrX:153926003-"` does not exclude an
open-ended range the way `-r` includes one, so 190 haploid records leaked into PAR2 and the
concatenated VCF carried 190 duplicated positions at contradicting ploidies. Nothing published was
affected — the T2T-Q100 confident regions end at 153,910,814, before PAR2 begins, so those records
were never scored — but it is the kind of error a splice invites and a ploidy BED cannot make.

### `--min-confidence`, and why there is no default

`GQN` exists so that one threshold can mean the same thing everywhere, and `--min-confidence X`
is that threshold: records below it are marked `FILTER=lowconf`. A site with no informative read at
all is marked `FILTER=noreads` instead -- the model has no evidence either way there, which is not
the same as low confidence. Both mark and never drop. A record whose genotype the linkage layer
moved off the per-site argmax has its `GQN` blanked to `.` and a stale `lowconf` cleared, because both
were measured on the genotype the reads alone preferred, which the record no longer carries. Its `GQ`
becomes the phred complement of the linkage posterior, discounted by the explained-read share and
capped at the per-site `GQI` -- the layer may lower confidence and may not raise it above what the
reads on their own supported.

A raw `GQ` threshold cannot do this job, and the failure is not subtle. Requiring `GQ ≥ 10` takes
a 5x diploid contig's F1 from **0.888 to 0.669**, because at that depth half the true calls sit
below it. `GQN ≥ 0.05` costs the same arm 0.009.

At `GQN ≥ 0.05`, precision rises on **every** arm measured:

| arm | precision → | recall → | F1 → |
|---|---|---|---|
| chr20 5x (diploid) | 0.9712 → 0.9759 | 0.8174 → 0.7993 | 0.8877 → 0.8788 |
| chr20 30x (diploid) | 0.9828 → 0.9863 | 0.9230 → 0.9139 | 0.9520 → 0.9487 |
| chrX 2.5x (haploid) | 0.8863 → 0.9256 | 0.8047 → 0.7936 | 0.8435 → 0.8545 |
| chrX 14.6x (haploid) | 0.9389 → 0.9840 | 0.9280 → 0.9184 | 0.9334 → **0.9501** |

**There is no good default, and that is a finding rather than an omission.** F1 weights precision
and recall equally, and under that weighting a threshold helps the haploid arms and hurts the
diploid ones — no single value wins everywhere. What it is worth depends on which error you would
rather make, which this code cannot know for you. So it ships off.

**It marks rather than drops**, and that is a correctness decision, not a convenience. Records are
buffered as text and the linkage layer rewrites their genotypes afterwards; a record withheld at
emission is withheld before linkage ever sees it, and a low-confidence site is exactly the kind
linkage exists to fix. `bcftools view -f PASS` is one step away for anyone who wants the records
gone, and cannot be undone by anyone who does not.

**`--depth-quality A`** scales `GQ` by `exp(−A · |ln DR|)` at records whose called alleles change
length by at least 50 bp, so a call whose read count is implausible for the sequence it claims
ranks lower. Off by default, and affects ranking only. See [`DR`, the depth
ratio](#dr-the-depth-ratio) for how to read the value it acts on.

## `--atomize-blocks`: one snarl, several variants

Symbolic collapsing answers one question per haplotype — is this the same *route* through the snarl
as the reference? A traversal that answers yes becomes allele 0 and its differences descend to the
child chains that own them. A traversal that answers no becomes a single ALT spanning the whole
snarl, whatever the shape of the difference.

`--atomize-blocks` — on by default under `--read-likelihood` — replaces that bit with an alignment. The reference and each called haplotype are
projected to symbolic alleles over the alphabet of *plain nodes plus one symbol per child chain*,
aligned, and the maximal runs of non-matching steps become difference blocks. Each block is one VCF
record.

The cost model is edit distance **with substitution at cost 1**, not the insert/delete-only model a
plain `diff` uses, and the reason is that "minimum edits" does not determine the answer without it.
Take σ(R) = `[a,b]` against σ(H) = `[b,b]`: matching `b` to the first `b` gives *delete a, match,
insert b* — two blocks separated by a spurious match — and matching it to the second gives one
substitution. Both cost 2 under insert/delete-only, so nothing chooses between them, and the choice
decides how many records the snarl emits. Substitution at 1 beats delete-plus-insert at 2, so the
one-block reading wins strictly. Remaining ties are broken by preferring the diagonal, then
deletion, then insertion, walked backwards from the end — an arbitrary rule, but a fixed one, which
is what output reproducibility needs.

### Reporting each difference exactly once

A chain symbol can appear *inside* a difference block, where the block's ALT spells out the route
through it — and descent would otherwise emit a record for that chain as well, reporting the same
variation twice. The rule that prevents it:

> a chain's ploidy at this level counts only its crossings inside **matched** steps

so a chain every called haplotype crosses only inside a block is not descended into. In practice
this population is tiny — 47 of chr20's 115,996 called ALTs — because the two larger cases are
already handled elsewhere: a chain the reference does not cross is never descended into at all
(there is no reference path to call it against), and a genotype carrying the reference allele matches
every chain by definition.

### What every block of a snarl shares, and why that matters

`AD`, `GL`, `GQ`, `GQI`, `GQN`, `GP`, `DP`, `DR` and `QUAL` are **the site's values, carried onto
every block**. Arity is correct — each block allele inherits the evidence of the site allele its
haplotype carries, so `AD` has one entry per allele and `GL` one per genotype — but the *values* are
replicated, because the likelihood was computed over whole-snarl traversals and there is no
principled per-block decomposition of it. A consumer that sums or averages evidence across records
will therefore double-count a split snarl.

`INFO/SB` gives block index and block count so that set is recoverable. The VCF ID stays the snarl
name, unsuffixed, and is shared by every block of a snarl: several things in `vg` parse the ID back
to a snarl and two of them abort rather than degrade on a name they do not recognise.

**This caveat is why the limitation is stated in `INFO/SB`'s own description rather than only here.**
Records carrying `SB` came from a decomposition and share an ID with their siblings; records without
it are the only record their snarl emitted and are unaffected. A consumer that aggregates evidence
must group by ID and count each snarl once. Per-block `AD` and `GL`, by a marginal fold over
`genotype_lls`, is the work that would retire the caveat, and it is not done.

### Measured effect

chr20, 32-haplotype graph, `readlik` arm, both sides called by the same binary in the same session
and scored under fresh labels. 487 snarls decompose, into 1,134 records where they previously
produced 487; 362 child chains are suppressed as already-spelled; net **115,038 → 115,427 records**.

| | off | on | delta |
|---|---|---|---|
| aardvark recall, all types | 0.965984 | 0.966206 | **+0.000222** |
| aardvark F1, all types | 0.972250 | 0.972322 | +0.000072 |
| aardvark F1, indels | 0.927672 | 0.928057 | **+0.000385** |
| aardvark F1, SNVs | 0.985193 | 0.985156 | −0.000037 |
| basepair F1, all types | 0.924542 | 0.924903 | **+0.000361** |
| basepair bases matched | 378,392 | 378,486 | **+94** |
| basepair bases claimed | 427,914 | 427,752 | **−162** |

Recall rises in every class. The caller claims fewer bases while matching more of the truth, which
is the effect the flag exists to produce, stated as directly as this evaluation can state it.

**Those are all-types figures, and they dilute the effect into noise.** Structural variants are about
765 of 94,691 truth variants here, so a real structural gain vanishes in the aggregate. Measured
directly on variants ≥50 bp, with truvari, both arms from the same binary:

| | off | on | delta |
|---|---|---|---|
| true positives | 433 | 444 | **+11** |
| false positives | 444 | 444 | **+0** |
| false negatives | 332 | 321 | **−11** |
| recall | 0.5660 | 0.5804 | **+0.0144** |
| **F1** | **0.5247** | **0.5346** | **+0.0099** |

Eleven more true structural variants with no additional false positives.

**chr6 gives a smaller and differently-shaped result, and the genome-wide figure is smaller still.**
On chr6-34hap — 1,547 truth SVs against chr20's 765, so the better-powered test — the same comparison
gives F1 0.582043 → 0.583696, **+0.0017**, driven by precision (14 fewer false positives) with recall
marginally *down* (−2 true positives).

Genome-wide on HG002, 23 scoreable contigs: **SV F1 0.5577 → 0.5620, +0.0043**, from 48 more true
structural variants and 266 fewer false ones — so both sides improve, unlike either single contig.
Autosomes alone give 0.5596 → 0.5642. That aggregate is the number to quote; chr20's +0.0099 is the
favourable tail of a 6× per-contig spread, not the typical case.

Under `truvari refine` the two agree far better — **+0.0233 on chr20 and +0.0121 on chr6** — and
chr6's ladder shows why. Loosening the distance tolerance alone does nothing there (+0.0017 at
refdist 500, −0.0002 at 1000, −0.0007 at 2000) while re-alignment gains +0.0121. That is the
signature of decomposed calls being the same events spelled differently: distance tolerance cannot
recognise them and sequence re-alignment can. So the record-matching metric systematically
understates this change, and understates it more where there are more records.

This figure is deliberately the *unrefined* truvari number, because that is what every other SV
figure in this project uses. As a diagnostic, `truvari refine` puts the same comparison at 0.6154 →
0.6387, **+0.0233**: the gain grows under re-alignment rather than shrinking, with false positives
falling 362 → 349. So the unrefined number is a floor, not an estimate — which is what you would
expect, since record-level matching penalises a change that splits records, and the 50 bp size floor
discards any block that falls below it.

So the honest summary is: **below resolution on small variants, and a real gain on structural
variants** — the population the design was aimed at. 487 snarls of 115,038 records decompose, 0.42%,
and their effect is concentrated where alleles are long.

Two things bound the size of the effect, and both are properties of the graph rather than of this
code. A single difference block spans every step but the two snarl boundaries, so an ALT that does
not split is barely changed — and `flatten_common_allele_ends` already trims those ends. And 9,279
sites on this contig had no symbolic projection at all, because `flip_snarl` reverses a snarl whose
reference path runs backwards and the reversed copy failed the identity check that finds the snarl
in the manager; symbolic collapsing was inert there, and so was this. Fixed separately — the check
now accepts the boundary pair in either order, since a snarl and its reversal are the same snarl.

Worth being precise about what that defect did and did not break, because the symbolic form itself
was never wrong. Every `SymbolicStep` carries its own orientation, and a chain symbol's identity is
the chain's boundary node ids, which are orientation-free. What failed was one guard asking "is this
the snarl I think it is", answered by comparing boundaries in a fixed order — so a snarl viewed from
the other end was rejected, no child was recognised, and the projection came out as the plain
oriented node path. Orientation intact; nothing symbolic in it. Descent was unaffected, because its
own lookup carries no such check, so those sites emitted child records *and* the parent's redundant
long ALT. Repairing it therefore removes a duplicate rather than recovering a missing call.

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

The flags introduced by this mode (the read sources, the model terms, the linkage and phasing
options, `--mosaic-out`, `--dump-likelihoods`, `--enumerate-support`, `--min-confidence`) require
`--read-likelihood`; `vg call` errors if one is passed without it, including in getopt's abbreviated
spellings. General options that this mode also uses -- `-d`/`--ploidy`, `-R`/`--ploidy-regex`,
`--ploidy-bed`, `--nested`/`--no-nested` -- work on the default caller too and are not gated.

### Reads in — exactly one source required

| Flag | Default | Effect |
|---|---|---|
| `--gam FILE` | — | Read alignments, held in memory. |
| `--gaf-reads FILE` | — | The same, as GAF. |
| `--gam-index FILE` | — | `.gai` for `--gam` (`vg gamsort -i`), so reads are fetched per site instead of all held in memory. |
| `--gaf-base FILE` | — | GAF-Base database, queried per site via `gbz-base`. |
| `--gbz-base FILE` | input graph | Graph to resolve `--gaf-base` queries against. |
| `--gaf-base-binary FILE` | `gbz-base` on `PATH` | The `gbz-base` executable to spawn for `--gaf-base` queries. |
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
| `--flat-mixture` | off | `w_h = 1/ploidy` instead of length-weighted. Restores the pre-correction mixture; the depth term keeps its read-length estimate and is unaffected. |

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
| `--min-confidence X` | 0 (off) | Mark records with `GQN < X` as `FILTER=lowconf`. 0.05 raises precision on every arm measured. Marks, never drops. No default: it helps haploid F1 and hurts diploid. |

### Ploidy — does change genotypes

| Flag | Default | Effect |
|---|---|---|
| `-d/--ploidy N` | 2 | Ploidy of the sample, `{1, 2}`. |
| `-R/--ploidy-regex RULES` | — | Comma-separated `REGEX:PLOIDY` rules assigning ploidy to contigs by name. Unmatched contigs fall back to `-d`. |
| `--ploidy-bed FILE` | — | BED of `CHROM START END PLOIDY` setting ploidy per *region*, overriding `-d`/`-R`. Linkage and the mosaic break at each boundary. |

### Phasing and mosaic output — needs `--linkage-weight` above 0

| Flag | Default | Effect |
|---|---|---|
| `--phased` | **on** | Emit phased genotypes (`0\|1`) and a `FORMAT/PS` phase set from the linkage decoding. Declines with the linkage layer when there is no panel to decode against; asking for it *explicitly* with the layer off is an error, not a silently unphased file. |
| `--no-phased` | — | Turn phasing off: unphased genotypes and no `FORMAT/PS`. |
| `--nested` | **on** under `--read-likelihood` | Symbolic collapsing and ploidy-propagating descent, so a variant inside a child chain gets its own record instead of being buried in a long ALT. Genome-wide it takes SNV F1 from 0.9752 to 0.9833 and SV F1 from 0.5134 to 0.5467 at no runtime or memory cost. Declines on the support-based caller, where it has never been measured; an explicit `--nested` still works there. |
| `--no-nested` | — | Genotype each snarl against its own full traversals, with no collapsing and no descent. |
| `--atomize-blocks` | **on** under `--read-likelihood` | Align the reference and each called haplotype as *symbolic* alleles and emit one record per difference block, so a snarl differing from the reference in two separated places reports two variants instead of one substitution spanning both. Worth +0.0043 SV F1 genome-wide (0.5577 → 0.5620 over 23 contigs; higher under `truvari refine`, and the per-contig spread is wide); small-variant F1 unchanged. Declines on any other calling path and with `-a`, whose record set must stay sample-independent; asking for it there by name is an error. Every block of a snarl repeats the site's `AD`/`GL`/`GQ` — see `INFO/SB` and [below](#--atomize-blocks-one-snarl-several-variants). |
| `--no-atomize-blocks` | — | One record per snarl, whatever the shape of the difference. |
| `--mosaic-out FILE` | — | Write the inferred genome as a run-length-encoded mosaic of panel haplotypes. Implies `--phased`. Format [above](#--mosaic-out-file). |

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

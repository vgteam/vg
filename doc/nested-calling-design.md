# Nested calling by symbolic alleles, with ploidy propagation

A design and implementation plan. Nothing here is implemented yet; the measurements motivating it
are, and are cited so each claim can be re-checked.

## Why

`vg call` emits one record per top-level snarl, and a `SnarlTraversal` runs from snarl start to snarl
end through every interior node. Nested variation is therefore already baked into each top-level
allele string, and two traversals differing only *inside* a nested chain become two long, nearly
identical ALTs. Measured on HG002 against T2T-Q100, autosomes:

| | |
|---|---|
| SNVs vg misses that sit inside a large allele vg itself emitted | **55,222** of 142,707 (38.7%) |
| the same, among variants PanGenie called | 71.2%, against a **0.6%** rate among variants vg calls correctly |
| indels likewise | 10,735 of 78,274 (19.6%) |
| vg SV "false positives" that are same-length substitutions | 2,219, against PanGenie's 29 |
| of those, differing at <=10 bases | **90.6%** — including a 4,710 bp allele differing at 3 |

One record can swallow many variants: across chr1+chr6+chr20, 349 large records account for 8,168
missed small variants, a mean of 23.4 each and a maximum of 435. The median swallowing record is
4,420 bp of REF; the largest is 210,215 bp.

The failure is not that these variants are unreported. It is **all-or-nothing at snarl scope**:
choosing one traversal commits to every nested variant inside it at once, so a slightly wrong long
allele loses all of them together. The 0.6% control rate proves the scorer credits nested variants
when the enclosing allele is right — aardvark compares by local haplotype — so the 71.2% are cases
where the allele is wrong.

Two further facts bound the work. **39.9% of vg's SNV false negatives are not in the panel at all**
and are irreducible. Of the rest, recovering the swallowed set would move SNV recall from 0.9570 to
about 0.9737, past PanGenie's 0.9659; the graph's own ceiling is 0.9828.

Full workings: `docs/small-swallowed.md`, `docs/fn-representable.md`, `docs/sv-nocall.md` and
`docs/sv-delta.md` in the companion evaluation repository.

## The idea

A **leaf snarl** contains no nested snarls; a **non-leaf snarl** contains nested chains.

For a non-leaf snarl, rewrite each traversal as a **symbolic allele**: the ordered sequence of steps
where every maximal excursion through a nested chain is replaced by a single symbol naming that
chain, rather than the concrete path taken through it. Allele identity at this snarl is then decided
by comparing symbolic sequences.

Two consequences follow:

- **A traversal symbolically equal to the reference traversal is the reference allele here.** Its
  difference lies entirely inside nested chains, and belongs to those chains' own records. The
  4,710 bp allele differing at 3 bases stops being a top-level substitution; the three SNVs become
  three nested calls.
- **A traversal that skips a chain, or crosses different chains, is symbolically distinct** and stays
  a genuine top-level call. A real deletion is still reported as a deletion.

The second half is **ploidy propagation**. Once the snarl is genotyped, each child chain's ploidy is
the number of called parent alleles that cross it: both alleles cross it, ploidy 2; one crosses and
the other deletes it, ploidy 1, and the chain is called haploid; neither crosses it, ploidy 0, no
descent, star allele. Recursion continues down the hierarchy with that ploidy.

**The read-likelihood model does not change.** What changes is how traversals are collapsed into
alleles for emission, and which snarls are visited at what ploidy.

## What already exists

This is not a new representation. The decomposition provides it:

| piece | where |
|---|---|
| a `Visit` carrying *either* a `node_id` *or* a `Snarl` — "each step is given as either a node or a child Snarl" | `deps/libvgio/deps/vg.proto:283` |
| `SnarlManager::chains_of(parent)` — the child chains to symbolise | `src/snarls.hpp:513` |
| `NetGraph`, in which "each chain and unary child snarl is treated as an ordinary node" — the symbolic view, already built | `src/snarls.hpp:245`, `net_graph_of` at `:521` |
| `into_which_snarl(visit)` — maps a visit to the child snarl it enters | `src/snarls.cpp:923` |
| a loop already projecting a traversal's visits onto child snarls | `src/graph_caller.cpp:3071` |

The deprecated `NestedFlowCaller` used Snarl-carrying Visits already — it had the representation and
lacked the ploidy discipline. Its Visit form is also why `--bottom-up -T` aborts in the GAF emitters.

### The defect this also fixes

Recursion today is a side-effect of emission. `FlowCaller::call_snarl_internal` ends with
`ret_val = trav_genotype.size() == ploidy && added` (`graph_caller.cpp:3749`), and `emit_variant`
returns `true` even when it emitted nothing, because the `genotype_snarls || !alt.empty()` gate at
`:1893` is simply skipped for a hom-ref call and control falls through to `return true` at `:1976`.
So **a snarl genotyped hom-ref writes no record, reports success, and its children are never
queued.** `RecurseOnFail` only fires where the caller could not decide at all — in practice, sites
with no reads. Under this design, recursion becomes an explicit decision driven by propagated
ploidy, not a by-product of whether a line was written.

## Decisions taken

**Arg-max before marginalisation.** Several concrete traversals can share one symbolic allele. Under
arg-max, the model runs exactly as it does now over concrete traversals, and collapsing to symbolic
happens strictly afterwards, for emission and propagation. Every quantity stays defined on a single
concrete traversal pair: the mixture weights `w_h`, the expected count `λ_G`, and hence
`DR = N_eff/λ_G`.

Marginalising over a symbolic class breaks that at the **depth term**. Class members differ inside
nested chains, which may contain indels, so they differ in length and therefore in `λ`. A logsumexp
over genotypes carrying different `λ` is a sound mixture but leaves `λ_G` undefined for the reported
call, so `DR`, the depth term and any gate on them have no single value. Marginalisation is
therefore a later stage with its own answer for `λ` — candidates are the arg-max member's `λ`, a
length-weighted average over the class, or a `λ` recomputed from the symbolic allele's
reference-projected length.

Note this is only lossy where it matters least. Symbolically equal traversals cross the same chains,
so ploidy propagation is identical whichever member arg-max picks; and where the call is symbolically
hom-ref nothing is emitted, so the parent's GQ is unused. The case marginalisation would improve is
the het one — a real deletion competing against a *class* of symbolically-reference traversals whose
probability is currently split across its members.

**Propagated ploidy is already safe for the depth model.** `local_read_rate` is reads/bp over the
read source's fetch window and is deliberately not divided by ploidy; the division happens at point
of use. A nested site called at ploidy 1 inside a heterozygous parent therefore gets the right `λ`
without further work.

**Off-reference nested content: opt-in, off by default.** A chain crossed only by a non-reference
parent allele has no reference path through it, so REF and POS for a nested record are ill-defined.
v1 descends only into chains the reference traversal also crosses — which is where the swallowed
variants are, since they are on-reference by construction. `--nested-pseudo-ref` enables descent
elsewhere against a pseudo-reference. Off-reference calling proper is intended later work.

**Repeated traversal of a chain by one allele** (cycles, tandem duplication) gives per-haplotype copy
number above 1. v1 caps at the {1, 2} the rest of the system assumes and logs the occurrence.

**Nested sites in the linkage layer and the mosaic: measured, not assumed.** Including them is more
correct for phasing but changes chain construction and the mosaic's site accounting, which the
harness asserts must equal the record count. Stage 5 runs it both ways and decides on the numbers.

**Shared-flank trimming: a separate, opt-in option.** `--trim-shared-flanks` removes the common
prefix and suffix of a REF/ALT pair before emission, so a record reports the edit rather than the
snarl. This is *not* the atomisation already tried and rejected (`--atomize-substitutions`, which
split a same-length record into one record per differing run and bought SV F1 0.4998 → 0.5015 on
chr20-4hap — real but small). Trimming is complementary and aimed at a different population:
a genuine 300 bp deletion in a snarl whose boundary anchors add shared flank is emitted at 400 bp,
and truvari's size-similarity test sees the wrong size. The specific hypothesis to test is that
trimming recovers part of the **561 vg-only false negatives that had a comparable call nearby and
were rejected on sequence or size similarity** (`docs/sv-unmatched.md`). Symbolic collapsing is
expected to subsume most of the substitution population on its own, so trimming must be measured
*after* it, not before, or it will be credited with work symbolic alleles did.

## Stages

Each stage has a gate. A stage that misses its gate stops the sequence rather than being carried
forward on the assumption a later one will rescue it.

### Stage 0 — Offline validation, no C++

Establish that the swallowing records are non-leaf snarls with child chains. If they are leaf snarls
holding one large complex bubble, symbolic collapsing buys nothing there and this whole design is
aimed at the wrong population.

Method, all from existing artefacts: take the snarl hierarchy from `vg snarls` on chr20, and the `AT`
INFO field of each emitted record, which already gives each allele's traversal as a node path.
Project each `AT` onto the child chains and compare symbolic sequences.

**Gate**: the records predicted to collapse to reference must account for the bulk of chr20's share
of the 2,219 same-length substitution false positives and of the 55,222 swallowed SNVs. Report the
fraction of swallowing records that are non-leaf.

### Stage 1 — Symbolic allele encoding

`symbolic_allele(trav, snarl, snarl_manager)` returning a symbol sequence, with equality and hashing.
Symbol granularity is the **child chain**, not the child snarl, and symbols carry orientation.

Unit tests on synthetic hierarchies: chain present, chain absent (deletion), chain crossed twice,
reversed orientation, trivial chains, and a nested chain inside a nested chain. Also an early probe
that panel enumeration behaves at a nested snarl under a propagated ploidy, since that is a
precondition for Stage 3 and is cheaper to find out now.

### Stage 2 — Reference-equivalence collapsing at emission

Map called traversals to symbolic alleles at emission; a traversal symbolically equal to the
reference collapses to allele 0, and only symbolically distinct ALTs are emitted.

**Gate**: chr20's same-length substitution false positives largely disappear; small-variant F1 does
not fall.

**Measured: the first half holds and the second does not, and the reason matters.** On chr20, SV
false positives fall 367 -> 282 and SV F1 rises 0.4944 -> 0.5106; small-variant F1 falls 0.9646 ->
0.9596, losing 912 variants, 585 of them (64%) inside the 660 records the collapse dropped.

Those long alleles were not pure noise. aardvark compares by local haplotype, so it was crediting
small variants carried *inside* them; collapsing the record removes that correct nested content along
with the wrong top-level allele. **Stage 2 is therefore not separately shippable**, as an earlier
draft of this plan claimed. It is demolition, and Stage 3 is the rebuild: the two are halves of one
change and `--nested` must not ship with only the first.

### Stage 3 — Ploidy propagation and an explicit recursion contract

Replace emission-driven recursion with explicit descent: per child chain, ploidy is the number of
called parent alleles crossing it; ploidy 0 gives a star allele and no descent. Fix
`emit_variant`'s `return true` on a no-emission path, so "nothing was written" can no longer read as
"this site is resolved".

New flag `--nested` enables the mode. It does not replace `-A` or `--top-down`, both of which stay
until Stage 6 decides.

**Gate**: nested SNVs appear as their own records; the swallowed-SNV count falls toward zero; record
count and peak memory stay within the scheduler's `2.25 + 11.2e-6 * records` budget, refitting the
coefficient if the record count rises materially.

### Stage 4 — Marginalised symbolic likelihood — *deferred to last, not yet started*

Only if Stages 2–3 land, which they now have. Deliberately sequenced after Stages 5 and 6: it is the
one stage that changes the model rather than the output, it needs an answer for `λ` before it can be
written at all, and everything measured so far is arg-max, so it should be A/B'd against a settled
baseline rather than a moving one. Requires a defined `λ` for a symbolic class, per the decision above. A/B
against arg-max on the four tier-2 arms.

**Gate**: het-site accuracy improves without `DR` or the depth term becoming unreportable.

### Stage 5 — Linkage, phasing and mosaic for nested sites

**Partly measured, and it found an open defect.** Nested sites already flow through the linkage
layer, because linkage is on by default under `--read-likelihood`; the question was never whether
they participate but whether they should, and what breaks.

What works: ploidy propagation reaches the output. On chr20, 2,135 records are single-allele
genotypes -- nested sites called at ploidy 1 because only one parent allele crosses them -- against
114,831 diploid and 81 unphased. None carry a missing or star allele.

What breaks: **the mosaic's site-total invariant.** The mosaic must account for every emitted record,
and the harness asserts it. On the default it holds exactly (105,251 sites, 105,251 records); under
`--nested` it does not (116,789 against 117,047, a gap of 258). So 258 nested records reach the VCF
without reaching a linkage chain, and the mosaic no longer describes the whole call set.

That has to be resolved before `--nested` can ship, either by bringing those records into the chains
or by defining precisely which records the mosaic covers and changing the assertion to match. It is
not a reason to change the design, but it is a reason not to call Stage 5 done.

Still to measure: whether including nested sites *helps* phasing accuracy, which needs a whatshap
comparison against the default arm.

#### Original plan

Run both ways — nested sites in the linkage chains and out of them — and decide on measured phasing
accuracy. Whichever wins, the mosaic's site-total invariant must still hold.

### Stage 6 — Shared-flank trimming, then full evaluation

`--trim-shared-flanks`, measured *after* symbolic collapsing so it is credited only with what it
adds. Then the four tier-2 arms, the whole genome, and the PanGenie comparison refreshed.

**Gates**: SNV recall past PanGenie's 0.9659, toward the 0.9737 the swallowed count implies; SV F1
up; small-variant F1 not down; runtime and memory acceptable for the 24-contig laptop run.

## Testing

Unit tests accompany each stage as described. Beyond those:

- **A regression test that fails on today's build**: in `test/t/18_vg_call.t`, a synthetic graph with
  a SNP nested inside a snarl that also carries a deletion, asserting the SNP gets its own record and
  the parent emits no compensating substitution. This is the whole design in one case, and it should
  be written before Stage 2 so it is red first.
- The existing 272 checks in `18_vg_call.t` must pass unchanged.
- Harness assertions in the evaluation repository on the swallowed-SNV count and the same-length
  substitution false-positive count, so a recovered win cannot silently regress.

## Risks

- **Record-count growth against the memory model.** The scheduler packs contigs under
  `2.25 + 11.2e-6 * records` GB. Nested emission raises the record count; the coefficient is fitted
  and will need refitting, and the per-contig worst case (6.1 GB on chr3 today) may move.
- **REF and coordinates for off-reference nested records**, deferred by keeping `--nested-pseudo-ref`
  off by default, but it is the thing that will make off-reference calling hard later.
- **Interaction with the linkage layer**, which currently assumes one site per position with a single
  ploidy; nested sites break both assumptions and Stage 5 exists for it.
- **Trimming being credited with symbolic alleles' work**, addressed by ordering.

## Working practice: detecting when a subcommand has finished

Three times in the investigation that produced this plan, a watcher polled for a condition that could
never become true, and the work sat finished and unnoticed — once for about an hour. The failures
were: `pgrep -f "schedule_wgs.py"` from a shell whose own command line contained that string, so it
matched itself; waiting on a `BENCH_DONE` marker the program does not emit; and `pgrep -cx`, where
`-c` is not supported by BSD `pgrep`, so the loop exited immediately and reported success.

Rules for this work:

1. **Prefer running the long job in the background directly** and letting the harness notify on exit.
   A separate watcher process is the exception, for external state the harness cannot see.
2. **Wait on a PID, never on a pattern.** Capture the pid at launch and poll `kill -0 "$PID"`. A
   pattern match can match the waiter.
3. **Never wait on a completion marker without first confirming the program emits it**, by reading
   the source or a prior run's output.
4. **Check the artefact, not the process**: a finished run is one whose output file exists and whose
   exit status was captured, not one whose process has vanished — those differ when a job is killed.
5. **Verify flags on this platform.** BSD and GNU differ in `pgrep`, `head -n -N` and `sed -i`, all
   three of which have already cost time in this project.

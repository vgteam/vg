# The shape of the read-likelihood caller

A companion to [read-likelihood-genotyping.md](read-likelihood-genotyping.md), which specifies the
*model*. This page is about the *code*: what the pieces are, which of them depend on which, and
where the coupling actually is. It exists because that question was asked in review of #4990 and
deserves an answer with numbers in it rather than an opinion.

Everything below was measured on the tree as it stands, by parsing `#include` lines and class
bodies. Where a count appears, it is a count.

## The pieces

| file pair | lines | what it is |
|---|---|---|
| `read_likelihood_caller.{hpp,cpp}` | 980 | the `SnarlCaller` itself: turns a likelihood matrix into a genotype |
| `allele_likelihood.{hpp,cpp}` | 2,440 | the matrix: `P(read \| allele)` for every read and candidate |
| `site_read_source.{hpp,cpp}` | 1,259 | fetching the reads that overlap a site |
| `linkage_model.{hpp,cpp}` | 3,525 | the Li–Stephens panel model, **and** the site store that drives it |
| `read_phasing.{hpp,cpp}` | 385 | phasing adjacent sites from reads that span both |
| `regenotype.{hpp,cpp}` | 755 | re-deciding a genotype from the settled phase |
| `anchor.{hpp,cpp}` | 1,218 | the anchor output: which reads pin to which allele |
| `symbolic_allele.{hpp,cpp}` | 557 | traversal → ALT sequence, and comparison of traversals |
| `alignment_scorer.{hpp,cpp}` | 909 | the per-base scoring model shared with the mapper |

Nine file pairs, 12,028 lines. `graph_caller.{hpp,cpp}` — 10,474 lines — orchestrates them, and
`subcommand/call_main.cpp` wires them to the CLI.

## The dependency graph, measured

Each line lists what that file `#include`s from inside the family. Unit tests are excluded; they
are noted separately, because they are what makes one of the answers below come out the way it
does.

```
call_main.cpp          -> read_likelihood_caller, site_read_source, gref
graph_caller.hpp       -> linkage_model, read_phasing, regenotype, anchor, gref
graph_caller.cpp       -> read_likelihood_caller, symbolic_allele, gref

read_likelihood_caller -> allele_likelihood                        (+ snarl_caller)
allele_likelihood      -> alignment_scorer, anchor, read_phasing, site_read_source
anchor                 -> site_read_source
regenotype             -> read_phasing            (.hpp)
                          anchor                  (.cpp)
read_phasing           -> nothing
linkage_model          -> nothing
site_read_source       -> nothing in the family
symbolic_allele        -> nothing in the family
```

Three facts fall straight out of it:

1. **No subsystem names `ReadLikelihoodSnarlCaller`.** Zero mentions in `read_phasing`,
   `regenotype`, `linkage_model` or `anchor`, in either the header or the implementation.
2. **`linkage_model` and `read_phasing` include nothing from vg at all** — not one `#include "..."`
   between them. They are self-contained algorithms over plain data.
3. **The hub is `allele_likelihood.hpp`, not the caller.** It includes four of the family;
   `read_likelihood_caller.hpp` includes exactly one of them, `allele_likelihood.hpp`.

## Q1 — how wedded is a subsystem to `ReadLikelihoodSnarlCaller`?

Taking regenotyping as the example asked about: **the code is not wedded at all; the data is.**

`regenotype.hpp` includes exactly one header from the family, `read_phasing.hpp`; its `.cpp` adds
`anchor.hpp`. Its eight entry points take `PhaseSite`, `PhaseReadEvidence` and `LambdaTable` —
plain structs — and never a caller. You could hand it evidence from any source and it would work.

But there is only one source. `PhaseReadEvidence` is constructed in exactly one place,
`AlleleReadLikelihoods`' builder (`allele_likelihood.cpp:1332`), moved into the `CallInfo` by
`ReadLikelihoodSnarlCaller` (`read_likelihood_caller.cpp:124`), and consumed in `graph_caller.cpp`.
So the *contract* is narrow and replaceable while the *supply* is singular.

That is the good kind of coupling and it should be left alone. Nothing is gained by breaking the
data dependency, because a phasing subsystem has to be handed per-read per-allele evidence by
somebody; what matters is that it is handed a struct rather than a caller, and it is.

The same holds for the other three. `anchor` additionally includes `site_read_source.hpp`, for
`SiteRead` — again a data type, not a caller.

## Q2 — should `AlleleLikelihoodCalculator` live inside `ReadLikelihoodSnarlCaller`?

The instinct behind the question is right: **it is already private to the caller in everything but
syntax.** No production file outside `read_likelihood_caller.hpp` includes
`allele_likelihood.hpp`. Its only other includers are three unit-test files.

Nesting it anyway costs two things and buys none of the usual benefits:

- `AlleleLikelihoodCalculator` is a one-method **abstract base** with a single implementation,
  `GraphAlignedAlleleLikelihoodCalculator`. It is an extension point, and an extension point
  nested inside one concrete consumer is awkward to extend.
- Those three unit-test files — 180 `REQUIRE`/`CHECK` statements covering the DP walks, the
  scoring window and the band — construct the calculator directly. Nesting makes every one of them
  go through the caller, which is a worse test.

**Recommendation: no.** Record the intent instead — the header already opens with a comment saying
what the class is for, and "constructed only by `ReadLikelihoodSnarlCaller` and its tests" is a
line of documentation, not a refactor.

If the goal is specifically to shrink what `read_likelihood_caller.hpp` exposes, the effective move
is different: `allele_likelihood.hpp` is 908 lines and pulls in `anchor.hpp`, `read_phasing.hpp`,
`site_read_source.hpp` and `alignment_scorer.hpp` transitively to anyone who includes the caller.
Forward-declaring `AlleleReadLikelihoods` in the caller header would cut that, and is a
twenty-line change.

## The coupling that was not asked about, and is the real one

`VCFOutputCaller` is the base class of four callers: `VCFGenotyper`, `LegacyCaller`, `FlowCaller`
and `NestedFlowCaller`. Its declaration is 787 lines.

**39 of its 64 data members — 61% — exist only for the read-likelihood path**, along with 13 of its
43 methods: the linkage collector and its GBWT caches, `render_phases`, `phase_sites`,
`render_lambda` and its temper, the regenotyping parameters and counters, the anchor writer, and
the mosaic output state. Three of the four subclasses touch none of it.

61% is a **floor**, not an estimate. The count matches members by name, and the rule misses at
least `phase_flips`, `phase_declined`, `emit_phasing`, `current_generation` and
`MOSAIC_WALK_LIMIT`, all of which are read-likelihood-only as well — `--phased` is itself one of
the options refused without `--read-likelihood`. The denominator is also generous: three of the 64
are nested `struct` declarations rather than fields.

So the subsystems are tidy and the orchestration is not. The free functions are in files named
after what they do; it is the *composition* that has been accumulating on a shared base class
because that is where `emit_variant` lives and everything needed to reach it.

It has a concrete fix, and the counts point straight at it. Of the four subclasses, only
`FlowCaller` names any of this machinery — 3 of its 27 member definitions, and 7 declarations in
its part of the header (`apply_read_phasing`, `regenotype_resettle`, `anchor_gqn_for` and the
pending-record plumbing). `VCFGenotyper`, `LegacyCaller` and `NestedFlowCaller` name none of it, in
either file. **So the 39 members belong on `FlowCaller`, not on the base class**, and nothing but
history put them where they are.

This is worth saying plainly because it changes what a reorganisation should aim at. Moving nine
file pairs into a folder tidies the listing and leaves the 61% exactly where it is.

## Q3 — folders, enclosing classes, or both?

### Free functions: a real but small problem

Counting declarations at namespace scope across the family:

| | free functions |
|---|---|
| `regenotype.hpp` | 8 |
| `symbolic_allele.hpp` | 5 |
| `anchor.hpp` | 4 |
| `read_phasing.hpp` | 2 |
| `linkage_model.hpp` | 1 |
| `site_read_source.hpp`, `allele_likelihood.hpp`, `read_likelihood_caller.hpp` | 0 |

Twenty functions in five files. The subsystems that already have a class — `SiteReadSource`,
`AlleleLikelihoodCalculator`, `LinkageCollector`, `AnchorWriter` — have essentially none, which is
the pattern to follow rather than invent.

Of the twenty, most are genuinely free: `phase_link(a, b, cap)` and `site_read_log_odds(q0, p)` are
pure functions of their arguments with no state to hold. The ones that would benefit from an
enclosing class are the four that thread the same three or four objects through every call —
`accumulate_lambda`, `fit_calibration`, `haploid_inclusion_correction`, `phase_aware_correction`
— which between them reconstruct the same `LambdaTable`/`PhaseSite`/params triple each time.

**So: yes to a `Regenotyper` widget holding that state, no to a blanket wrapping.** Wrapping a
pure function in a class to satisfy a naming scheme makes it harder to test, not easier.

### `linkage_model` is two subsystems wearing one name

3,525 lines, split:

| | lines |
|---|---|
| `LinkageModel` — the Li–Stephens HMM, forward/backward, the transition kernel | ≈2,167 |
| `LinkageCollector` — the site store, generation resolution, allele map, retract/rescore | ≈1,275 |

`graph_caller` names `LinkageCollector` 29 times and `LinkageModel` only for `WILDCARD` (24 uses)
and `Params` (1). The HMM is already internal to the collector in everything but file layout.
**Splitting this pair is the single highest-value move in the whole reorganisation** and it is
independent of any folder decision.

### Folders

Folders are wanted, and `src/` earns them: **297 `.cpp`/`.hpp` files, 206,800 lines, flat**, with
`algorithms/`, `io/`, `config/`, `subcommand/` and `unittest/` as existing precedent. The
convention those set is `#include "algorithms/foo.hpp"` from outside and `"../bar.hpp"` from
inside.

Each folder is not free: `algorithms/` costs 23 lines of `Makefile`, `io/` and `config/` 34 each,
`subcommand/` 16 — a `_SRC_DIR`, an `_OBJ_DIR`, a `_SHARED_OBJ_DIR`, a `.d` include, a pattern rule
and a `mkdir`. That is an argument for **one flat folder, not a nested hierarchy.**

Between the two candidates in review:

- `caller/` or `lib/call/` with a `subsystems/` inside it — two folders, two sets of Makefile
  plumbing, and it asserts a hierarchy the dependency graph does not have. `linkage_model` and
  `read_phasing` do not depend on the caller; calling them its subsystems is backwards.
- `callers/read_likelihood/` — one folder, one set of plumbing, and it is honest about what the
  grouping is: everything #4990 adds.

**Recommended: one folder, `src/caller/`, flat, holding the nine pairs above.** Not
`callers/read_likelihood/`, because `graph_caller`, `snarl_caller`, `traversal_finder` and
`traversal_support` belong in the same place and are not read-likelihood-specific — and because a
`read_likelihood` folder makes the *next* caller a sibling directory rather than a sibling file,
which is more structure than two callers justify.

Deferred deliberately: whether `graph_caller` itself moves. Seven files include it, four of them
outside the calling code (`deconstructor.hpp`, `traversal_support.hpp`, `mcmc_caller.{hpp,cpp}`),
so moving it reaches past #4990's footprint for little gain; it can follow once the folder exists.

## Order of work

1. **Split `linkage_model` into `linkage_model` + `linkage_collector`.** Independent, mechanical,
   and the biggest readability win. No folder decision needed.
2. **Derive the CLI refusal set from the option table** (`call_main.cpp`) — the defect the review
   actually found: a hand-maintained list of 56 flag strings against 106 options, plus a
   hand-rolled reimplementation of `getopt_long`'s prefix resolution. **Done**; the list had
   already drifted on four live flags.
3. **Give regenotyping its enclosing class** for the four state-threading functions.
4. **Create `src/caller/` and move.** Late, because it renames the most and teaches the least.
5. **Move the 39 read-likelihood members from `VCFOutputCaller` down to `FlowCaller`**, which is
   the only subclass that uses them. Mechanical in principle — nothing else references them — but
   it touches the largest file in the calling code and is best done once the review has settled,
   so it is recorded here rather than started.

# Read-likelihood caller: a guide to the source

`vg call --read-likelihood` is built from eight modules and one driver class, `FlowCaller` in
`graph_caller.hpp`. A module is a header in `src/` and its `.cpp` file, and is named here by its
header. The source can be read in the order of the method. The method, and the terms of the
method used below (site, allele, genotype, ploidy, direct call, linkage model, strand, phase set),
are described in [read-likelihood-genotyping.md](read-likelihood-genotyping.md).

The method has four **steps**: site likelihood computation, genotyping, phasing and output.
`FlowCaller` carries them out in five **passes** over the sites.

## How a run is organised

With `--read-likelihood`, nested calling is on by default. The linkage model runs under haplotype
enumeration (the default with a GBZ graph, or `-z` or `-g`), when the panel holds at least two
haplotypes and `--linkage-weight` is not 0. When either is on, `FlowCaller` stages its records and
makes these passes. The calls that start them are near the end of `main_call` in
`subcommand/call_main.cpp`.

1. **Sweep** (`call_top_level_snarls`, then `call_snarl_internal` for each site). Genotype every
   site directly from its reads, and **stage** a record for it: keep what the record will be built
   from, but write nothing. Nested sites are genotyped here too, by recursion from their parent, at
   both ploidy 1 and ploidy 2. Each genotyped site is also filed with the linkage model
   (`record_site`), except `retain_only` chains (see below). This pass does step 1 and the direct
   call of step 2.
2. **Barrier** (`run_barrier`). Settle the genotypes with the linkage model and phase them from the
   panel, one generation at a time, parents before children. Between generations, set each child's
   ploidy from its parent's settled, phased pair, and swap in the child's answer at that ploidy. No
   reads are fetched in this pass or later ones: later passes use the per-read evidence the sweep
   kept. This pass does the rest of step 2 and the panel part of step 3.
3. **Read phasing and re-genotyping** (`phase_and_regenotype`), with `--read-phasing`. Re-decide the
   phase from reads that span several sites (`apply_read_phasing`). With `--regenotype`, then repeat
   rounds of: correct the likelihoods from the phase (`apply_regenotyping`), settle again
   (`regenotype_resettle`, which rescores the records and runs the barrier), and phase from the
   reads again. Records that get no line of their own (off-reference chains, and chains a parent's
   block record spells out) keep their sweep likelihoods. This pass does the rest of step 3.
4. **Render** (`render_retained_records`). Build each staged site's records from its settled
   genotype, add their lines to an output buffer, and collect its anchors. This pass and the next do
   step 4.
5. **Write** (`write_anchors`, then `write_variants`). Write the anchor file, then the VCF: add the
   nesting tags, sort the lines, write the mosaic file, and write each line, rewriting the quality
   fields of records whose genotype the linkage model changed.

Without nested calling and without the linkage model, nothing is staged: the sweep adds each
record's line to the buffer as soon as its site is genotyped, and the write pass follows.

The barrier keeps each site's settled pair in `linkage_phased`, which it rebuilds on every barrier
pass. It fills it only when phasing is written, and it needs it to set a child's ploidy, so where
the linkage model runs, `--no-phased` turns nested calling off. The barrier also drops a child that
no allele of its parent's settled genotype crosses, and can bring back a child that an earlier
barrier pass dropped. A child whose parent has too many candidate alleles to tell which cross it is
never revised or dropped on its parent's account.

Re-genotyping always corrects the likelihoods the sweep computed, at both ploidies, rather than
those of the previous round. It stops when no record's corrected direct call changes, when a barrier
pass changes no settled genotype, when the settled genotypes return to an earlier state, or after
`--regeno-passes` barrier passes.

## Module map

An arrow points from a module to a module that includes its header, in its own header or in its
`.cpp` file. `call_main.cpp` also builds objects of other modules, reaching their headers through
these includes. Colour marks the step each module mainly serves; the grey modules are general vg
modules that other callers use too.

```mermaid
flowchart TD
    snarls["snarls.hpp<br/>SnarlManager"]:::vg
    scorer["alignment_scorer.hpp"]:::vg
    finder["traversal_finder.hpp<br/>GBWTTraversalFinder"]:::vg
    snarlcaller["snarl_caller.hpp<br/>SnarlCaller"]:::vg
    gref["gref.hpp"]:::vg

    reads["site_read_source.hpp<br/>SiteReadSource"]:::likelihood
    phasing["read_phasing.hpp<br/>PhaseSite, read_phase_flips"]:::phase
    anchor["anchor.hpp<br/>AnchorWriter"]:::output
    likelihood["allele_likelihood.hpp<br/>AlleleReadLikelihoods"]:::likelihood
    caller["read_likelihood_caller.hpp<br/>ReadLikelihoodSnarlCaller"]:::genotype
    linkage["linkage_model.hpp<br/>LinkageModel, LinkageCollector"]:::genotype
    symbolic["symbolic_allele.hpp<br/>symbolic_allele, symbolic_diff"]:::genotype
    regenotype["regenotype.hpp<br/>accumulate_lambda, phase_aware_correction"]:::phase
    graphcaller["graph_caller.hpp<br/>FlowCaller"]:::driver
    main["subcommand/call_main.cpp"]:::driver

    reads --> likelihood
    scorer --> likelihood
    snarls --> likelihood
    phasing --> likelihood
    anchor --> likelihood
    reads --> anchor
    phasing --> anchor
    phasing --> regenotype
    snarls --> symbolic
    likelihood --> caller
    snarlcaller --> caller
    caller --> graphcaller
    linkage --> graphcaller
    symbolic --> graphcaller
    regenotype --> graphcaller
    phasing --> graphcaller
    anchor --> graphcaller
    finder --> graphcaller
    snarls --> graphcaller
    snarlcaller --> graphcaller
    gref --> graphcaller
    graphcaller --> main
    caller --> main
    reads --> main
    gref --> main

    classDef vg fill:#eeeeee,stroke:#888888,color:#000000
    classDef likelihood fill:#dbe9f6,stroke:#3a6ea5,color:#000000
    classDef genotype fill:#e3f1dc,stroke:#4a8a3a,color:#000000
    classDef phase fill:#f6e7d2,stroke:#b07a2a,color:#000000
    classDef output fill:#eadcf2,stroke:#7a4a9a,color:#000000
    classDef driver fill:#ffffff,stroke:#000000,color:#000000
```

Blue: site likelihood computation. Green: genotyping. Orange: phasing. Purple: output. White: the
driver and the command line. `gref.hpp` handles the gRef cover, a set of extra reference paths.

`linkage_model.hpp` and `read_phasing.hpp` include no other vg header: they are algorithms over
plain data, and can be read on their own. `graph_caller.cpp` reaches `ReadLikelihoodSnarlCaller`'s
results only through a `dynamic_cast`, so `FlowCaller` also runs vg's other genotypers, and skips
the linkage model and phasing when the cast fails.

## Modules by step

| Step | Module | What it holds |
|---|---|---|
| Site likelihood computation | `site_read_source.hpp` | `SiteRead`, a read as a site sees it, and `SiteReadSource`, which delivers the reads of a site from a GAM or GAF file in memory, an indexed GAM file, or a GAF-base database |
| | `allele_likelihood.hpp` | `GraphAlignedAlleleLikelihoodCalculator`, which scores each read against each candidate allele by pairing their node visits, and `AlleleReadLikelihoods`, the resulting matrix of relative likelihoods, which computes $\mathcal{L}(G)$ for each genotype |
| Genotyping | `traversal_finder.hpp` | `GBWTTraversalFinder`, which gives the candidate alleles from the panel, or `FlowTraversalFinder`, under support enumeration |
| | `read_likelihood_caller.hpp` | `ReadLikelihoodSnarlCaller`, which makes one site's direct call from its likelihoods and computes its quality fields |
| | `linkage_model.hpp` | `LinkageModel`, the hidden Markov model over the panel, and `LinkageCollector`, which holds each site's entry and runs the model over linkage chains, one generation at a time |
| | `symbolic_allele.hpp` | symbolic alleles and difference blocks, which let a nested site's variation be reported once |
| | `graph_caller.hpp` | the sweep (`call_snarl_internal`, with nested descent), the barrier (`run_barrier`) and the staged records (`PendingRecord`) |
| Phasing | `linkage_model.hpp` | the Viterbi phase (`LinkageModel::phasing`), recorded by `LinkageCollector` as each site's settled pair (`PhaseCall`) |
| | `read_phasing.hpp` | per-read evidence (`PhaseReadEvidence`, `PhaseSite`), links between sites, and the decision of each site's order (`read_phase_flips`) |
| | `regenotype.hpp` | strand log-odds, tempering, and the per-read likelihood correction |
| | `graph_caller.hpp` | `phase_and_regenotype`, which applies the two modules above and runs the barrier again |
| Output | `graph_caller.hpp` | records (`render_retained_records`, `emit_variant`, `emit_block_records`), sorting and writing (`write_variants`), and the mosaic (`write_mosaic`) |
| | `read_likelihood_caller.hpp` | the VCF fields of a record (`update_vcf_info`) and their header lines |
| | `anchor.hpp` | the anchor file (`build_site_anchors`, `AnchorWriter`) |
| Command line | `subcommand/call_main.cpp` | the options, grouped by subsystem, the `--preset` table, and the construction of every object above |

## Words the headers use

- **Record key.** A hash of the record's VCF `ID` column, which holds the printed snarl. Most tables
  that follow a site from pass to pass are keyed by it.
- **Generation.** How deep in the descent a site was genotyped: 0 for a site the sweep calls as
  top-level (including the children of a snarl it could not genotype), and its parent's generation
  plus 1 for a site reached by descent.
- **Group.** A linkage chain below the top level: the sites of one child chain, at one ploidy, under
  one parent. When the group's ploidy equals its parent's, the parent is the group's first site,
  held at its settled state. A haploid group under a diploid parent instead starts from the panel
  haplotype that the parent's carrying strand copies.
- **Settled pair.** A site's settled genotype as the linkage model phased it, a `PhaseCall` in
  `linkage_phased`. `LinkageCollector`'s entry keeps the settled genotype sorted, without phase.
- **Nested strand.** The strand of its parent on which a ploidy-1 child lies
  (`PhaseCall::nested_strand`).
- **Crossing mask.** For each child, which of the parent's candidate alleles cross it, one bit per
  allele; it is unknown when the parent has more candidates than the mask has bits.
- **Moved.** A site whose settled genotype differed from its direct call in any barrier pass. Its
  records get the moved quality fields (see the state table).
- **Staged, pending, deferred, retained, render records.** One idea, a record kept to be built
  later, in different containers: the sweep stages records into `pending_records` (nested) and
  `render_records` (top-level); the barrier moves nested ones to `deferred_pending`; the render's
  **hand-off** (`hand_off_deferred_records`) moves the surviving ones into `render_records`, and
  collects anchors for chains that get no line.
- **`retain_only`.** Under the linkage model, a child site that no allele of its parent's direct
  call crosses, and the sites nested in it. It is genotyped and staged, but not filed with the
  linkage model or written, unless the barrier finds that the parent's settled genotype crosses it.
  (Without the linkage model such a site is not genotyped at all.) A site with no reference path is
  filed anyway. The barrier does not drop a `retain_only` site on its parent's account when the
  parent has no settled pair, when the crossing mask is unknown or empty, or when the sweep has no
  answer at the settled ploidy; such a site is written. A dropped ancestor still removes it.
- **Barrier pass and round.** A barrier pass settles every generation once. A re-genotyping round
  corrects the likelihoods, runs one barrier pass, and runs read phasing again.
- **Nested.** Four senses:
  - **nested calling**, the descent into child chains in `call_snarl_internal` (`--nested`, turned
    on in the code by `VCFOutputCaller::set_symbolic_collapsing`);
  - `FlowCaller`'s constructor flag `nested`, which is `--top-down`;
  - `set_nested` and `include_nested`, which write the nesting INFO tags (`LV`, `PS`, `CH`, ...),
    under `-A`, `--top-down`, `--bottom-up` or off-reference nesting;
  - `SiteContext::nested`, which marks a site with one copy under its parent.

  `PS` names both a nesting INFO tag and `FORMAT/PS`, the phase set.

## State carried between passes

| State | Written by | Read by |
|---|---|---|
| `PendingRecord` (in `render_records`, `pending_records`, `deferred_pending`) | the sweep | the barrier, the render |
| `ReadLikelihoodCallInfo`, one per site: every genotype's likelihood, the quality fields, the per-read evidence, and the answer at the other ploidy | `ReadLikelihoodSnarlCaller::genotype`, in the sweep; the barrier exchanges it with the other ploidy's answer; re-genotyping corrects both | the barrier, pass 3, the render |
| `LinkageCollector` entries | the sweep (`record_site`); the barrier retracts an entry and records it again when it changes a child's ploidy, and records a child it brings back; re-genotyping replaces the likelihoods and direct call of each corrected record's entry (`rescore`); the render records whether each site got a line and, for a whole-site record, how each allele was written (`set_allele_map`) | the barrier; the render (`settled_traversals`); the write pass |
| `linkage_phased`, the settled pair (`PhaseCall`) of each site | the barrier, rebuilt on every barrier pass; read phasing swaps pairs in place | the barrier (parents), re-genotyping (nested strands), the render, through `render_phases`, a copy keyed by record that `build_render_phases` makes, and `write_mosaic` |
| `phase_sites`, `phase_flips` | `apply_read_phasing` | re-genotyping; `build_render_lambda`, for the anchors |
| the temper, and `render_lambda` | re-genotyping fits the temper once, on its first round; `build_render_lambda` uses it, or fits its own when re-genotyping did not, and computes the strand log-odds the anchors use | the anchors |
| `moved_quality`, the linkage model's posterior of each moved site's settled genotype, and its explained share | the barrier | `write_variants`, which rewrites `GQ` from it and recomputes `GQN` and `FILTER` from the line's `GL`; the anchors, to test whether a record moved |

## Reading order

Read the headers in this order. The glossary above covers the words that appear before the header
that explains them.

1. **`graph_caller.hpp`, the class comment of `FlowCaller` only,** and the end of `main_call` in
   `subcommand/call_main.cpp`, where the passes are called. They give the rest a frame. The rest of
   `graph_caller.hpp` is item 9.
2. **`site_read_source.hpp`.** Which reads a site gets, and the three places they come from. The
   details of the backends and their caches can be skimmed.
3. **`allele_likelihood.hpp`.** The site likelihood. Read `AlleleReadLikelihoods` first: one site's
   matrix of relative likelihoods, the mismapping probabilities, the mixture weights and the depth
   term, from which it computes $\mathcal{L}(G)$ for a genotype. Then
   `GraphAlignedAlleleLikelihoodCalculator`, which builds the matrix by pairing each read's node
   visits with each allele's, and `AlleleLikelihoodParams`, its settings. The header includes
   `anchor.hpp` and `read_phasing.hpp`, because the calculator also keeps each read's evidence for
   those two modules while the read is at hand; their types can be skipped until items 6 and 8.
4. **`snarl_caller.hpp` (`SnarlCaller` and `SupportBasedSnarlCaller` only), then
   `read_likelihood_caller.hpp`.** `ReadLikelihoodSnarlCaller` subclasses `SupportBasedSnarlCaller`,
   because `FlowCaller` is built with one. It turns one site's matrix into a direct call. Its
   `genotype` returns the genotype, and beside it a `ReadLikelihoodCallInfo` (see the state table).
5. **`linkage_model.hpp`.** `LinkageModel` first: sites, states, emissions and transitions, then
   `posteriors` and `phasing`. Then `LinkageCollector`, which files an entry for each site during
   the sweep and, in the barrier, runs the model over top-level linkage chains and groups one
   generation at a time, recording each site's settled pair. `LinkageCounters` last.
6. **`read_phasing.hpp`.** The evidence each read gives at a phaseable site, links between sites,
   and the four stages that decide each site's order. It also declares `allele_length_weights`,
   which `regenotype.cpp` and `anchor.cpp` use too.
7. **`regenotype.hpp`.** How a read's evidence at the other phaseable sites becomes per-read
   weights, built from the allele-length weights, and the correction to the likelihoods.
8. **`anchor.hpp`.** Pins, slots and anchors, and the file that records them.
9. **`symbolic_allele.hpp`, then the rest of `graph_caller.hpp`.** `symbolic_allele.hpp` compares a
   parent site's allele with the reference allele when its differences lie inside child chains, and
   cuts it into difference blocks. In `graph_caller.hpp`, read `VCFOutputCaller` (the output, and
   the linkage, phasing and anchor settings) and `FlowCaller` (the passes); skip `VCFGenotyper`,
   `LegacyCaller`, `NestedFlowCaller`, `SnarlGraph` and `GAFOutputCaller`, which serve vg's other
   genotypers. Entry points, in pass order:
   - the sweep: `call_snarl_internal`, `record_site`, `stage_render_record`, `PendingRecord`;
   - the barrier: `run_barrier`, `resolve_linkage_generation`;
   - pass 3: `phase_and_regenotype`, `apply_read_phasing`, `apply_regenotyping`,
     `regenotype_resettle`;
   - the render: `render_retained_records`, which first calls `build_render_lambda`,
     `build_render_phases` and `hand_off_deferred_records`, then, for each record,
     `settled_genotype_for`, `collect_anchors_for_record` and `emit_variant`;
   - the write pass: `write_anchors`, then `write_variants`, which calls `finalise_linkage_outputs`
     and through it `write_mosaic`.
10. **`subcommand/call_main.cpp`.** The option table, and how each option reaches the objects above.

The `.cpp` files follow the same order, and `traversal_finder.hpp` can be read whenever the
candidate alleles need explaining. Each module has unit tests in `src/unittest/<module>.cpp`, and
`allele_likelihood.hpp` has a second file of them, `allele_likelihood_scoring.cpp`.
`test/t/18_vg_call.t` tests the whole command.

## Where the code departs from this order

- `allele_likelihood.hpp` fills two types that belong to later modules, the anchor evidence and
  the phase evidence, because only it has each read's alignment in hand. Under `--anchors-out` it
  keeps only the anchor evidence, and the phase evidence is derived from it by `phase_evidence_of`,
  a static function in `graph_caller.cpp` that no header declares.
- `allele_length_weights` is declared in `read_phasing.hpp` but is also used by `regenotype.cpp`
  and `anchor.cpp`.
- `linkage_model.hpp` holds two concepts, the model (`LinkageModel`) and the code that drives it
  over the snarl tree (`LinkageCollector`).
- `graph_caller.hpp` holds the passes and the output. Most of the read-likelihood state lives on
  `VCFOutputCaller`, a base class that vg's other callers share, although only `FlowCaller` uses it.
- `call_main.cpp` constructs `GraphAlignedAlleleLikelihoodCalculator` and `LinkageCollector` itself.
- `--anchors-out` changes more than the output. It turns on the genotyping of off-reference chains,
  which adds groups and sites to read phasing, so the written phase can change, and under
  `--regenotype` genotypes too.

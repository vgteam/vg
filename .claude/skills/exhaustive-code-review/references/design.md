# Design lens

Ask only: is the change built from parts a reader can name, each in the place
it belongs, using what vg already has, with interfaces whose contracts can be
stated? Read the `class-identity` skill in `.claude/skills/` first; it is the
standard for what a class must be.

Code generators are prone to every problem below. They work a few files at a
time, so a change can be reasonable where it sits and wrong for the codebase.
They add code rather than restructure or reuse it. They satisfy tools and tests
literally: a linter's complaint becomes a pile of helpers, a test becomes a
reason for a flag. They err in both directions on structure, writing one
function that does four jobs in one place and a layer of one-implementation
abstractions in the next, so look for both. Review as if the author could not
see the whole project.

## First, write down the design

Before looking for problems, write a short inventory, for the report and for
yourself:

- each new class, struct, free function group, file and concept the change
  adds, with one sentence on what it is (not what it does);
- where each lives, and what layer that is (graph library, snarl
  decomposition, an algorithm under `src/algorithms/`, a caller, a subcommand);
- how data and control flow between them for one ordinary input.

Count nested classes, inner structs and new blocks inside existing functions as
parts too.

Then test what you wrote, part by part. Each of these in your own sentence is a
finding about the part, not a flaw in your sentence:

- it lists what the part holds or does instead of saying what it is;
- it defines the part by negation or by contrast ("nodes that bound no
  snarl", "the other case", "everything not handled above");
- it names the part by the code that uses it ("the stuff the X needs");
- it needed a word that vg does not use for that idea.

If a part has no one-sentence identity, or the flow cannot be traced, that is
the first finding. Changes that add many parts and no overview to read them in
("which files do I read, in what order?") get that as a finding too.

## Parts

- **A class that is not a thing.** A cluster of related members with no
  identity: a "pair" of other things, a bag of settable fields with several ways
  to build it, a "bench" or "helper" or "manager" that only groups functions.
  Ask what the class *is*; if the answer is a list of what it holds, say so.
- **A concept defined by what it is not.** "Bounding-no-snarl", "non-trivial
  other case". A real concept has a positive definition.
- **An invented concept where vg has an idiom.** "Along/against the chain"
  where vg says `is_reverse` relative to the chain; a new word for a glossary
  term in CLAUDE.md. Different names for the same phase of an algorithm in the
  code and its docs.
- **Machinery for exactly two.** A test or adapter built around comparing two
  implementations, where the natural design handles any number (0, 1, or n).
  Two implementations of the same interface are not two roles: a "pair" of them
  is a bundle, even when one is treated as the reference. Each implementation
  should be checked against stated expectations, or all of them against each
  other.
- **Names that hide the concept.** `data`, `result`, `process()`, `info`,
  `handle_case()`; void no-argument methods whose names don't say what state
  they change (`resolve_linkage()` beside `finalise_linkage_outputs()`).

## Placement and layering

- **A specific concept in a generic layer.** Algorithm-specific state or
  counters added to a base class or a helper class meant for many users; a
  helper class that quietly stops being one.
- **A lower layer reaching back up, or the reverse.** Code above an
  abstraction (such as `SnarlDecomposition`) that brings back the types it hides
  (`Snarl`, `SnarlManager`).
- **A free function that knows too many classes' concepts.** It belongs to one
  of them.
- **General code buried in a feature.** A graph algorithm belongs in
  `src/algorithms/` (or libhandlegraph's algorithms); a string or numeric
  helper in `src/utility.hpp`.
- **Files in the wrong place.** Generated code, write-ups and scripts that no
  person has cleared belong in `artifacts/` (see CLAUDE.md); stray data or
  scripts committed under `test/`.

## Reuse and fit

- **Reinventing what exists.** For every new function, class or block that
  does a recognizable job, look for existing code that does it, as described in
  "Finding existing code" below. Shapes seen:
  a number-to-snarl map where `SnarlManager::translate_snarl_num()` exists;
  hand-built paths-to-GBWT code; a second connected-components routine.
- **A facade that doesn't abstract.** A new layer that reshapes an adequate
  API, such as a vector-returning wrapper over a callback (`for_each_*`) API.
  vg prefers the callback style; use it directly.
- **Departing from the conventions of sibling implementations.** A new
  implementation of an interface that represents things differently from the
  existing ones (cyclic chains, trivial snarls, orientation flags) without a
  reason that would convince the owner of the other implementation. Where the
  existing implementation is itself inconsistent, the fix is to state the rule
  in the interface's documentation and make both follow it; say so.
- **Keeping a legacy case alive.** Code paths for inputs the project is
  removing (unary snarls, `CactusSnarlFinder`) instead of a stated precondition
  that excludes them.
- **A parallel way of doing the same job.** A new class, module, error-handling
  style or logging style that does what an existing one does, in a different
  shape. Each change reads well alone; together they leave the codebase with
  two of everything.
- **Hand-written loops where vg or the standard library has the operation.**
  A raw loop that `<algorithm>`, a handle graph algorithm or a project helper
  already provides.
- **Semantic duplication.** The same logic written again in a different shape,
  so duplicate-code checkers miss it: a second parser for the same format, a
  second scoring rule for the same quantity, a copy of an existing loop with
  one condition changed.

## Finding existing code

Grep finds code with the same names; reimplementations usually have different
ones. For each new part that does a recognizable job:

1. Write the job as one sentence without the new code's names ("build a GBWT
   from a graph's paths"), and list a few verbs and nouns for it, including the
   CLAUDE.md glossary terms.
2. Run a similarity search with Semble, which finds functions and classes that
   do a similar job whatever they are called. Use the first of these that
   works:
   - Semble's MCP tools, if this session has them: search for "semble" with
     ToolSearch.
   - The `semble` command, if `command -v semble` finds it.
   - Otherwise install the `semble` package from PyPI. Use `pipx install
     semble` if pipx is present. If not, make a virtual environment outside the
     checkout and install it there:

     ```
     venv="${XDG_DATA_HOME:-$HOME/.local/share}/semble-venv"
     python3 -m venv "$venv" && "$venv/bin/pip" install semble
     ```

     Then run it as `"$venv/bin/semble"`. If installing fails, say so in the
     report and go on without it.

   Run `semble search "<job sentence>" src -k 20 --format json` for each job;
   this finds the most. Run `semble find-related <file> <line> src -k 20
   --format json` on a line inside each new class and each block that does a
   recognizable job, not on every function; it returns many near-copies of
   test setup code and few real matches. Give `<file>` relative to `src`, as
   results print it (`unittest/foo.cpp`, not `src/unittest/foo.cpp`). Add the
   public header directories of the dependencies the change uses as further
   paths (for example `deps/libhandlegraph/src/include`); if a submodule
   directory is empty in this checkout, Semble returns nothing for it without
   saying so, so check that it has files. The first run builds an index, in
   about 10 seconds for `src/`; later queries take under a second.

   The change's own files crowd the top of the results. Drop them using the
   `file_path` field of each result in the JSON, for example with `jq` or
   `python3`. The rest are candidates, many of them wrong. A hit that is not
   itself the existing code often calls it, so open each likely one, including
   tests, and follow what it uses.

   Semble finds code with similar purpose and different names. It misses
   functions whose names and bodies share little wording with the new code,
   which the next two steps find; do them too.
3. Search by signature: grep headers for functions that take the same input
   types and return the same type as the new code.
4. Read the public methods of the classes whose objects the new code takes
   apart; the job is often already a method there.
5. Look in the places shared code goes: `src/algorithms/`, `src/utility.hpp`,
   the `*_helper.hpp` headers, libhandlegraph's `algorithms/`, libbdsg, and the
   subcommands that already do the same job.

## Size and shape

- **Piling on with no plan.** New machinery bolted onto an existing class or a
  file that is already thousands of lines, with no boundary drawn around it, or
  a feature added to a function that needed restructuring first. Ask how the
  complexity is meant to be managed. A large change with no answer is
  `blocking`: it needs redesign, not fixes.
- **A function or class with several jobs.** One function that parses, decides,
  does I/O and formats; a class whose methods use disjoint sets of its fields.
  Split along the jobs.
- **Fragmentation.** Many small helpers that exist to satisfy a linter's length
  or complexity limit and are not reusable jobs in their own right. One
  coherent 100-line function beats ten fragments a reader must reassemble.
- **New API with no caller.** Groundwork for a later change is fine, but the
  intended caller and flow must be written down (in the header's file comment or
  the commit message), and the pieces must fit each other: what one returns,
  the next should take.
- **Speculative generality.** An interface or base class with one
  implementation and no second in sight; a factory, provider or wrapper that only
  delegates; options, template parameters, hooks or virtual methods nothing uses
  yet. Ask whether deleting the layer would change any behavior.
- **Dead parts.** New functions, overloads, members or headers that nothing
  calls, and old names left behind after a rename.
- **Logic that can only be tested through the whole pipeline.** A decision
  (which anchors to keep, which allele wins) buried in a long function, so tests
  need flags or elaborate setup to reach it. Split out the decision as a
  function with inputs and an answer.

## Interfaces and contracts

- **One function, two jobs.** If the doc comment cannot describe it as one
  coherent operation, it is two functions.
- **Requirement or detail?** A doc comment that mixes what callers may rely on
  with how it happens to work, so neither can be told apart.
- **Unstated or unenforced preconditions.** Say what inputs are allowed; then
  either check them (and throw) or name them as preconditions.
- **Relaxed invariants.** An implementation that loosens an interface's
  guarantee (ordering, completeness, orientation) instead of meeting it.
- **Many hooks.** A function taking several callbacks to customize its steps
  is a class with virtual methods to override. Header or output text built by
  callers and passed in as strings is an interface (`get_format_headers()`).
- **State that is filled in lazily or held for later.** Members marked
  `mutable` and filled on first use are a synchronization problem under OpenMP;
  rank them `important`, and build the table up front instead. An object that
  stores a graph reference only to hand it back, or that both stores a graph
  and takes one per call and mixes them up, should take the graph per call as
  vg's methods do. An adapter may keep a reference to the object it adapts.

## Errors and edge cases as design

- **Silent fallbacks.** A clamp, default object or empty result where bad
  input should be reported; a trivial substitute algorithm run when the real one
  fails. Garbage in one field suggests garbage in others; fail or warn instead
  of repairing. Comments saying "fallback", "for now" or "simple approach" are a
  lead.
- **Swallowed errors.** A `catch` that logs and carries on as if it worked, or a
  broad `catch (...)` meant for one case that hides every other.
- **Stubs presented as finished.** A function with a realistic signature whose
  body only logs, returns a constant or ignores its parameters.
- **Defensive code for impossible inputs.** Guards, `max()` clamps and
  try/catch around things that cannot happen, which make a reader wonder when
  they can. State the invariant once instead.
- **Checks in the wrong place.** Validation done by the Makefile, the CLI or a
  caller that belongs in the component, or the reverse.

## Options and configuration

- **Hacks on shared structures.** Extending a standard option table or struct
  with tags for one subsystem, instead of giving the subsystem its own.
- **Inputs an installed binary cannot see.** Data files read from the source
  tree at run time; build them in or install them.
- **Names that promise something else.** An option or setting whose name says
  one thing (`--realign`) while the code does another.

## Code shaped by its tests

- Production code that checks for a particular test input, contains a test's
  literal expected values, or has a flag or path that only tests use.
- Tests loosened in the same change as the code they test.

## The author's process showing through

- Comments, names or docs that refer to how the code was developed: earlier
  versions, the runs used to tune it, the conversation, "now", "new approach".
- Documentation that denies interpretations no reader would form.
- A model or algorithm described one way in the docs and implemented another,
  with the gap explained away instead of closed.

## Scope

- Leave test structure to the tests lens, which judges whether shared test code
  still reads clearly.
- A bug you find while tracing the design (wrong coordinates, a crash) still
  goes in the report; mark it "correctness, found by the design lens".
- Report a pre-existing design problem only when the change extends or builds
  on it. Say that it is pre-existing.

## Reporting

Name the part and the principle in one line, then the concrete change: which
class to split or delete, where code should move, which existing function to
call. Rank a design finding:

- `blocking` when the change needs redesign rather than fixes: many new parts
  with no identity, or hundreds of lines added to an already large class or
  file with no boundary drawn;
- `important` when it will be copied or built on: a class or concept in a
  header, a base class, a new subsystem's shape, a nested class with no
  identity, a duplicated algorithm, lazily filled shared state;
- `minor` when it stays inside one function.

Leave out matters of taste that none of the patterns above covers.

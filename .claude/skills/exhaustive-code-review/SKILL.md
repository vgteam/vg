---
name: exhaustive-code-review
description: >-
  Run a thorough, multi-lens code review of a diff, branch, PR, or file. Use
  this whenever the user asks to review code, "look over" or "check" a change,
  audit a patch before merge, do a pre-commit or PR review, or asks whether a
  change is safe to land — even if they don't say the word "review." Especially
  use it for C/C++ and other systems code, file parsers, and any change where
  correctness, memory safety, portability, concurrency, or security matter. Do
  not hand-wave a quick skim; this skill exists to make the review deep,
  consistent, and structured every time. In vg, run it before code is pushed,
  before a pull request is opened, and before handing code to the user for
  their own review.
---

# Exhaustive Code Review

## The idea

A single reviewer asked to watch for everything at once tends to catch the
obvious and miss the subtle, because attention diffuses. This skill instead runs
several **focused lenses**, each with one mental model, then **triages** their
combined output into a single ranked report. Two things make that work instead
of drowning the author in noise:

1. **Route before you review.** Not every lens applies to every change. Pick the
   lenses the diff actually warrants (see the router below) so you go deep where
   it matters and skip lenses that would only generate filler.
2. **Let deterministic tools do the mechanical work.** A sanitizer or static
   analyzer is more reliable than eyeballing for the things it can detect. Run
   the tool, then have the review *interpret and prioritize* its output — don't
   re-derive by hand what ASan/TSan/clang-tidy already know, and don't claim a
   class of bug is absent if the tool that would find it wasn't run.

The output is judged on signal: a short, ranked list of real issues beats an
exhaustive list of nitpicks the author will learn to ignore.

## Workflow

Follow these steps in order. Don't skip step 1 or step 5.

### 1. Gather context

Establish what you're reviewing and why before judging anything:

- Get the diff. Use the range the caller names, if any; otherwise follow the
  `vg-code-review` skill. For a branch/PR: `git diff <base>...<head>`
  (three-dot, so you see only what this branch added). For uncommitted work: `git diff` and
  `git diff --staged`. For a file, read the whole file, not just changed lines —
  context determines whether a change is correct.
- Read the intent: PR description, linked issue, commit messages. A change can be
  flawless and still wrong if it doesn't do what it set out to do.
- Note the blast radius: is this a leaf function, a widely-called library
  routine, a subcommand parsing user files, an on-disk format? This drives both
  routing and severity.
- **Load the `vg-code-review` skill** for what to review and skip in vg, the
  patterns that have turned up real bugs, and how to settle whether a bug is
  new. With no target named, review the current branch against its merge base
  with `origin/master` (or `master` when there is no `origin`), plus any
  uncommitted changes.

### 2. Run deterministic tools first

Before reading for issues a tool could find, run what's available and feed the
results into the relevant lenses. See `references/tooling.md` for specifics. In
brief, for C/C++:

- Build with `-Wall -Wextra` (and ideally `-Werror` on the diff) and read every
  new warning — many "review findings" are just unread compiler warnings.
- `clang-tidy` / `cppcheck` for static analysis, where installed, and the
  `lint-new-code` skill (SonarQube) on new classes and functions.
- ASan + UBSan for memory and undefined-behavior bugs; TSan if the code is
  threaded; Valgrind where sanitizers aren't an option. A sanitizer build of vg
  is a full rebuild: ask before starting one.
- Run the existing tests; note coverage gaps over the changed lines.

If a tool can't be run in this environment, say so explicitly in the report
rather than silently implying that class of bug was checked. When there is no
vg build, follow "Without a build" in the `vg-code-review` skill instead of the
build steps above.

### 3. Route: choose the lenses

Use this table. Run a lens if any trigger matches; skip it otherwise and say you
skipped it (and why) so the author knows it wasn't an oversight.

| Lens | Run it when the diff… | Reference |
|------|-----------------------|-----------|
| Correctness | …contains any logic, control flow, arithmetic, parsing, or state change (almost always) | `references/correctness.md` |
| Memory & safety | …is C/C++/unsafe, handles buffers/pointers/allocation, or parses untrusted/external input | `references/memory-and-safety.md` |
| Concurrency | …touches threads, shared mutable state, signals, locks, atomics, or async callbacks | `references/concurrency.md` |
| Portability | …uses sizes/endianness/alignment assumptions, platform APIs, compiler extensions, or builds for multiple arches | `references/portability.md` |
| Performance | …adds work to a hot path, allocates in loops, changes algorithmic complexity, or touches I/O patterns | `references/performance.md` |
| Design | …adds or changes a class, struct, file, function signature, option, or concept (almost always for new features) | `references/design.md` |
| Maintainability | …is non-trivial in size or will be read/extended by others (almost always, but lowest severity) | `references/maintainability.md` |
| Tests | …adds or changes any test | `references/tests.md` |
| Comments | …adds or changes any comment or doc comment, or adds a concept a reader would need explained | `references/comments.md` |

A pure mechanical refactor with green tests might warrant only Correctness +
Maintainability. A new subcommand option that reads a file warrants Correctness +
Memory & safety at minimum. Reach for the others by the triggers above.

### 4. Run each selected lens

For each lens, **read its reference file**, then pass over the diff with only
that lens's questions in mind. Keep passes separate — don't blend them — so each
gets undiluted attention. If you have subagents available, spawn one per lens in
parallel, each loading its own reference file and every skill that file names;
if not, do them sequentially. Brief each lens subagent with the diff, the
intent, its reference file and the constraints, not with what you suspect, so
its read is not anchored on yours. Give Design, Tests and Comments their own subagents
even when you do the other lenses yourself: a reader hunting for bugs stops seeing
whether tests and comments can be understood. Collect raw findings per lens;
don't filter yet.

### 5. Synthesize and triage

This step is what keeps the review usable. Take all raw findings and:

- **Deduplicate.** The same underlying bug often shows up in several lenses
  (a buffer overrun is both a memory and a security finding). Merge into one.
- **Score** each finding by **severity** and **confidence** (see below).
- **Rank**, and **cut the low-confidence chatter.** A `suggestion`-level item you
  are 40% sure about is noise; drop it or fold it into a single "minor notes"
  line. Spend the author's attention on what's real.
- **Weigh readability findings fairly.** A test that cannot be checked by
  reading it, or a new concept that no comment explains, keeps the next person
  from telling whether the code is right. Rank such findings `important`, not
  `minor`, and do not cut them as chatter.
- **Verify before asserting.** If a finding depends on behavior elsewhere in the
  codebase, check it rather than assuming. A confident wrong finding costs more
  trust than a missed nit.

### 6. Emit the report

Use the template in `references/report-template.md`. Lead with the verdict and
the blocking issues; put nits and praise at the end. Write the report to be read
once by a busy person: apply the `humanizer` skill's pattern removal to it, but
not its "personality" section.

### 7. Fix, then record the review

Do this step only when the review covered the code in the current checkout
that is about to be pushed or handed over, not when reviewing someone else's
pull request. Fix what the user agrees should be fixed, then run the tests that
cover the fixes. Then record that this state of the code has been reviewed:

```
python3 .claude/hooks/review_state.py mark
```

The hooks in `.claude/settings.json` use that record to stop urging a review
before a push or a pull request. Any later edit makes the state unreviewed
again; a commit of the reviewed changes does not. Do not mark a state whose
blocking findings are still open unless the user says to send it anyway.

## Severity and confidence

**Severity:**

- `blocking` — must fix before merge: wrong results, memory corruption, security
  hole, data loss, crash on reachable input.
- `important` — should fix; may block depending on context: missing error
  handling, a real but narrow bug, a significant perf regression on a warm path.
- `minor` — style, naming, small clarity issues, optional improvements.
- `praise` — call out genuinely good work; it's signal too, and it builds trust
  in the rest of the review.

**Confidence** — state it when below certain. "Possible use-after-free if `b` can
be freed by the callback (didn't trace all callers)" is far more useful than
either silent omission or false certainty. Prefer phrasing findings as the
specific condition under which the bug triggers.

## Output discipline

- Every finding gets a `file:line`, a one-line statement of the problem, *why*
  it's a problem, and a concrete suggested fix or the question to resolve.
- Quote the smallest snippet needed; point to the line, don't reproduce blocks.
- Be collaborative: a reviewer's job is to surface risk, not to gatekeep. Frame
  judgment calls as questions where the author may know something you don't.
- If the change is clean, say so plainly and stop. Length is not thoroughness.

## Reference files

- `references/tooling.md` — deterministic tools to run and how to read them
- `references/correctness.md` — logic, edge cases, error handling
- `references/memory-and-safety.md` — memory lifetime, UB, untrusted input, injection
- `references/concurrency.md` — races, deadlock, signal/async safety
- `references/portability.md` — sizes, endianness, platform/compiler assumptions
- `references/performance.md` — complexity, allocation, I/O, hot paths
- `references/design.md` — parts, layering, reuse, contracts; patterns code generators produce
- `references/maintainability.md` — structure, naming, testability
- `references/tests.md` — whether tests are worth their cost and readable
- `references/comments.md` — whether comments are true, needed and readable
- `references/report-template.md` — the output format

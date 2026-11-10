---
name: vg-code-review
description: Scope rules and bug patterns for reviewing changes to vg (C++ variation graph toolkit). Use alongside the exhaustive-code-review skill on any vg diff, branch, or pull request, to know what to review, what to skip, which patterns have turned up real bugs, and how to settle whether a bug is new.
---

# Reviewing vg changes

## Scope

- Review C++ under `src/` (including `src/subcommand/` and `src/unittest/`), the
  Bash-TAP tests under `test/t/`, scripts under `scripts/`, and the build
  (`Makefile`, `source_me.sh`, CI files).
- Leave `deps/` out of the review, except for a change to which submodule
  commit vg uses. Reading a dependency at its pinned commit, to check what vg
  relies on, is fine.
  Then check that the new commit exists upstream and that vg still builds and
  passes the tests that exercise that dependency.
- Leave out generated or bulky test data (`test/graphs/`, `test/small/` data
  files) unless the change is the data. Check that a new test file is small and
  that a test reads only files that are committed or that the test makes.
- Get the diff against where the branch started: `git diff $(git merge-base
  origin/master HEAD)` for committed work, plus `git diff HEAD` for uncommitted
  work. Read per commit (`git show <hash>`) when a commit history is messy;
  `git diff A..B -- file` mixes in other people's merged changes.

## Patterns worth looking for

- **A doc comment that states a precondition or a rule nothing enforces.** This
  is the most productive single thing to grep a diff for. Shapes seen: "callers
  must pass a non-empty path name" with a caller that passes an option's empty
  default; "only called on a snarl with no unary children" with no check and a
  caller that builds its `SnarlManager` from `CactusSnarlFinder`; a "handles are
  always forward" claim over code that flips them. Fix one of two ways: make the
  code enforce it (an assertion, a check that throws, a filter) or rewrite the
  comment to say what the code does. Leaving the claim standing costs most,
  because the next caller believes it.
- **A comment narrower than the code it guards.** If the comment says "skip
  reference paths" and the test only checks that the path name has a `#`, the
  test also catches haplotype paths. Reconcile them.
- **A check that compares a value with its own parts.** `kept + dropped !=
  total`, where `total` was built by adding those two, can never fire. Compare
  what came out with what went in: the number of input records, minus each skip
  counter, against what was written. Every `continue` in the loop then has to
  bump exactly one counter.
- **A changed branch condition that re-routes cases nobody meant to touch.** A
  condition usually serves more input cases than the one in the change. List the
  combinations it routes (forward or reverse handle × start or end of a snarl ×
  trivial or not) and check that each still goes the same way, except the one
  meant to change. Prefer a guard built from a variable that already says the
  real intent over a coarser stand-in.
- **Orientation mistakes.** `get_sequence` on a reverse handle is already
  reverse-complemented; doing it again is a bug. Offsets into a node are in the
  handle's orientation. `follow_edges(h, true, ...)` goes left of `h` *as
  oriented*. A traversal of a snarl read from the end needs both its order and
  each visit flipped. A child that is reversed in its parent chain flips again.
- **Snarl decomposition assumptions.** New code may require a `SnarlManager`
  with no unary snarls; check it says so and that callers use
  `IntegratedSnarlFinder`. Check for looping chains (first and last node the
  same), trivial chains, tips inside snarls, and snarls that are not start-end
  connected.
- **Distance index sentinel values.** Minimum distances are the largest
  `size_t` when there is no walk; adding to one overflows to a small number.
  Look for arithmetic on a distance that is not first checked against
  `std::numeric_limits<size_t>::max()`.
- **Path name parsing that ignores PanSN parts.** Code that splits on `#`, or
  compares a whole path name, must handle the phase block and the subrange
  (`[start-end]`) the path metadata API gives. Prefer `PathMetadata` calls over
  string handling.
- **Shared state under OpenMP.** `for_each_handle(..., true)` and `#pragma omp
  parallel for` loops that push to a shared vector, bump a shared counter, or
  write to a stream need a critical section, an atomic, or per-thread buffers
  (`omp_get_thread_num()`). Reading a `HandleGraph` in parallel is fine;
  mutating one is not.
- **An output that claims success it did not check.** A stream write that is
  never checked, a `--progress` message printed before the work is done, a
  subcommand that exits 0 after an error message. vg subcommands report a usage
  or input error through the `Logger` in `src/log.hpp`: `logger.error() << ...`
  prints `error[vg <subcommand>] ...` and exits non-zero once the statement
  ends. Code that goes on after it is dead; code that writes to `cerr` itself
  and keeps going may be a bug.
- **Tests that pass without testing the change.** A Bash-TAP test that only
  checks the exit code, a `REQUIRE` inside a loop that never runs, an expected
  value copied from the program's own output with the bug in it. Back the fix
  out and check the test goes red.
- **Tests that cost more than they are worth.** Reviewers push back on these as
  hard as on bugs. Tests longer or harder to follow than the code they test.
  Setup steps written as tests (an `is` on the exit code of `vg convert` or of
  writing an input file), which test some other subcommand or nothing. A test
  that an option got the short letter it was given. A long test aimed at one
  imagined mistake, where a few small tests of the real cases (each orientation
  of the read against the path, say) would cover it. A test whose logic, with
  many moving parts, is not explained. For each test, ask what bug it would
  catch that no other test would; suggest cutting the ones with no answer.
- **Test inputs built at run time.** GFA, GAM or JSON written out from strings
  inside a `.t` file or a unit test is hard to read and review. Prefer small
  input files committed under `test/`, or a heredoc with real newlines and tabs.
  Scratch files go in the current directory with fixed names, as in the rest of
  the `.t` files, not in a temporary directory of their own.
- **A subcommand option added without its help text, or a help text with no
  option.** `scripts/lint.py` checks that the help text, `long_options`, the
  short-option string and the `switch` agree; run it rather than checking by
  hand, and do not ask for a test of it.

## Settling "regression or pre-existing"

For each bug that matters, build the merge base in a scratch worktree (`git
worktree add --detach <dir> $(git merge-base origin/master HEAD)`), run the same
input against both builds, and say which it is. A vg build takes a long time;
ask before starting one, and reuse an existing build of the base if there is one.

## Without a build

A fresh worktree has no vg build, and a build takes a long time. Unless the
user asks for one, review by reading: run `scripts/lint.py` and anything else
that needs no build, compile a small standalone program to settle a question
about the language, and trace each finding through the code. Mark such findings
"traced, not run", and say in the report what could not run.

## Running sanitizers

`make clean && make jemalloc=off asan=on` builds vg with AddressSanitizer. That
rebuilds everything, so do it in a scratch worktree, not in the checkout the
user is working in, and ask first. Run the unit tests for the changed code
(`bin/vg test "[tag]"`) and the Bash-TAP tests that use the changed subcommand.

# Deterministic tools for a vg review

Run what applies, read every line of new output, and say in the report which
tools were not run. Check what is installed (`command -v clang-tidy cppcheck
valgrind`) instead of assuming.

## Compiler warnings on the changed files

vg builds with few warnings turned on. Rebuild just the changed objects with
more. `CXXFLAGS` from the environment is added to the Makefile's own flags:

```
rm -f obj/<file>.o
CXXFLAGS="-Wall -Wextra -Wshadow" make obj/<file>.o 2>&1 | grep -A3 '<file>'
```

Subcommand sources build to `obj/subcommand/<name>_main.o`, unit tests to
`obj/unittest/<name>.o`. Report only warnings in lines the change touched or
in code it made reachable; the rest are pre-existing.

## Static analysis

- `clang-tidy <file> -- -std=c++17 -Isrc -Iinclude` with the
  `bugprone-*,performance-*,clang-analyzer-*` checks, if installed.
- `cppcheck --enable=warning,performance --std=c++17 -Isrc <file>`, if
  installed. It does not follow vg's includes well; expect noise.
- The `lint-new-code` skill for new or redesigned classes and functions.

## Tests

- Unit tests: `bin/vg test "[tag]"` for the tags of the changed code; find them
  with `grep -n TEST_CASE src/unittest/<area>*.cpp`.
- Command tests: from `test/`, `prove -v t/<NN>_<subcommand>.t` for each
  subcommand the change can affect.
- Note changed lines that no test reaches. To show that a test covers a fix,
  back the fix out, rebuild, and watch the test fail.

## Memory and undefined behavior

- AddressSanitizer: `make clean && make jemalloc=off asan=on` in a scratch
  worktree. It is a full rebuild; ask first. Run the tests above under it.
- UBSan: add `-fsanitize=undefined` through `CXXFLAGS` in the same build.
- Valgrind: `valgrind --error-exitcode=1 bin/vg <subcommand> ...` on a small
  input from `test/`. Slow but needs no rebuild. jemalloc can hide leaks from
  it; trust its invalid read and write reports more than its leak report.
- ThreadSanitizer only when the change adds sharing between threads, and only
  with `OMP_NUM_THREADS` small: OpenMP runtimes give false positives.

## Reading tool output

A tool finding is a lead, not a verdict. Confirm each one against the code and
fold it into the lens it belongs to.

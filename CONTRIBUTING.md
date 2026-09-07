# Contributing guidelines for `vg`

Thank you for helping us make `vg` better! Below are some ways to help:

## I want to report a bug

1. Search our open [issues](https://github.com/vgteam/vg/issues) to make sure it hasn't been reported before.
    - If you find the same bug, please don't just post "I have this problem too".
    That does not help us with fixing the bug. We need details.
    Adding more/updated debugging details is helpful, however.
    - Many bugs have similar results, such as causing `vg` to crash.
    If the command being run is significantly different, then you have a new bug.
2. Collect debugging details. If we can reproduce the problem ourselves,
    then it is much easier for us to help you.
    - What **version** are you running? Find this out with `vg version`.
    - How much **memory** does your system have? Both disk space and RAM.
    This is often relevant during indexing operations, especially for GCSA.
    - What **commands** did you run? This includes both the command that led to the error,
    and the pipeline which led up to that. Since many errors are caused by bad input
    we can check if your previous commands did something unexpected.
    - What **data** can we use to reproduce the error? If you can't share data
    full data for legal or file size reasons, that's fine. We prefer a minimal example anyways.
    Ideally you can copy/paste or upload directly to GitHub.
        - For graphs, can you remove chromosomes or samples to create a small example?
        - For reads, can you find a specific read/pair, e.g. by binary search?
        - For VCFs, GFFs, etc. can you narrow down which records are problematic?
    - What did you **expect** to happen? Sometimes this is obvious (e.g. "not crash")
    but sometimes it isn't (e.g. "find a different alignment"). Please be specific.
    - If you have a **stacktrace/logfile**, please provide the entire log.
3. [Open a bug](https://github.com/vgteam/vg/issues/new?template=report-a-bug.md).

## I want to contribute code/functionality

You'll need a full `vg` installation in order to build locally. Instructions are in the
[README](https://github.com/vgteam/vg/blob/master/README.md#installation).

**You are responsible for the code you commit.** If you use a generative AI tool such as Claude Code,
you must inspect the output and make sure that it actually does what you want it to do. 
Ideally you would at least write your own unit tests.

### File organization

- Dependencies go in `deps/`. If you want to add a new one, you'll need to:
    - Add it to [`.gitmodules`](https://github.com/vgteam/vg/blob/master/.gitmodules).
    - Edit the [`Makefile`](https://github.com/vgteam/vg/blob/master/Makefile) with necessary instructions.
- Helper scripts (not part of the main build) go in `scripts/`
- Command-line tests go in `test/`
    - Most subdirectories here are for test inputs/outputs
    - There is a small [`Makefile`](https://github.com/vgteam/vg/blob/master/test/Makefile)
    for the `build_graph.cpp` test.
    - Testing scripts go in the `test/t` directory. See "Command-line unit tests".
- Documentation (including a version of the wiki) is in `doc/`.
- Source code for the main build goes in `src/`
    - The `src/algorithms/` directory is for drop-in modules/functions
    that were factored out into their own files.
    - The `src/io/` directory is like `algorithms` but specifically for I/O code.
    Each different kind of filetype generally gets its own I/O file.
    - The `src/subcommand/` directory is for files defining subcommands. See "Subcommands".
    - The `src/unittest/` directory is for direct C++ unit tests. See "C++ unit tests".
- If you're touching other directories, I really hope you know what you're doing.

New `.cpp` files will be auto-discovered by the build system.

### Command-line unit tests

We have a suite of Bash-TAP tests within `test/t/`. Run with `prove -v` for proper error handling.

```bash
cd test
prove -v t/02_vg_construct.t
```

The top of the test file has boilerplate:

```bash
#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="C" # force a consistent sort order 

plan tests 33
```

The last line is very important. It must be the actual number of unit tests run.
Otherwise, even if all tests pass, the file overall will fail for incorrect test count.

Each unit test takes the form of

```bash
is $query $correct "name"
```

You can run any commands you like before/around/within these `is` checks.
The first argument will be compared to the second. If they differ, the test fails.
The name should be descriptive, e.g. `"The 1mb graph has the expected number of nodes"`.

Ideas for unit tests:
- Make sure that bad input leads to an error. Use `$?` to check exit status of the last command
- Make sure that output is well-formed
- Make sure that two outputs that should be equivalent are the same (`diff`)
- Make sure that logfiles have necessary elements (e.g. error messages)

### C++ unit tests

You can call individual functions for direct C++ testing. These live in `src/unittest/`.
Run with

```bash
vg test
```

Or, to run tests for a particular tag:

```bash
vg test [chain_items]
```

Our unit tests use the Catch v2 framework. Put `#include "catch.hpp"` at the top of your `.cpp` test file.
The tests themselves should be inside of the `vg::unittest` namespace.
Tests are split up into test cases, each of which can have any number of sections.
Use `REQUIRE()` for asserts. For example, here's a unit test from `src/unittest/aligner.cpp`:

```c++
TEST_CASE("Full-length bonus can hold down the left end", "[aligner][alignment][mapping]") {
    VG graph;
    
    TestAligner aligner_source_1;
    aligner_source_1.set_alignment_scores(1, 4, 6, 1, 0);
    const Aligner& aligner_1 = *aligner_source_1.get_regular_aligner();
    
    TestAligner aligner_source_2;
    aligner_source_2.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner_2 = *aligner_source_2.get_regular_aligner();
    
    Node* n0 = graph.create_node("AGTGCTGAAGT");
    
    string read = string("AATGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph, true);
    aligner_2.align(aln2, graph, true);
    
    SECTION("left end is detached without bonus") {
        REQUIRE(aln1.path().mapping_size() == 1);
        REQUIRE(aln1.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln1.path().mapping(0).position().offset() == 2);
        REQUIRE(aln1.path().mapping(0).edit_size() == 2);
        REQUIRE(aln1.path().mapping(0).edit(0).from_length() == 0);
        REQUIRE(aln1.path().mapping(0).edit(0).sequence() == "AA");
    }
    
    SECTION("left end is attached with bonus") {
        REQUIRE(aln2.path().mapping_size() == 1);
        REQUIRE(aln2.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln2.path().mapping(0).position().offset() == 0);
        REQUIRE(aln2.path().mapping(0).edit_size() == 3);
        REQUIRE(aln2.path().mapping(0).edit(0).from_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(0).to_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(0).sequence() == "");
    }
}
```

This test is named `"Full-length bonus can hold down the left end"`
and has three tags, `[aligner]`, `[alignment]`, and `[mapping]`.
It also has two sections, `"left end is detached without bonus"` and `"left end is attached with bonus"`.
The code is run from top to bottom. To run the entire test case:

```bash
vg test "Full-length bonus can hold down the left end"
```

To run only the first section (plus the setup code before the sections):

```bash
vg test "Full-length bonus can hold down the left end" --section "left end is detached without bonus"
```

Unit test files can also have helper functions which are called by the unit tests.
For an example of this see `src/unittest/chain_items.cpp`.

### Subcommands

Functionality within `vg` is split up over many subcommands, e.g. `vg convert` converts between file types.
Subcommands are the command-line user interface. Thus, we have strict rules about how they are set up.
These rules are enforced by a unit test defined by `scripts/lint.py`.
You can run this script yourself to check in advance of opening a PR.

Each subcommand has a file in `src/subcommand/`. Make sure to `#include "subcommand.hpp"`.

#### Register a subcommand

- All subcommands must be in the autocomplete files `autocomplete.sh` and `autocomp.bash`.
Each has a line `opts="cmd1 cmd2"` etc. When adding or removing a subcommand, update the list.
- Non-deprecated subcommands must be included in the manpage generators:
    - The lists `doc/vgmanmd.desc.md` must have at least one instance of the subcommand.
    Make sure to format the links properly or they won't be accepted.
    - The `cmds` list at the top of `doc/vgmanmd.py` must have the subcommand.

#### Helptext

Each subcommand has helptext which can be printed by running the command with no options,
or by using the `--help`/`-h`/`-?` command-line options.
(Having all of those options work is required by `scripts/lint.py`.)

```
usage: vg gamcompare aln.gam truth.gam >output.gam

options:
  -d, --distance-index FILE  use distances from this distance index
                             instead of path position annotations
  -r, --range N              distance within which to consider reads correct
  -n, --rename Q=T           treat query contig Q as truth contig T (may repeat)
  -I, --ignore T             ignore the given truth contig name (may repeat)
  -o, --output-gam FILE      output GAM to FILE instead of standard output
  -T, --tsv                  output TSV (correct, mq, aligner, read)
                             compatible with plot-qq.R to standard output
  -a, --aligner STR          aligner name for TSV output ["vg"]
  -s, --score-alignment      get alignment correctness score (higher is better)
  -t, --threads N            number of threads to use
  -h, --help                 print this help message to stderr and exit
```

There are several requirements for helptext:

- It must be printed by a `help_<name>()` function.
- There must be a `usage:` line.
- Options must be prefaced by exactly two spaces.
- Shortform options (e.g. `-r`) must be followed by a comma, a space, and then the longform.
- Longform options (e.g. `--range`) must all be vertically aligned, even if no shortform is present.
- Argument hints (e.g. `N`) must be separated from the longform by a single space, and must be all-caps.
- Options in the helptext must be composed of alphanumeric characters and hyphens.
- Option descriptions must be vertically aligned.
- Lines with option descriptions must be no more than 80 characters long.
Descriptions may be multiline. This check is a little naive and can be tricked,
including by adding newlines in the code which aren't reflected in the output.

#### Command-line options

In order to correctly set up a command-line option it must be put in four separate places.
(If you think that's annoying, I agree, but at least `scripts/lint.py` will catch mistakes.)

- **Helptext**: not all options have to appear here. See "Helptext" for requirements.
- **`long_options[]` array**: all options must appear here, in the form `{"longform", arg_type, 0, shortform}`.
    - `"longform"` must be a string.
    - `arg_type` must be either `no_argument` or `required_argument`.
    - `shortform` must be a single character (NOT a bare non-quoted integer) or an `ALL_CAPS` variable name.
    Note that `ALL_CAPS` variables must be before the `long_options[]` array as `constexpr int`s.
    Shortforms may repeat, and the first longform is retained. This is to allow longform aliases.
- **`getopt_long` string**: all options with a shortform must appear here. Order does not matter.
Options which take an argument must have a `:` after their shortform.
The `vg gamcompare` one is `"h?d:r:I:n:o:Ta:st:"`.
- **`switch(c)` block**: all options must appear here. The default case must exit/abort. Fallthroughs are handled.
    - Each case must be of the form `case <shortform>:` or `default:`.
    - The shortform must be a single character or an `ALL_CAPS` variable name.
    - If the case takes an argument, it must use `optarg`.
    - If it *intentionally* doesn't work (e.g. is deprecated, calls `exit`) it may or may not use `optarg`.
    - If the case takes a filename argument, it must use either `require_exists()` (for an input)
    or `ensure_writable()` (for an output).
    - If the option is `--threads`, it must use `set_thread_count()`.

Information must match between all of these places. If `long_options[]` says an option takes an argument,
then it must have an argument in the helptext (if it appears there),
a `:` in the `getopt_long` string, and `optarg` usage or an error in the `switch(c)` block.
If an option has a shortform somewhere then it must have a shortform everywhere it appears, etc.

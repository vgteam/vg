---
name: lint-new-code
description: Lint the classes and functions you added or redesigned with SonarQube for IDE, run headless by extension-problems-cli, and fix what it reports in them, keeping any classes and functions it makes you split into things a reader can picture. Use before finishing a task that adds or redesigns a class or adds a function.
---

# Lint the code you designed

A static analyzer finds unused variables, needless copies, missing `const` and functions too
tangled to follow, which are cheap to fix while the code is new. Run it on what you wrote, and
fix what it finds there. Problems in code you did not touch belong to another change.

## Get the tools

1. **Check for Node.js 20.12 or newer and Java 21 or newer** (`node --version`,
   `java -version`). The extension's analyzer runs on Java. If either is missing, ask the person
   you are working for to install it.
2. **Keep one copy of the extension, in the tool's cache directory,** and download it only if
   it is not there. The file is 170 MB, and a copy in a repository, a temporary directory or a
   scratch directory gets downloaded again by every session:
   ```bash
   EXT_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/extension-problems-cli"
   EXT="$EXT_DIR/sonarlint-vscode-6.1.0.vsix"
   [ -f "$EXT" ] || { mkdir -p "$EXT_DIR" && curl -L -o "$EXT" \
       https://github.com/SonarSource/sonarlint-vscode/releases/download/6.1.0%2B91156/sonarlint-vscode-6.1.0.vsix; }
   ```
   The tool keeps its unpacked copy of the extension, and the C++ analyzer the extension
   downloads, in its own cache directory, so they are shared by every session too.
3. **For C++, bring the compilation database up to date** with `make compile_commands.json`
   from the repository root. It needs the dependencies built, and it picks up new source
   files. The analyzer reads each file's compiler flags from it, and reports nothing useful for a
   file it has no flags for.

## Run it

4. **Analyze the files that hold your new code,** headers and sources, in one run. New code
   lives in `artifacts/`, or `artifacts/unittest/` for tests, until a person clears it for
   `src/`:
   ```bash
   npx -y github:adamnovak/extension-problems-cli --extension "$EXT" \
       --compile-commands compile_commands.json --format pretty \
       artifacts/new_thing.hpp artifacts/new_thing.cpp
   ```
   Leave out `--compile-commands` for files other than C and C++, such as Python scripts. A run
   on a few files takes about a minute. The first run on a machine also downloads the C++
   analyzer and takes longer, so run it in the background. Exit status 2 means some files got no
   report; run those again.
5. **Pick out the findings that are yours.** The tool reports on whole files, and an old file
   can have hundreds of findings. Each finding names one line, and "[+N locations]" at the end
   of its message means it involves other lines that the report does not list, so a finding
   that your change caused can name a line you did not touch.
   - For a new file, or one you mostly rewrote, every finding is yours.
   - For a file you added to, read the findings on and around your code, and the findings about
     the size of a class or function you grew, which name its first line.
   - Before a large change to existing files, lint them first and save the report. Afterwards,
     compare the two reports for findings that are new.

## Fix what it finds

6. **Fix small problems directly:** unused variables, copies that can be `const` references,
   missing `const` on methods, `contains()` instead of `count()`, and the like.
7. **Split a class or function that is too big along real parts,** as the `class-identity`
   skill says. Find a part of the machine the big object can hand work to, and give each new
   class its organizing metaphor. Give each function you split out one
   operation that one sentence can explain, and name it for that operation. Dealing methods into
   two classes, or cutting a function into `step1()` and `step2()`, satisfies the analyzer but
   makes code nobody can picture, which is worse than the warning.
8. **Never suppress a finding.** When no real part exists and the code cannot
   be made simpler, or a rule conflicts with how the surrounding code is
   written, leave the finding and tell the person you are working for which
   findings remain, and why you did not fix each one.
9. **Run it again** after fixing, until each finding in your code is fixed or on that list. Then
   rewrite any comments you wrote or changed while fixing, as the `intelligible-comments` skill
   says.

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
2. **Keep the extension in the tool's cache directory,** so every session reuses one download:
   ```bash
   EXT_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/extension-problems-cli"
   EXT="$EXT_DIR/sonarlint-vscode-6.1.0.vsix"
   [ -f "$EXT" ] || { mkdir -p "$EXT_DIR" && curl -L -o "$EXT" \
       https://github.com/SonarSource/sonarlint-vscode/releases/download/6.1.0%2B91156/sonarlint-vscode-6.1.0.vsix; }
   ```
3. **For C++, bring the compilation database up to date** with `make compile_commands.json`
   from the repository root. It needs the dependencies built, and it picks up new source
   files. The analyzer reads each file's compiler flags from it, and reports nothing useful for a
   file it has no flags for.

## Run it

4. **Analyze the files that hold your new code,** headers and sources, in one run:
   ```bash
   npx -y github:adamnovak/extension-problems-cli --extension "$EXT" \
       --compile-commands compile_commands.json src/new_thing.hpp src/new_thing.cpp \
       > lint.json
   ```
   Leave out `--compile-commands` for files other than C and C++, such as Python scripts. A run
   on a few files takes about a minute. The first run on a machine also downloads the C++
   analyzer and takes longer, so run it in the background. Exit status 2 means some files got no
   report; run those again.
5. **Keep the problems in your code** by piping the report through the script next to this
   file: `python3 .claude/skills/lint-new-code/new_problems.py < lint.json`. It prints the
   problems on lines your change added since the merge base with `origin/master` (give another
   base commit as its argument). Size and shape problems, such as too many methods in a class or
   too complex a function, are reported at the first line of the class or function, so it lists
   those for every changed file in a separate section. Act on those only for classes and
   functions you added or grew.

## Fix what it finds

6. **Fix small problems directly:** unused variables, copies that can be `const` references,
   missing `const` on methods, `contains()` instead of `count()`, and the like.
7. **Split a class or function that is too big along real parts,** as the `class-identity`
   skill says. Find a part of the machine the big object can hand work to, and give each new
   class its "A `Name` is a ___ that ___" sentence. Give each function you split out one
   operation that one sentence can explain, and name it for that operation. Dealing methods into
   two classes, or cutting a function into `step1()` and `step2()`, satisfies the analyzer and
   makes code nobody can picture, which is worse than the warning.
8. **Suppress a warning only with a reason.** If no real part exists and the code cannot be made
   simpler, or a rule conflicts with how the surrounding code is written, put `// NOSONAR` on
   the reported line with a comment saying why, and tell the person you are working for. Without
   the reason, the next reader cannot tell a judgment from a shortcut.
9. **Run it again** after fixing, until your new code has no problems left but suppressed ones.
   Then rewrite any comments you wrote or changed while fixing, as the `intelligible-comments`
   skill says.

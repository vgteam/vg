#!/usr/bin/env python3
"""
Print the problems that extension-problems-cli found on lines a change added.

Usage, from the repository root:
    extension-problems-cli ... | python3 .claude/skills/lint-new-code/new_problems.py [BASE]

Reads the tool's JSON output on standard input. BASE is the commit the change starts from, by
default the merge base with origin/master; lines added since it, committed or not, count as
new. Prints one line per problem as FILE:LINE:COLUMN: SEVERITY: MESSAGE [CODE], then the
files the tool did not analyze, if any.

A problem with the size or shape of a whole class or function is reported at its first line,
which a change that only grew it did not add. These STRUCTURE_RULES problems are printed for
every changed file, in a section of their own.
"""

import json
import os
import re
import subprocess
import sys


# Too many methods in a class, too complex a function, too deeply nested, too long a function
# or lambda, too many parameters.
STRUCTURE_RULES = {"cpp:S1448", "cpp:S3776", "cpp:S134", "cpp:S138", "cpp:S1188", "cpp:S107"}


def git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def added_lines(base):
    """Map each changed file, relative to the repository root, to the line numbers it added."""
    added, path = {}, None
    for line in git("diff", "-U0", "--no-color", base, "--").split("\n"):
        if line.startswith("+++ "):
            path = line[6:] if line.startswith("+++ b/") else None
        elif line.startswith("@@") and path:
            start, count = re.match(r"@@ -\S+ \+(\d+)(?:,(\d+))? @@", line).groups()
            count = 1 if count is None else int(count)
            added.setdefault(path, set()).update(range(int(start), int(start) + count))
    # A new file that is not yet added to git has every line new.
    for path in git("ls-files", "--others", "--exclude-standard").split("\n"):
        if path and path.endswith((".cpp", ".hpp", ".h", ".cc")):
            with open(path, errors="replace") as f:
                added[path] = set(range(1, sum(1 for _ in f) + 2))
    return added


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else \
        git("merge-base", "HEAD", "origin/master").strip()
    root = git("rev-parse", "--show-toplevel").strip()
    report = json.load(sys.stdin)
    added = added_lines(base)
    structure = []
    for p in report.get("problems", []):
        path = os.path.relpath(p["file"], root)
        line = (f"{path}:{p['line']}:{p['column']}: {p['severity']}: {p['message']} "
                f"[{p.get('code')}]")
        if p["line"] in added.get(path, ()):
            print(line)
        elif path in added and p.get("code") in STRUCTURE_RULES:
            structure.append(line)
    if structure:
        print("\nSize and shape problems elsewhere in changed files:")
        print("\n".join(structure))
    for path in report.get("unanalyzedFiles", []):
        print(f"{os.path.relpath(path, root)}: not analyzed")


if __name__ == "__main__":
    main()

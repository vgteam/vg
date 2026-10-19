#!/usr/bin/env python3
"""
Print the problems in an extension-problems-cli report that lie on lines a change added.

Usage, from the repository root:
    python3 .claude/skills/lint-new-code/new_problems.py [BASE] < REPORT.json

BASE is the commit the change starts from, by default the merge base with origin/master. A line
added since BASE is new, whether or not it is committed, and so is every line of a source file
that git does not track yet.

Print each new problem as FILE:LINE:COLUMN: SEVERITY: MESSAGE [CODE], then each file the tool
did not analyze.

The analyzer reports a problem with a whole class or function, such as too many methods, at the
first line of the class or function. A change that only grew the class did not add that line,
so print the problems of STRUCTURE_RULES for every changed file, in a section of their own.
"""

import json
import os
import re
import subprocess
import sys


# The rules that judge a whole class or function: too many methods, too high a cognitive
# complexity, too deep nesting, too long a function or lambda, and too many parameters.
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
    # Every line of a source file that git does not track yet is new.
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

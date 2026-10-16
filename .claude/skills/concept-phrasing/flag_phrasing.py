#!/usr/bin/env python3
"""
Print comment lines changed since a base commit that may introduce an idea by negation, by a
headless noun phrase, by an unbounded quantifier, by ellipsis, by a vague relation, or by
borrowing another component's behaviour as a justification, as described in SKILL.md next to this script.

Usage, from the repository root:
    python3 .claude/skills/concept-phrasing/flag_phrasing.py [BASE]
    python3 .claude/skills/concept-phrasing/flag_phrasing.py --files FILE...

BASE defaults to origin/master, and only the changed lines of C++ files under src/ are checked.
With --files, every comment in the named files is checked. Each changed comment run
(consecutive comment lines) is joined into one string, so a phrase broken across lines is found.
Each match prints as `path:line: shape: text`. A match is a candidate for a person or agent to
judge, so the exit status is always 0.
"""

import re
import subprocess
import sys

SHAPES = [
    ("negation", r"\b(that|which|who)\s+(do|does|did|can|could|will|would|is|are)\s*n[o']t\b"),
    ("negation", r"\b(that|which|who)\s+\w+\s+no\s+\w+"),
    ("negation", r"\bnon-\w+"),
    ("negation", r"\bor nothing if\b"),
    ("headless noun phrase", r"\bwhat\s+(a|an|the|we|it|this|they|you|each|its)\b"),
    ("headless noun phrase", r"\b(where|whichever|whatever)\s+(a|an|the|we|it|this)\s+\w+\s+(lies|is|goes|reads)\b"),
    ("definition through a user", r"\b(will|would)\s+be\s+(built|made|used|written)\s+from\b"),
    ("definition through a user", r"\bthe\s+(stuff|things|data)\s+(for|we need|the caller needs)\b"),
    ("unbounded quantifier", r"\b(every|all|any)\s+(comparison|check|law|test|difference|problem)s?\b"),
    ("unbounded quantifier", r"\bthe\s+known\s+\w+"),
    ("contract word", r"\bonce\b"),
    ("ellipsis", r"\b(gives?|has|have|lists?|with)\s+(\w+\s+)?none\b"),
    ("ellipsis", r"\b\w+'s\s*[.,;:)]"),
    ("ellipsis", r"\b(the|a)\s+(first|last|other|next|previous)\s+one's\b"),
    ("vague relation", r"\b(meets?|meeting|touch(es)?|sees?)\s+(it|them|the|a|an|each|its|their)\b"),
    ("borrowed justification", r"\b(in the (same )?order|the way|where|as)\s+(\w+\s+){0,4}?\w+\(\)\s+(finds?|produces?|returns?|does|do|starts?|lists?)\b"),
    ("borrowed justification", r"\b(to match|matching|matches|which is (also )?(what|where|how))\b"),
]

COMMENT = re.compile(r"^\s*(?://+!?<?|/\*+|\*+/?)\s?(.*)$")
TRAILING = re.compile(r"\S.*?(?://+!?<?)\s?(.*)$")
STRING = re.compile(r'"(?:\\.|[^"\\])*"')


def changed_lines(base):
    """Map each changed C++ file under src/ to the set of its changed line numbers."""
    diff = subprocess.run(["git", "diff", "-U0", base, "--", "src"], capture_output=True, text=True,
                          check=True).stdout
    changed, path = {}, None
    for line in diff.split("\n"):
        if line.startswith("+++ "):
            path = line[6:] if line.startswith("+++ b/") else None
            if path is not None and not path.endswith((".cpp", ".hpp", ".h")):
                path = None
        elif path is not None and line.startswith("@@"):
            match = re.match(r"^@@ -\S+ \+(\d+)(?:,(\d+))? @@", line)
            start, count = int(match[1]), int(match[2] or 1)
            changed.setdefault(path, set()).update(range(start, start + count))
    return changed


def comment_runs(text):
    """Yield (line numbers, joined text, offsets) for each run of consecutive comment lines."""
    lines, parts, offsets, length = [], [], [], 0
    in_block = False
    for number, line in enumerate(text.split("\n"), 1):
        stripped = line.strip()
        if in_block:
            # Inside a /* */ block, every line is comment text, with or without a leading *.
            body = re.sub(r"^\*+/?\s?", "", stripped).replace("*/", "").strip()
            in_block = "*/" not in stripped
        else:
            match = COMMENT.match(line)
            if not match:
                match = TRAILING.match(STRING.sub('""', line)) if "//" in line else None
            body = match.group(1).strip() if match else ""
            in_block = stripped.startswith("/*") and "*/" not in stripped
            body = body.replace("*/", "").strip()
        if body:
            lines.append(number)
            offsets.append(length)
            parts.append(body)
            length += len(body) + 1
        elif parts:
            yield lines, " ".join(parts), offsets
            lines, parts, offsets, length = [], [], [], 0
    if parts:
        yield lines, " ".join(parts), offsets


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--files":
        to_check = {path: None for path in sys.argv[2:]}
    else:
        to_check = changed_lines(sys.argv[1] if len(sys.argv) > 1 else "origin/master")
    for path, changed in sorted(to_check.items()):
        with open(path, encoding="utf-8", errors="replace") as source:
            text = source.read()
        for lines, run, offsets in comment_runs(text):
            if changed is not None and not changed & set(lines):
                continue
            for shape, pattern in SHAPES:
                for match in re.finditer(pattern, run, re.IGNORECASE):
                    line = max(n for n, o in zip(lines, offsets) if o <= match.start())
                    if changed is None or line in changed:
                        context = run[max(0, match.start() - 40):match.end() + 40]
                        print(f"{path}:{line}: {shape}: ...{context}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Check prose in doc/*.md and in source comments for phrasings the project does not accept.

Run with scripts/lint_prose.py from the repository root; `make lint` runs it.

Each rule is a regular expression matched against one paragraph of a Markdown file, or one run of
consecutive comment lines in a .cpp or .hpp file, joined with spaces so that a phrase broken
across lines is still found. A match prints `path:line: message` and fails the check.
"""

import os
import re
import sys

RULES = [
    # A thing named only by what will later be made from it ("keep what the record will be built
    # from"). Name the data, or its type, instead.
    (r"\bwhat\b[^.;:]{0,60}?\b(will|would)\s+be\s+(built|rendered|made|written|computed)\s+from\b",
     "names data by what will be built from it; name the data or its type"),
    (r"\bwhat\b[^.;:]{0,60}?\b(is|are)\s+(built|rendered|made)\s+from\b",
     "names data by what is built from it; name the data or its type"),
    (r"\b(kept|held|stored)\s+to\s+be\s+(built|rendered|written|made)\b",
     "defines a thing by its later use; say what it holds"),
    (r"\bto\s+be\s+(built|rendered|written|made)\s+later\b",
     "defines a thing by its later use; say what it holds"),
    # The same through what some code would do with the data ("what the caller would write").
    (r"\bwhat\b[^.;:]{0,40}\bwould\s+(write|have\s+returned)\b",
     "names data by what would be done with it; name the data or its type"),
    (r"\bkept\s+for\s+later\b",
     "defines a thing by its later use; say what it holds"),
]

# A line that is only a comment: //, ///, //!, /*, /**, or a * continuation line.
COMMENT = re.compile(r'^\s*(?://+!?<?|/\*+|\*+/?)\s?(.*)$')


def units(path, text):
    """Yield (line numbers, joined text, offsets) for each paragraph or comment run."""
    buf, lines, offsets, length = [], [], [], 0
    for number, line in enumerate(text.split('\n'), 1):
        if path.endswith('.md'):
            body = line.strip() or None
        else:
            match = COMMENT.match(line)
            body = match.group(1).strip() if match else None
        if body:
            offsets.append(length)
            lines.append(number)
            buf.append(body)
            length += len(body) + 1
        elif buf:
            yield lines, ' '.join(buf), offsets
            buf, lines, offsets, length = [], [], [], 0
    if buf:
        yield lines, ' '.join(buf), offsets


def main():
    paths = []
    for top, exts in (('doc', ('.md',)), ('src', ('.cpp', '.hpp'))):
        for root, _, files in os.walk(top):
            paths += [os.path.join(root, f) for f in files if f.endswith(exts)]
    ok = True
    for path in sorted(paths):
        with open(path, encoding='utf-8', errors='replace') as f:
            text = f.read()
        for lines, para, offsets in units(path, text):
            for pattern, message in RULES:
                for match in re.finditer(pattern, para, re.IGNORECASE):
                    line = max(n for n, o in zip(lines, offsets) if o <= match.start())
                    print(f"{path}:{line}: {message}")
                    ok = False
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

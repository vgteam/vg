#!/usr/bin/env python3
"""Check prose in doc/*.md and in source comments for phrasings the project does not accept,
and that library comments do not name command-line options.

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

# A command-line option, such as --read-likelihood. Library code describes a setting by what it
# means; only the subcommands in src/subcommand/ name options. Comments in src/subcommand/ and
# src/unittest/ may name them.
OPTION = re.compile(r'(?<![\w-])--[a-z][a-z0-9]*(?:-[a-z0-9]+)*')
OPTION_MESSAGE = "names a command-line option in a library comment; describe the setting instead"
# Library comments that named an option before this check existed, by their text after the
# comment marker.
OPTION_COMMENTS_ALLOWED = {
    '--progress output rather than emitted as a zero-length path.',
    '--progress output.',
    'Load a translation file (created with vg gbwt --translation) and return a backwards mapping',
    'Load a translation file (created with vg gbwt --translation) and return a mapping',
    'NestedFlowCaller, and --bottom-up is rejected with -L.  It is kept because the consequence of',
    'is this variant" behind --cluster-min-len in BOTH vg call and vg deconstruct.  It is',
    'post-genotyping ALT merging (vg call -L / --cluster-min-len).  Deliberately NOT named',
    'sampled graph. The same floor as `vg paths --min-gref-len`.',
    'similarity metric and the same core-length gate as "vg deconstruct -L/--cluster-min-len" (a',
    'uncomment to make vg map --debug very interesting',
}
# A string literal, so that a // inside one does not start a comment.
STRING = re.compile(r'"(?:\\.|[^"\\])*"')
MARKER = re.compile(r'^\s*(?://+!?<?|/\*+|\*+/?)\s?')


def comment_lines(text):
    """Yield (line number, comment text) for every comment in C++ source, including one that
    follows code on its line, with the comment marker removed."""
    in_block = False
    for number, line in enumerate(text.split('\n'), 1):
        if in_block:
            if '*/' in line:
                in_block = False
            yield number, MARKER.sub('', line).strip()
            continue
        # Blank out string literals, keeping their length, to find where a comment starts.
        code = STRING.sub(lambda m: '"' + 'x' * (len(m.group()) - 2) + '"', line)
        line_comment, block_comment = code.find('//'), code.find('/*')
        if line_comment >= 0 and (block_comment < 0 or line_comment < block_comment):
            yield number, MARKER.sub('', line[line_comment:]).strip()
        elif block_comment >= 0:
            if '*/' not in code[block_comment:]:
                in_block = True
            yield number, MARKER.sub('', line[block_comment:]).strip()


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
        library = path.endswith(('.cpp', '.hpp')) and not path.startswith(
            (os.path.join('src', 'subcommand') + os.sep, os.path.join('src', 'unittest') + os.sep))
        if library:
            for line, comment in comment_lines(text):
                if OPTION.search(comment) and comment not in OPTION_COMMENTS_ALLOWED:
                    print(f"{path}:{line}: {OPTION_MESSAGE}: {OPTION.search(comment).group()}")
                    ok = False
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

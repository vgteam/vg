#!/usr/bin/env python3
"""Check prose in doc/*.md, in source comments and in unit test names for phrasings the project
does not accept, and check that source comments and unit test names outside src/subcommand/ name
no long command-line option.

Run with scripts/lint_prose.py from the repository root; `make lint` runs it.

Each rule is a regular expression matched against one paragraph of a Markdown file, the comments
on one run of consecutive lines of a .cpp or .hpp file, or the name of one test case, scenario or
section in src/unittest/. Lines are joined with spaces, so a phrase broken across lines is still
found. A match prints `path:line: message` and fails the check.

Unit tests exercise library code, so their comments and names follow the library's rule: they
describe a setting, and do not name the option that sets it.

The option check finds long options such as --min-mapq only. In a comment, a short option such
as -L looks like a minus sign, another program's flag (grep -v) or part of a diagram, so the
writer and the reviewer check for short options themselves.
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
    # A condition called a gate ("the discount is gated on the length change"). The word names no
    # condition, so say what must be true for the behaviour to happen ("the discount applies when
    # the length change is at least 50 bp").
    (r"\bgat(?:e|ed|es|ing)\b",
     "says something is gated; say what must be true for it to happen"),
]

# A command-line option, such as --read-likelihood. Only the subcommands in src/subcommand/ have
# options, so only their comments name them. Library code and its unit tests describe a setting
# by what it means.
OPTION = re.compile(r'(?<![\w-])--[a-z][a-z0-9]*(?:-[a-z0-9]+)*')
OPTION_MESSAGE = "names a command-line option outside src/subcommand/; describe the setting instead"
TEST_OPTION_MESSAGE = ("names a command-line option in a unit test name; "
                       "name the library setting or describe the behaviour")
SUBCOMMAND = os.path.join('src', 'subcommand') + os.sep
UNITTEST = os.path.join('src', 'unittest') + os.sep
# Comments outside src/subcommand/ that may name an option, each by its text after the comment
# marker.
OPTION_COMMENTS_ALLOWED = {
    'Respects catch.hpp --rng-seed option to `vg test`.',
    '--progress output rather than emitted as a zero-length path.',
    '--progress output.',
    'Load a translation file (created with vg gbwt --translation) and return a backwards mapping',
    'Load a translation file (created with vg gbwt --translation) and return a mapping',
    'NestedFlowCaller, and --bottom-up is rejected with -L.  It is kept because the consequence of',
    'is this variant" behind --cluster-min-len in BOTH vg call and vg deconstruct.  It is',
    'post-genotyping ALT merging (vg call -L / --cluster-min-len).  Deliberately NOT named',
    'sampled graph. The same floor as `vg paths --min-gref-len`.',
    'uncomment to make vg map --debug very interesting',
}
# A string literal, so that a // inside one does not start a comment.
STRING = re.compile(r'"(?:\\.|[^"\\])*"')
# The marker that opens or continues a comment: //, ///, //!, //<, /*, /**, or the * at the start
# of a block comment's continuation line.
MARKER = re.compile(r'^\s*(?://+!?<?|/\*+|\*+/?)\s?')
# The opening of a Catch TEST_CASE, SCENARIO or SECTION, up to the end of its name: one string
# literal, or several that the compiler concatenates.
TEST_NAME = re.compile(r'\b(?:TEST_CASE|SCENARIO|SECTION)\s*\(\s*((?:"(?:\\.|[^"\\])*"\s*)+)')
LITERAL_BODY = re.compile(r'"((?:\\.|[^"\\])*)"')


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
    """Yield (line numbers, joined text, offsets) for each paragraph of a Markdown file, or for the
    comments on each run of consecutive lines of C++ source. An empty line, or an empty comment,
    ends a run."""
    if path.endswith('.md'):
        bodies = ((number, line.strip()) for number, line in enumerate(text.split('\n'), 1))
    else:
        bodies = comment_lines(text)
    buf, lines, offsets, length = [], [], [], 0
    for number, body in bodies:
        if buf and (not body or number != lines[-1] + 1):
            yield lines, ' '.join(buf), offsets
            buf, lines, offsets, length = [], [], [], 0
        if body:
            offsets.append(length)
            lines.append(number)
            buf.append(body)
            length += len(body) + 1
    if buf:
        yield lines, ' '.join(buf), offsets


def test_names(text):
    """Yield (line number, name) for each test case, scenario and section in a unit test file."""
    for match in TEST_NAME.finditer(text):
        name = ''.join(LITERAL_BODY.findall(match.group(1)))
        yield text.count('\n', 0, match.start(1)) + 1, ' '.join(name.split())


def problems(path, text):
    """Yield (line number, message) for each problem in one file."""
    for lines, para, offsets in units(path, text):
        for pattern, message in RULES:
            for match in re.finditer(pattern, para, re.IGNORECASE):
                yield max(n for n, o in zip(lines, offsets) if o <= match.start()), message
    if not path.endswith(('.cpp', '.hpp')) or path.startswith(SUBCOMMAND):
        return
    for line, comment in comment_lines(text):
        option = OPTION.search(comment)
        if option and comment not in OPTION_COMMENTS_ALLOWED:
            yield line, f"{OPTION_MESSAGE}: {option.group()}"
    if not path.startswith(UNITTEST):
        return
    for line, name in test_names(text):
        for pattern, message in RULES:
            if re.search(pattern, name, re.IGNORECASE):
                yield line, message
        option = OPTION.search(name)
        if option:
            yield line, f"{TEST_OPTION_MESSAGE}: {option.group()}"


def main():
    paths = []
    for top, exts in (('doc', ('.md',)), ('src', ('.cpp', '.hpp'))):
        for root, _, files in os.walk(top):
            paths += [os.path.join(root, f) for f in files if f.endswith(exts)]
    ok = True
    for path in sorted(paths):
        with open(path, encoding='utf-8', errors='replace') as f:
            text = f.read()
        for line, message in sorted(problems(path, text), key=lambda problem: problem[0]):
            print(f"{path}:{line}: {message}")
            ok = False
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

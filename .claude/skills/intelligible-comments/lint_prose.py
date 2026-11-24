#!/usr/bin/env python3
"""Check the prose you wrote in Markdown files, source comments and unit test names for phrasings
the project does not accept, and check that source comments and unit test names outside
src/subcommand/ name no long command-line option.

Usage, from anywhere in the repository:

    lint_prose.py [--base REF] [PATH ...]

With paths, check those files in full. Without, check the lines added since REF (HEAD by default,
so your uncommitted work) in the Markdown files under doc/ and artifacts/, and in the .cpp and .hpp
files under src/ and artifacts/. Use `--base master` to check a whole branch. Old prose is not
checked, since much of it predates the rules.

Each rule is a regular expression matched against one paragraph of a Markdown file, the comments
on one run of consecutive lines of a .cpp or .hpp file, or the name of one test case, scenario or
section in src/unittest/ or artifacts/unittest/. Lines are joined with spaces, so a phrase broken
across lines is still found. A match prints `path:line: message` and fails the check.

Unit tests exercise library code, so their comments and names follow the library's rule: they
describe a setting, and do not name the option that sets it.

The option check finds long options such as --min-mapq only. In a comment, a short option such
as -L looks like a minus sign, another program's flag (grep -v) or part of a diagram, so the
writer and the reviewer check for short options themselves.
"""

import argparse
import os
import re
import subprocess
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

# Rules checked one comment line, paragraph or test name at a time.
LINE_RULES = [
    # A thing defined by likeness to another ("oriented like the start and end of a Snarl"). Say
    # what the thing is.
    (re.compile(r"\b(oriented|ordered|sorted|spelled|numbered|placed|named)\s+(just\s+)?like\b",
                re.IGNORECASE),
     "defines a thing by likeness to another; say what it is"),
]
# A command-line option, such as --read-likelihood. Only the subcommands in src/subcommand/ have
# options, so only their comments name them. Library code and its unit tests describe a setting
# by what it means.
OPTION = re.compile(r'(?<![\w-])--[a-z][a-z0-9]*(?:-[a-z0-9]+)*')
OPTION_MESSAGE = "names a command-line option outside src/subcommand/; describe the setting instead"
TEST_OPTION_MESSAGE = ("names a command-line option in a unit test name; "
                       "name the library setting or describe the behaviour")
SUBCOMMAND = os.path.join('src', 'subcommand') + os.sep
UNITTEST = (os.path.join('src', 'unittest') + os.sep, os.path.join('artifacts', 'unittest') + os.sep)
# The files checked without paths given: documentation, and C++ source.
CHECKED = ['doc/*.md', 'artifacts/*.md', 'src/*.cpp', 'src/*.hpp', 'artifacts/*.cpp',
           'artifacts/*.hpp']
# The top of the repository, which the checked paths are relative to.
ROOT = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True,
                      check=True).stdout.strip()
# The header of a diff hunk, whose first group is the start and length of its lines in the new file.
HUNK = re.compile(r'^@@ -\S+ \+(\S+) @@')
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
        if path.endswith('.md'):
            for pattern, message in LINE_RULES:
                for match in pattern.finditer(para):
                    yield max(n for n, o in zip(lines, offsets) if o <= match.start()), message
    if not path.endswith(('.cpp', '.hpp')):
        return
    for line, comment in comment_lines(text):
        for pattern, message in LINE_RULES:
            if pattern.search(comment):
                yield line, message
        option = OPTION.search(comment)
        if option and not path.startswith(SUBCOMMAND):
            yield line, f"{OPTION_MESSAGE}: {option.group()}"
    if not path.startswith(UNITTEST):
        return
    for line, name in test_names(text):
        for pattern, message in RULES:
            if re.search(pattern, name, re.IGNORECASE):
                yield line, message
        for pattern, message in LINE_RULES:
            if pattern.search(name):
                yield line, message
        option = OPTION.search(name)
        if option:
            yield line, f"{TEST_OPTION_MESSAGE}: {option.group()}"


def added_lines(base):
    """Return, for each checked file that changed since `base`, the set of line numbers it adds."""
    result = subprocess.run(['git', 'diff', '-U0', '--no-color', base, '--', *CHECKED],
                            stdout=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit(result.returncode)
    diff = result.stdout
    added, path = {}, None
    for line in diff.split('\n'):
        if line.startswith('+++ '):
            path = line[len('+++ b/'):] if line != '+++ /dev/null' else None
        elif line.startswith('@@') and path is not None:
            start, _, count = HUNK.match(line).group(1).partition(',')
            first, length = int(start), int(count) if count else 1
            added.setdefault(path, set()).update(range(first, first + length))
    return added


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--base', default='HEAD',
                        help='check the lines added since this commit (default: HEAD)')
    parser.add_argument('paths', nargs='*', help='check these files in full instead')
    args = parser.parse_args()
    if args.paths:
        paths = [os.path.relpath(os.path.abspath(p), ROOT) for p in args.paths]
    os.chdir(ROOT)
    lines_to_check = None if args.paths else added_lines(args.base)
    if lines_to_check is not None:
        paths = sorted(lines_to_check)
    ok = True
    for path in paths:
        with open(path, encoding='utf-8', errors='replace') as f:
            text = f.read()
        for line, message in sorted(problems(path, text), key=lambda problem: problem[0]):
            if lines_to_check is None or line in lines_to_check[path]:
                print(f"{path}:{line}: {message}")
                ok = False
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

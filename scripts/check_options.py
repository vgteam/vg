#!/usr/bin/env python3
"""Check that command line options are correctly registered.

Run with scripts/check_options.py

## Summary

Reads the options within the helptext (help_<command>() function),
the long_options[] array, the getopt_long() string, and the switch(c)
block, and checks that they are consistent with each other.

There are a lot of rules; this script has *strong* opinions about how
options should appear, though I've been flexible about valid
variations. If something looks off, a detailed message will be
printed to stdout. If that message doesn't make sense, something
might've been ignored due to a wrong format, so check "Format".
If it still doesn't make sense, ping @faithokamoto
(if I'm still here) and/or open an issue on GitHub.

## Checks

### File as a whole

The file is processed four times in an attempt to extract
the four option sets: helptext, long_options[], getopt_long() string,
and switch(c) block. If any of these processes run into completely
unusable formatting, such as some cases of duplicate options that
simply can't be stored in my data structures, a ValueError is raised
and printed to stdout. No other checks are performed.

If you get one of these errors, then, you might have to
fix multiple things in the file before the script can run well.

Checks performed on the *file as a whole*.

1. The help options (-h, -?) must be present in the helptext,
   getopt string, and the switch block, and must not take an
   argument. In the switch block, they must exit/return/abort.
   (It has at least one of the keywords in `NON_WORK_WORDS`)
2. The default case in the switch block must exit/return/abort.
3. All options in the getopt string must be in long_options[].
4. All options in the switch block must be in long_options[].

While long_options[], the getopt string, and the switch block
must have the same option set, the helptext is allowed to be
missing some options (e.g. deprecated/developer-only options).

Also -? should appear only in the getopt string and the switch block,
and not in long_options[] or the helptext.

### Each option separately

Longform options are pulled from the helptext and long_options[]

1. All longform options must be in long_options[]
2. The shortform option must be in the getopt string
   (if not an ALL_CAPS variable name), and the switch block.
3. The long_options[] entry must be either `no_argument`
   or `required_argument`.
4. The long_options[] entry must match the helptext entry, both
    in whether it takes an argument and in its shortform, if the
    latter is a single character (not ALL_CAPS variable name).
5. If the longform option has a single-character shortform,
   (i.e. not an ALL_CAPS variable name), it must be in the
   getopt string, which must correctly indicate whether it
   takes an argument via using a trailing colon.
6. The long_options[] entry must match the switch block entry,
   in whether it takes an argument.

Note that only the first of these which fails will be printed,
so if you see a message about a missing long_options[] entry,
you may have to do multiple runs/fixes to see all the problems.

## Format

- helptext: within the `help_<command>()` function,
  options must be printed with `"  -<short>, --<long> <arg>  <desc>"`
  (the description is ignored, and the shortform and argument are optional).
  The longform option must be composed of alphanumeric characters
  and hyphens, and the argument must be in all-caps.
  In addition, there must be a "usage:" line somewhere.
  
  Option helptext is not allowed to be over 80 characters long,
  and all descriptions must line up properly. If there are more
  than 80 characters between the `"`s, but that's because you
  have some long default value, just stick that bit on a new line
  of code without changing the printed text. Check out
  `filter_main.cpp`'s `--batch-size` for an example.

- long_options[]: within the `long_options[]` array,
  must be an array of `{"longform", arg_type, 0, shortform}`
  where `arg_type` is either `no_argument` or `required_argument`.
  The longform option must be a string, and the shortform
  must be a single character (NOT a bare non-quoted integer)
  or an ALL_CAPS variable name. Note that ALL_CAPS variables
  must be before the long_options[] array as constexpr ints.

  Shortforms may repeat, and the first longform is retained.
  This is to allow longform aliases.

- getopt_long() string: within the short option string,
  a string of single-character shortform options.
  If and only if a shortform option takes an argument,
  it must be followed by a colon.

- switch(c) block: within the switch block handling options,
  each case must be of the form `case <shortform>:` or `default:`.
  The shortform must be a single character or an ALL_CAPS variable name.
  If the case takes an argument, it must use `optarg`.
  If it *intentionally* doesn't work (e.g. is deprecated, calls `exit`)
  it may or may not use `optarg`. Fallthroughs are handled.

  An attempt is made to enforce the usage of `require_exists()`,
  `ensure_writable()`, and `set_thread_count()`.

For all checks, commented-out lines are ignored.
(Though multiline comments aren't handled correctly.)

In addition, the logic to figure out when we're inside of a block
marked by {} assumes that there is only one curly brace of each
type maximum per line. }} might compile but this script will complain.

## Attribution

The base of this script was written by ChatGPT.
It has since been developed by @faithokamoto.
Some GitHub Copilot autocompletions were used.
"""

import os
import re
import sys
from typing import Dict, List, Optional
from dataclasses import dataclass

SUBCOMMAND_DIR = 'src/subcommand'
"""Where to search for subcommand files."""
SKIP_FILES = {'test_main.cpp', 'help_main.cpp'}
"""Files to skip in the consistency check."""
HELPTEXT_FUNCTION = r'void\shelp_\w+\s*\(char\*\* argv.*\) \{'
"""Regex to match the start of a helptext function."""
HELP_DESC = 'print this help message to stderr and exit'
"""Expected description for the --help option."""
ANNOTATE_EXCEPTIONS = {'xg-name', 'bed-name'}
"""annotate_main.cpp lets these appear twice in helptext."""
NON_WORK_WORDS = {'exit', 'return', 'deprecated', 'abort', 'throw', 'error'}
"""Keywords that indicate something is meant to not work in the switch block."""
FILENAME_VAR_ENDS = {'file', 'filename', 'file_name', 'filepath', 'file_path'}
"""Suffixes that indicate a variable is definitely a filename"""

@dataclass
class OptionInfo:
    """Information about a command line option.
    
    Doesn't store a longform option because it expects
    to be a value in a {longform: OptionInfo} dictionary.
    (Though the switch block parsers uses {shortform: OptionInfo})
    """

    takes_argument: bool = None
    """Whether the option takes an argument."""
    shortform: Optional[str] = None
    """Short option name, or None if not present."""
    errors: List[str] = None
    """Extra things to complain about at the end."""

    def is_unset(self) -> bool:
        """Check if this option is unset."""
        return self.takes_argument is None and self.shortform is None
    
    def is_compatible(self, other) -> bool:
        """Check if two OptionInfo objects are compatible.
        
        Whether they take an argument must match.
        The `other` is assumed to have a shortform, and
        if that shortform is a char (i.e. not an ALL_CAPS variable),
        then it must match the shortform of this object.

        Raises
        ------
        ValueError
            If `other` has a None shortform, which is not allowed.
        """

        if not isinstance(other, OptionInfo):
            return False
        if other.shortform is None:
            raise ValueError("`other` OptionInfo cannot have None shortform; "
                             "please guarantee that it has a shortform")
        if self.shortform is None and len(other.shortform) != 1:
            # If own shortform is unset and other's shortform is long
            # probably due to being all-caps, only compare takes_argument
            return self.takes_argument == other.takes_argument
        else:
            return (self.takes_argument == other.takes_argument
                    and self.shortform == other.shortform)

    def __str__(self) -> str:
        """String representation of this OptionInfo."""
        if self.is_unset():
            return "no info"
        
        if self.shortform is None:
            shortform_str = "no shortform"
        else:
            # If shortform is a single character, print it
            # Otherwise, print it as a string
            if len(self.shortform) == 1:
                shortform_str = f"shortform -{self.shortform}"
            else:
                shortform_str = f'shortform "{self.shortform}"'

        if self.errors is None or not self.errors:
            error_string = "no errors"
        else:
            error_string = "; ".join(self.errors)

        return (f"{shortform_str}, "
                + ("does not take" if not self.takes_argument else "takes")
                + f" an argument; errors: {error_string}")

def extract_help_options(text: str) -> Dict[str, OptionInfo]:
    """Extract options from help_<command>()
    
    Looks within lines the helptext function,
    starting with a line with `void help_<command>(`
    and ending with the outermost closing curly brace `}`.
    (i.e. nested curly braces are handled correctly)

    Ignores lines with single-line comments.

    Looks for `"  -<short>, --<long> <arg>  <desc>"`
    The shortform option is optional, as is the argument.

    If the shortform option is present, it must be followed by
    a comma and a single space before the longform option.
    
    If the argument is present, it is expected to be in all-caps,
    have nothing surrounding it (e.g. no <>), and be exactly
    one space after the long option.

    Helptext lines should have all the options line up properly.

    Parameters
    ----------
    text : str
        The text of the file to search for help options.
    
    Returns
    -------
    Dict[str, OptionInfo]
        A dictionary mapping long option names to OptionInfo.
    
    Raises
    ------
    ValueError
        Gross formatting issues, such as duplicate longform options.
    """

    if not re.search(HELPTEXT_FUNCTION, text):
        raise ValueError("Helptext function not found in file")

    help_opts = dict()
    has_usage = False

    # Match the line's prefix, which should be two spaces
    prefix = r'\s+'
    # Match a shortform option: `-<short>`
    shortform_patten = r'-[^\s]'
    # Match the space between the forms: `, `
    between_forms_pattern = r',*\s+'
    # Match a longform option: `--<long>`
    longform_patten = r'--[a-zA-Z0-9\-]+'
    # Match an optional argument in all-caps, with a preceeding space
    arg_pattern = r'\s[A-Z_\-=:.,;\[\]]+'
    # Match the rest of the description
    desc_pattern = r'\s+.+'
    # Match what is printed
    quoted_pattern = r'"(.+?)"'

    # Certain exceptions can only be used in annotate_main.cpp
    is_annotate = 'help_annotate' in text
    # giraffe_main.cpp has a different --help description
    is_giraffe = 'help_giraffe' in text

    help_pattern = re.compile(
        rf'"({prefix})({shortform_patten})({between_forms_pattern})'
        rf'({longform_patten})({arg_pattern})?({desc_pattern})"'
        )
    # Same pattern but without the shortform option
    long_only_pattern = re.compile(
        rf'"({prefix})({longform_patten})({arg_pattern})?({desc_pattern})"'
        )

    inside_help = False
    curly_brace_nesting = 0
    description_start = None

    def save_option(match: re.Match, has_shortform: bool) -> None:
        """Helper function to save an option."""
        errors = []
        nonlocal description_start

        prefix = match.group(1)
        shortform = match.group(2)[1:] if has_shortform else None
        between = match.group(3) if has_shortform else None
        longform = match.group(4 if has_shortform else 2)[2:]
        takes_arg = match.group(5 if has_shortform else 3) is not None
        description_num = 6 if has_shortform else 4
        description = match.group(description_num)

        if len(prefix) != (6 if shortform is None else 2):
                errors.append("Shortforms must have a two-space indent; "
                              "longforms must have a six-space indent")
        if between is not None:
            if between[0] != ',':
                errors.append(f"There must be a comma between "
                              "shortform and longform")
                spaces = len(between)
            else:
                spaces = len(between) - 1 # Ignore comma
            if spaces != 1:
                errors.append(f"There must be one space after the comma "
                              "between shortform and longform")
        if not description.startswith('  '):
            errors.append(f"There must be at least two spaces between"
                          "option and description")
        if (longform == 'help' and description.strip() != HELP_DESC
            and not is_giraffe):
            errors.append(f"--help must have the description '{HELP_DESC}'")
        
        # Count of blank spaces between option and description
        text_offset = re.match(r'\s+(.+)', description).start(1) - match.start()
        # Column within printed text which starts the description
        text_start = match.start(description_num) + text_offset
        # Use first line's text_start as the description start
        if description_start is None:
            description_start = text_start
        # All other lines must match
        elif description_start != text_start:
            errors.append(f"Help option line '{stripped}' does not "
                          "line up with the previous option's description")
            
        if longform not in help_opts:
            help_opts[longform] = OptionInfo(takes_arg, shortform, errors)
        elif not (is_annotate and longform in ANNOTATE_EXCEPTIONS):
            raise ValueError(f"Duplicate longform option '{longform}' found in "
                             "helptext")

    lines = text.splitlines()
    for i in range(len(lines)):
        stripped = lines[i].strip()
        # Are we inside the helptext printing function?
        if re.search(HELPTEXT_FUNCTION, stripped):
            inside_help = True
            curly_brace_nesting = -1
        # End at end of outermost curly braces
        elif inside_help and stripped == '}' and curly_brace_nesting == 0:
            break

        # Ignore comments and lines outside the helptext
        # Multiline comments are not handled correctly,
        # but they shouldn't be here anyway.
        if not inside_help or stripped.startswith('//'):
            continue

        # Keep track of curly brace nesting
        if '{' in stripped:
            curly_brace_nesting += 1
        if '}' in stripped:
            curly_brace_nesting -= 1

        if '"usage:' in stripped:
            if has_usage:
                raise ValueError("Multiple 'usage:' lines found in helptext")
            has_usage = True
            # Don't check usage lines for length
            continue

        printed_line = re.search(quoted_pattern, stripped)
        if printed_line:
            printed_len = printed_line.end(1) - printed_line.start(1)
            if printed_len > 80:
                raise ValueError(f"line {i+1} has helptext string `{stripped}` "
                                 f"which is {printed_len}>80 characters long")
            if r"\n" in printed_line.group(1):
                raise ValueError(f"line {i+1} has helptext string `{stripped}` "
                                 "with a newline character; use endl instead")

        # Parse out sections of full pattern
        match = help_pattern.search(stripped)
        if match:
            save_option(match, True)
            continue

        # Parse out sections of long-only pattern
        match = long_only_pattern.search(stripped)
        if match:
            save_option(match, False)

    if not has_usage:
        raise ValueError("Helptext is missing a 'usage:' line")
    return help_opts

def extract_long_options(text: str) -> Dict[str, OptionInfo]:
    """Extract options from long_options[].

    Looks within the long_options[] array,
    starting with a line with `struct option long_options`
    or `std::vector<struct option> long_options`
    and ending with `};`. (i.e. assuming no nested curly braces)

    Ignores lines with comments.

    Looks for: `{"longform", arg_type, 0, shortform}`
    Ignores all-zeros, which tends to end long_options[].
    Smart enough to ignore comments at the end of a line.

    The argument type is either `no_argument` (False)
    or `required_argument` (True). If the argument type is neither,
    it is set to None, marking an unknown argument type.

    Note that the shortform can be a number, which is
    preserved as an integer, or an all-caps variable name.
    All-caps variable name must be defined as constexpr ints.

    Shortforms may repeat, and the first longform is retained.

    Parameters
    ----------
    text : str
        The text of the file to search for long_options.

    Returns
    -------
    Dict[str, OptionInfo]
        A dictionary mapping long option names to OptionInfo.
    
    Raises
    ------
    ValueError
        Gross formatting issues, such as duplicate longform options.
    """

    options = dict()
    all_shortforms = set()
    inside_longopts = False
    
    # Track constexpr int definitions
    const_pattern = re.compile(r'constexpr int ([A-Z_]+)\s+=\s+(\d+);')
    all_caps_names = set()
    all_caps_values = set()

    for line in text.splitlines():
        # Found start of long_options[]
        if ('struct option long_options' in line
            or 'std::vector<struct option> long_options' in line):
            inside_longopts = True
        # End of long_options[]
        elif inside_longopts and '};' in line:
            break
        
        match = const_pattern.search(line)
        if match:
            # This is now a legal all-caps variable name
            all_caps_names.add(match.group(1))
            cur_value = match.group(2)
            if cur_value in all_caps_values:
                raise ValueError(f"Duplicate constexpr int value '{cur_value}'")
            all_caps_values.add(cur_value)

        if not inside_longopts:
            continue

        stripped = line.strip()
        if stripped.startswith('{') and not stripped.startswith('{0'):
            parts = stripped.split('//')[0].strip('{} \t,\n').split(',')

            if len(parts) == 4:
                longform = parts[0].strip().strip('"')
                arg_type = parts[1].strip()
                shortform = parts[3].strip()
                errors = []

                if longform in options:
                    raise ValueError(f"Duplicate longform option '{longform}' "
                                     "found in long_options[]")

                try:
                    # Check if shortform is a number
                    shortform = int(shortform)
                    errors.append(f"--{longform} has int shortform {shortform} "
                                  "in long_options[]; use a char or ALL_CAPS "
                                  "variable instead")
                except ValueError:
                    # Otherwise, it is a character or string
                    shortform = shortform.strip("'")
                    if len(shortform) > 1 and not shortform in all_caps_names:
                        errors.append(f"--{longform} has non-char shortform "
                                      f"'{shortform}' in long_options[] "
                                      " which isn't a constexpr int")
                
                # Skip duplicate shortforms (aliases)
                if shortform in all_shortforms:
                    continue

                # Ignore placeholder line with all zeros
                if longform == '0':
                    continue

                if arg_type == 'no_argument':
                    takes_arg = False
                elif arg_type == 'required_argument':
                    takes_arg = True
                else:
                    # Marker for an unknown argument type
                    takes_arg = None

                options[longform] = OptionInfo(takes_arg, shortform, errors)
                all_shortforms.add(shortform)
    return options

def extract_getopt_string(text: str) -> Dict[str, bool]:
    """Extract short options from getopt_long() string.

    Looks for the line with `getopt_long()` or `short_options =`.
    For the former, take the third argument.
    For the latter, take the string assigned.

    Shortform options followed by : take an argument.
    If an option appears multiple times, it is assigned None
    as a marker value.

    Parameters
    ----------
    text : str
        The text of the file to search for long_options.

    Returns
    -------
    Dict[str, bool]
        A dictionary mapping shortform names to
        whether they take an argument.
    """

    # Grab third argument of getopt_long()
    match = re.search(r'getopt_long\s*\([^,]+,[^,]+,\s*"([^"]+)"', text)
    if match:
        opts_str = match.group(1)

    # Second try to catch vg giraffe/augment weirdness
    if not match:
        match = re.search(r'short_options = "([^"]+)"', text)
        if match:
            opts_str = match.group(1)

    # Give up
    if not match:
        return dict()
    
    options = dict()
    i = 0
    while i < len(opts_str):
        shortform = opts_str[i]
        next_is_colon = (i+1 < len(opts_str) and opts_str[i+1] == ':')

        # Is this a duplicate?
        if shortform in options:
            options[shortform] = None
            i += 2 if next_is_colon else 1
            continue

        # Add shortform with its `:`
        if next_is_colon:
            options[shortform] = True
            # Skip `:`
            i += 2
        # This shortform option has no `:`
        else:
            options[shortform] = False
            i += 1
    return options

def extract_switch_optarg(text: str) -> Dict[str, OptionInfo]:
    """Compile optarg usage within switch(c) block.

    Looks within the switch block handling options,
    starting with a line with `switch(c)` (plus optional spaces)
    and ending with the outermost closing curly brace `}`.
    (i.e. nested curly braces are handled correctly)

    Ignores lines with single-line comments. Handles fallthroughs,
    breaking cases on either a `break;` or a keyword (see
    `NON_WORK_WORDS`) which is meant to not work.

    Looks for whether each case uses `optarg` or not.
    Things that don't work are marked with None.

    Parameters
    ----------
    text : str
        The text of the file to search for long_options.

    Returns
    -------
    Dict[str, OptionInfo]
        A dictionary mapping shortform names to a more
        detailed OptionInfo struct

    Raises
    ------
    ValueError
        Gross formatting issues, such as duplicate shortform options.
    """

    optarg_usage = dict()
    extra_errors = []
    shortforms = set()
    inside_switch = False
    # Current cases being processed
    current_cases = []
    has_optarg = False

    def save_case_info(force_none: bool = False) -> None:
        """Set optarg usage & errors for a list of cases."""
        # Look up the current state that we're saving
        nonlocal current_cases
        nonlocal extra_errors
        nonlocal has_optarg
        nonlocal optarg_usage

        # What usage state are we saving?
        usage = None if force_none else has_optarg

        for case in current_cases:
            if case in shortforms:
                raise ValueError(f"Duplicate shortform option '{case}' found "
                                 "in switch block")
            shortforms.add(case)
            optarg_usage[case] = OptionInfo(usage, case, extra_errors)

        # Reset the current cases, optarg usage, and errors
        has_optarg = False
        current_cases = []
        extra_errors = []

    for line in text.splitlines():
        stripped = line.strip()

        # Are we inside the switch statement?
        if re.search(r'switch\s*\(\s*c\s*\)', line):
            inside_switch = True
            curly_brace_nesting = -1
        elif inside_switch and line.strip() == '}' and curly_brace_nesting == 0:
            break

        # Ignore comments and lines outside the switch
        # Multiline comments are not handled correctly,
        # but they shouldn't be here anyway.
        if not inside_switch or line.strip().startswith('//'):
            continue
        if '{' in stripped:
            curly_brace_nesting += 1
        if '}' in stripped:
            curly_brace_nesting -= 1

        if 'omp_set_num_threads' in stripped:
            # Extra check for thread count options
            extra_errors.append("Parse thread count using set_thread_count() "
                                "for standardized error messages")
    
        # This won't catch all file variables (e.g. if called `xg_name`),
        # but it'll catch some at least
        for suffix in FILENAME_VAR_ENDS:
            if f'{suffix} = optarg;' in stripped:
                # Extra check for file-existance functions
                extra_errors.append("Use require_exists() or ensure_writable() "
                                    "for standardized file checks: " + stripped)
                break

        # Detect new case
        case_match = re.match(r'case\s+(.+)\s*:', stripped)
        if case_match or stripped == 'default:':
            # Flush old cases to avoid optarg
            # being present after fallthrough
            if has_optarg:
                save_case_info()
            
            # Remove end-of-line comments
            stripped = stripped.split('//')[0].strip() 
            case_value = case_match.group(1) if case_match else 'default'

            if case_value.startswith("'") and case_value.endswith("'"):
                # Single-character variable
                case_key = case_value.strip("'")
            else:
                # All-caps variable
                case_key = case_value
            current_cases.append(case_key)
            continue

        # If it's not a new case, it belongs to the current block
        has_optarg = has_optarg or 'optarg' in stripped

        # On intentional non-working cases, flush the current case group
        if (any(keyword in stripped for keyword in NON_WORK_WORDS)
            and curly_brace_nesting == 0):
            # None is a marker for not working
            save_case_info(force_none=True)

        # On break, flush the current case group
        if stripped == 'break;' and curly_brace_nesting == 0:
            save_case_info()

    # Handle trailing block without break (not common but valid)
    if current_cases:
        save_case_info()

    return optarg_usage

def check_file(filepath: str) -> bool:
    """Run all consistency checks on a single file.
    
    Any problems are printed to stdout.
    See file docstring for details.

    Parameters
    ----------
    filepath : str
        The path to the file to check.

    Returns
    -------
    bool
        True if the file is OK, and False if it
        contains problems.
    """

    is_ok = True
    def problem(description: str) -> None:
        """
        Log a problem and remember that there was one.

        Parameters
        ----------
        description : str
            The description of the problem to diplay.
        """
        nonlocal is_ok
        is_ok = False
        print(f"{filepath}: {description}")

    # Save file contents
    with open(filepath) as f:
        text = f.read()

    # Process whole file four times to look up four option sets
    try:
        help_opts = extract_help_options(text)
        long_opts = extract_long_options(text)
        getopt_opts = extract_getopt_string(text)
        switch_opts = extract_switch_optarg(text)
    except ValueError as e:
        problem(f"{e}")
        return is_ok

    all_longform = set(help_opts) | set(long_opts)

    # Overall check 1: are help options present?
    for help_alias in ['?', 'h']:
        if help_alias not in getopt_opts:
            problem(f"help alias -{help_alias} is missing from "
                    "getopt string")
        elif getopt_opts[help_alias]:
            problem(f"help alias -{help_alias} should not take an "
                    "argument")
        
        if help_alias not in switch_opts:
            problem(f"help alias -{help_alias} is missing from "
                    "switch block")
        elif switch_opts[help_alias].takes_argument is not None:
            problem(f"help alias -{help_alias} should exit/return/abort "
                    "switch block")
        
        # Change help alias to be non-argument
        # in order to match long_options[]
        switch_opts[help_alias] = OptionInfo(False, help_alias, [])
    
    if 'help' not in help_opts:
        problem("help option --help is missing from helptext")
    elif help_opts['help'].takes_argument:
        problem("help option --help should not take an argument")
    elif help_opts['help'].shortform != 'h':
        problem("help option --help should have shortform -h")
    
    # Overall check 2: does the default case fail to work?
    if 'default' not in switch_opts:
        problem("switch block is missing a default case")
    else:
        if switch_opts['default'].takes_argument is not None:
            problem("switch block's default case should exit/return/abort")
        
        # Get rid of default once it's checked; will appear nowhere else
        switch_opts.pop('default')

    # Check longform options between the four sources
    for longform in all_longform:
        # Look up the longform option in the helptext, if it is there
        if longform in help_opts:
            cur_help = help_opts[longform]
            for cur_error in cur_help.errors:
                problem(f"{filepath} has error in helptext for --{longform}: "
                        f"{cur_error}")
        else:
            cur_help = OptionInfo()

        # Check 1: all options must appear in long_options[] with shortform
        if longform not in long_opts:
            problem(f"--{longform} is not in long_options[]")
            continue

        # Get the long_options[] entry, which is treated as truth
        cur_long = long_opts[longform]
        should_str = 'should' if cur_long.takes_argument else "shouldn't"
        
        # Check 2: this option appears in the getopt string & switch block
        if len(cur_long.shortform) == 1 :
            if cur_long.shortform not in getopt_opts:
                problem(f"--{longform}'s -{cur_long.shortform} "
                        f"should be in getopt string")
                continue
            else:
                cur_getopt = getopt_opts.pop(cur_long.shortform)
        if cur_long.shortform not in switch_opts:
            problem(f"--{longform}'s -{cur_long.shortform} is "
                    "not used in the switch block")
            continue
        else:
            cur_switch = switch_opts.pop(cur_long.shortform)
            for e in cur_switch.errors:
                problem(f"{e}")

        # Check 3: long_options uses only no_argument or required_argument
        if cur_long.takes_argument is None:
            problem(f"--{longform} has an unknown argument type "
                    "in long_options[]; use no_argument or required_argument")
            continue

        # Any other errors associated with this long_options entry,
        # before we compare it to the other sections?
        if cur_long.errors:
            for e in cur_long.errors:
                problem(f"{e}")
            continue

        # Check 4: help vs long_options[]
        if not cur_help.is_unset() and not cur_help.is_compatible(cur_long):
            problem(f"--{longform} has mismatch between helptext "
                    f"({cur_help}) and long_options[] ({cur_long})")
            continue

        # Check 5: long_options[] vs getopt string
        # This check only runs for single-char shortforms
        if len(cur_long.shortform) == 1:
            if cur_getopt is None:
                problem(f"--{longform}'s -{cur_long.shortform} "
                        "appears multiple times in getopt string")
                continue

            if cur_long.takes_argument != cur_getopt:
                problem(f"--{longform}'s -{cur_long.shortform} "
                        f"{should_str} have a : after it in getopt string")
                continue

        # Check 6: long_options[] vs switch
        if (cur_switch.takes_argument is not None 
            and cur_long.takes_argument != cur_switch.takes_argument):
            problem(f"--{longform}'s -{cur_long.shortform} "
                    f"{should_str} use optarg")
            continue
    
    # Overall check 3: getopt is allowed to have ?
    # without that being in long_options[]
    if '?' in getopt_opts:
        getopt_opts.pop('?')
    if getopt_opts:
        problem("getopt string has option(s) not in long_options[]: "
                f"{', '.join(getopt_opts)}")
    
    # Overall check 4: similarly for switch block
    if '?' in switch_opts:
        switch_opts.pop('?')
    if switch_opts:
        problem("switch block has option(s) not in long_options[]: "
                f"{', '.join(switch_opts)}")

    return is_ok

if __name__ == "__main__":
    is_ok = True
    for fname in os.listdir(SUBCOMMAND_DIR):
        if not fname.endswith('_main.cpp') or fname in SKIP_FILES:
            continue
        is_ok = check_file(os.path.join(SUBCOMMAND_DIR, fname)) and is_ok
    sys.exit(0 if is_ok else 1)

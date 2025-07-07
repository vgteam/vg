"""Check that command line options are correctly registered.

## Summary

Reads the options within the helptext (help_<command>() function),
the long_options[] array, the getopt_long() string, and the switch(c)
block, and checks that they are consistent with each other.

There are a lot of rules; this script has *strong* opinions about how
options should appear, though I've been flexible about valid
variations. If something looks off, a detailed message will be
printed to stdout. If that message doesn't make sense, something
might've been ignored due to a wrong format, so check "Format".
If it still doesn't make sense, ping me.

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

1. The help options (-h, -?) must be present in the getopt string
   and the switch block, and must not take an argument. In the
   switch block, they must crash.
2. The default case in the switch block must crash.
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
2. The long_options[] entry must be either `no_argument`
   or `required_argument`.
3. The long_options[] entry must not have a numeric shortform.
4. The long_options[] entry must match the helptext entry, both
    in whether it takes an argument and in its shortform.
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

- long_options[]: within the `long_options[]` array,
  must be an array of `{"longform", arg_type, 0, shortform}`
  where `arg_type` is either `no_argument` or `required_argument`.
  The longform option must be a string, and the shortform
  must be a single character or an ALL_CAPS variable name.
  Shortforms may repeat, and the first longform is retained.

- getopt_long() string: within the short option string,
  a string of single-character shortform options.
  If and only if a shortform option takes an argument,
  it must be followed by a colon.

- switch(c) block: within the switch block handling options,
  each case must be of the form `case <shortform>:` or `default:`.
  The shortform must be a single character or an ALL_CAPS variable name.
  If the case takes an argument, it must use `optarg`.
  If a case crashes or is otherwise deprecated, it may or may not
  use `optarg`; I'm flexible about that. Fallthroughs are handled.

For all checks, commented-out lines are ignored.
(Though multiline comments aren't handled correctly.)

## TODO

- Check that all helptext lines are no more than 80 characters
- Check that all helptext descriptions start at the same column
- Check that the second line of getopt_long()'s arguments matches up
- Check that the order of options in the helptext
  matches the order in long_options[] and getopt_long() string
- Require "usage" line in helptext
- Require "options" line in helptext
- Require helptext to have a "help" option

## Attribution

The base of this script was written by ChatGPT.
It has since been developed by @faithokamoto.
Some GitHub Copilot autocompletions were used.
"""

import os
import re
from typing import Dict, List, Optional
from dataclasses import dataclass

SUBCOMMAND_DIR = 'src/subcommand'
"""Where to search for subcommand files."""
SKIP_FILES = {'test_main.cpp', 'help_main.cpp'}
"""Files to skip in the consistency check."""
ANNOTATE_EXCEPTIONS = {'xg-name', 'bed-name'}
"""annotate_main.cpp lets these appear twice in helptext."""
# Some giraffe arguments are processed outside of the switch block
# (at least I think so? I'm not entirely sure)
GIRAFFE_EXCEPTIONS = {'max-multimaps', 'batch-size'}
"""Options for giraffe_main.cpp that are not in the switch block."""
CRASHES = {'exit', 'return', 'deprecated', 'abort', 'throw'}
"""Keywords that indicate a crash/deprecation in the switch block."""

@dataclass
class OptionInfo:
    """Information about a command line option."""

    takes_argument: bool = None
    """Whether the option takes an argument."""
    shortform: Optional[str | int] = None
    """Short option name, or None if not present."""
    errors: List[str] = None
    """Some special errors to print about this option."""

    def is_unset(self) -> bool:
        """Check if this option is unset."""
        return self.takes_argument is None and self.shortform is None
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, OptionInfo):
            return False
        # If either shortform is unset, only compare takes_argument
        if self.shortform is None or other.shortform is None:
            return self.takes_argument == other.takes_argument
        else:
            return (self.takes_argument == other.takes_argument
                    and self.shortform == other.shortform)

def extract_help_options(text: str) -> Dict[str, OptionInfo]:
    """Extract options from help_<command>()
    
    Looks within lines the helptext function,
    starting with a line with `void help_<command>(`
    and ending with the outermost closing curly brace `}`.
    (i.e. nested curly braces are handled correctly)

    Ignores lines with comments.

    Looks for `"  -<short>, --<long> <arg>  <desc>"`
    The shortform option is optional, as is the argument.

    If the shortform option is present, it must be followed by
    a comma and a single space before the longform option.
    
    If the argument is present, it is expected to be in all-caps,
    have nothing surrounding it (e.g. no <>), and be exactly
    one space after the long option.

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

    help_opts = dict()
    # Match the line's prefix, which should be two spaces
    prefix = r'\s+'
    # Match a shortform option: `-<short`
    shortform_patten = r'-[\w]'
    # Match a longform option: `--<long>`
    longform_patten = r'--[a-zA-Z0-9\-]+'
    # Match an optional argument in all-caps, with a preceeding space
    arg_pattern = r'\s[A-Z_\-=:.,;\[\]]+'
    # Match the rest of the description
    desc_pattern = r'\s+.+'

    # Certain exceptions can only be used in annotate_main.cpp
    is_annotate = 'help_annotate' in text

    help_pattern = re.compile(
        rf'"({prefix})({shortform_patten}),\s({longform_patten})'
        rf'({arg_pattern})?({desc_pattern})"'
        )
    # Same pattern but without the shortform option
    long_only_pattern = re.compile(
        rf'"({prefix})({longform_patten})({arg_pattern})?({desc_pattern})"'
        )

    inside_help = False
    curly_brace_nesting = 0

    def save_option(shortform: Optional[str | int], 
                    longform: str, takes_arg: bool, errors: List[str]) -> None:
        """Helper function to save an option."""
        if longform not in help_opts:
            help_opts[longform] = OptionInfo(takes_arg, shortform, errors)
        elif not (is_annotate and longform in ANNOTATE_EXCEPTIONS):
            raise ValueError(f"Duplicate longform option '{longform}' found in "
                             "helptext")

    for line in text.splitlines():
        stripped = line.strip()
        # Are we inside the helptext printing function?
        if re.search(r'void\shelp_\w+\s*\(char\*\* argv\) \{', stripped):
            inside_help = True
            curly_brace_nesting = -1
            continue
        # End at end of outermost curly braces
        elif inside_help and stripped == '}' and curly_brace_nesting == 0:
            break

        # Ignore comments and lines outside the helptext
        if not inside_help or stripped.startswith('//'):
            continue

        # Keep track of curly brace nesting
        if '{' in stripped:
            curly_brace_nesting += 1
        if '}' in stripped:
            curly_brace_nesting -= 1

        # Parse out sections of full pattern
        match = help_pattern.search(stripped)
        if match:
            errors = []
            if len(match.group(1)) != 2:
                errors.append(f"Help option line '{stripped}' "
                              "does not start with two spaces")
            if not match.group(5).startswith('  '):
                errors.append(f"Help option line '{stripped}' has less than "
                              "two spaces between option and description")
            save_option(match.group(2)[1], match.group(3)[2:],
                        match.group(4) is not None, errors)
            continue

        # Parse out sections of long-only pattern
        match = long_only_pattern.search(stripped)
        if match:
            errors = []
            if len(match.group(1)) != 6:
                errors.append(f"Longform option {match.group(2)} "
                              "does not line up with the others")
            if not match.group(4).startswith('  '):
                errors.append(f"Help option line '{stripped}' has less than "
                              "two spaces between option and description")
            save_option(None, match.group(2)[2:],
                        match.group(3) is not None, errors)

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

    for line in text.splitlines():
        # Found start of long_options[]
        if ('struct option long_options' in line
            or 'std::vector<struct option> long_options' in line):
            inside_longopts = True
        # End of long_options[]
        elif inside_longopts and '};' in line:
            break

        if not inside_longopts:
            continue

        stripped = line.strip()
        if stripped.startswith('{') and not stripped.startswith('{0'):
            parts = stripped.split('//')[0].strip('{} \t,\n').split(',')

            if len(parts) == 4:
                longform = parts[0].strip().strip('"')
                arg_type = parts[1].strip()
                shortform = parts[3].strip()

                if longform in options:
                    raise ValueError(f"Duplicate longform option '{longform}' "
                                     "found in long_options[]")

                try:
                    # Keep shortform as a number if it is one
                    shortform = int(shortform)
                except ValueError:
                    # Otherwise, it is a character or string
                    shortform = shortform.strip("'")
                
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

                options[longform] = OptionInfo(takes_arg, shortform)
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

def extract_switch_optarg(text: str) -> Dict[str, Optional[bool]]:
    """Compile optarg usage within switch(c) block.

    Looks within the switch block handling options,
    starting with a line with `switch(c)` (plus optional spaces)
    and ending with the outermost closing curly brace `}`.
    (i.e. nested curly braces are handled correctly)

    Ignores lines with comments. Handles fallthroughs,
    breaking cases on either a crash or a `break;`.

    Looks for whether each case uses `optarg` or not.
    Crashes/deprecation are marked with None.

    Parameters
    ----------
    text : str
        The text of the file to search for long_options.

    Returns
    -------
    Dict[str, Optional[bool]]
        A dictionary mapping shortform names to
        whether they use optarg (i.e. take an argument).
        If the case crashes, it is set to None.

    Raises
    ------
    ValueError
        Gross formatting issues, such as duplicate longform options.
    """

    optarg_usage = dict()
    shortforms = set()
    inside_switch = False
    # Current cases being processed
    current_cases = []
    has_optarg = False

    def set_optarg_usage(cases: list[str], usage: Optional[bool]):
        """Set optarg usage for a list of cases."""
        for case in cases:
            if case in shortforms:
                raise ValueError(f"Duplicate shortform option '{case}' found "
                                 "in switch block")
            shortforms.add(case)
            optarg_usage[case] = usage

        # Reset the current cases and optarg usage
        nonlocal has_optarg 
        has_optarg = False
        nonlocal current_cases
        current_cases = []

    for line in text.splitlines():
        stripped = line.strip()

        # Are we inside the switch statement?
        if re.search(r'switch\s*\(\s*c\s*\)', line):
            inside_switch = True
            curly_brace_nesting = -1
        elif inside_switch and line.strip() == '}' and curly_brace_nesting == 0:
            break

        if not inside_switch or line.strip().startswith('//'):
            continue
        if '{' in stripped:
            curly_brace_nesting += 1
        if '}' in stripped:
            curly_brace_nesting -= 1

        # Detect new case
        case_match = re.match(r'case\s+(.+)\s*:', stripped)
        if case_match or stripped == 'default:':
            # Flush old cases to avoid optarg
            # being present after fallthrough
            if has_optarg:
                set_optarg_usage(current_cases, True)
            
            # Remove comments
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

        # On crash/deprecation, flush the current case group
        if (any(keyword in stripped for keyword in CRASHES)
            and curly_brace_nesting == 0):
            # None is a marker for a crash
            set_optarg_usage(current_cases, None)

        # On break, flush the current case group
        if stripped == 'break;' and curly_brace_nesting == 0:
            set_optarg_usage(current_cases, has_optarg)

    # Handle trailing block without break (not common but valid)
    if current_cases:
        set_optarg_usage(current_cases, has_optarg)

    return optarg_usage

def check_file(filepath: str) -> None:
    """Run all consistency checks on a single file.
    
    Any problems are printed to stdout.
    See file docstring for details.

    Parameters
    ----------
    filepath : str
        The path to the file to check.
    """

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
        print(f"{filepath}: {e}")
        return

    all_longform = set(help_opts) | set(long_opts)

    # Overall check 1: are help options present?
    for help_alias in ['?', 'h']:
        if help_alias not in getopt_opts:
            print(f"{filepath}: help alias -{help_alias} is missing from "
                  "getopt string")
        elif getopt_opts[help_alias]:
            print(f"{filepath}: help alias -{help_alias} should not take an "
                  "argument")
        
        if help_alias not in switch_opts:
            print(f"{filepath}: help alias -{help_alias} is missing from "
                  "switch block")
        elif switch_opts[help_alias] is not None:
            print(f"{filepath}: help alias -{help_alias} should crash in "
                  "switch block")
        
        # Change help alias to be non-argument instead of crash
        # in order to match long_options[]
        switch_opts[help_alias] = False
    
    # Overall check 2: does the default case crash?
    if 'default' not in switch_opts:
        print(f"{filepath}: switch block is missing a default case")
    else:
        if switch_opts['default'] is not None:
            print(f"{filepath}: switch block's default case should crash")
        
        # Get rid of default once it's checked; will appear nowhere else
        switch_opts.pop('default')

    # Check longform options between the four sources
    for longform in all_longform:
        # Look up the longform option in the helptext, if it is there
        if longform in help_opts:
            cur_help = help_opts[longform]
            for cur_error in cur_help.errors:
                print(f"{filepath} has error in helptext for --{longform}: "
                      f"{cur_error}")
        else:
            cur_help = OptionInfo()

        # Check 1: all options must appear in long_options[] with shortform
        if longform not in long_opts:
            print(f"{filepath}: --{longform} is not in long_options[]")
            continue

        # Get the long_options[] entry, which is treated as truth
        cur_long = long_opts[longform]
        should_str = 'should' if cur_long.takes_argument else "shouldn't"

        # Check 2: long_options uses only no_argument or required_argument
        if cur_long.takes_argument is None:
            print(f"{filepath}: --{longform} has an unknown argument type "
                  "in long_options[]; use no_argument or required_argument")
            continue

        # Check 3: long_options[] don't use raw numbers
        if isinstance(cur_long.shortform, int):
            print(f"{filepath}: --{longform} has an int ({cur_long.shortform}) "
                  "in long_options[]; use a char or ALL_CAPS variable instead")
            continue

        # Check 4: help vs long_options[]
        if not cur_help.is_unset() and cur_help != cur_long:
            print(f"{filepath}: --{longform} has mismatch between helptext "
                  f"{cur_help} and long_options[] {cur_long}")
            continue

        # Check 5: long_options[] vs getopt string
        # This check only runs for single-char shortforms
        if len(cur_long.shortform) == 1:
            if cur_long.shortform not in getopt_opts:
                print(f"{filepath}: --{longform}'s -{cur_long.shortform} "
                      f"should be in getopt string")
                continue
        
            # Mark that this shortform has been used in getopts
            cur_getopt = getopt_opts.pop(cur_long.shortform)
            
            if cur_getopt is None:
                print(f"{filepath}: --{longform}'s -{cur_long.shortform} "
                      "appears multiple times in getopt string")
                continue

            if cur_long.takes_argument != cur_getopt:
                print(f"{filepath}: --{longform}'s -{cur_long.shortform} "
                      f"{should_str} have a : after it in getopt string")
                continue

        # Check 6: long_options[] vs switch
        if cur_long.shortform not in switch_opts:
            if 'giraffe_main' in filepath and longform in GIRAFFE_EXCEPTIONS:
                # Giraffe exceptions are allowed to not have a switch
                continue
            print(f"{filepath}: --{longform}'s -{cur_long.shortform} is "
                    "not used in the switch block")
            continue

        # Mark that this shortform has been used in the switch block
        cur_switch = switch_opts.pop(cur_long.shortform)

        if cur_switch is not None and cur_long.takes_argument != cur_switch:
            print(f"{filepath}: --{longform}'s -{cur_long.shortform} "
                    f"{should_str} use optarg")
            continue
    
    # Overall check 3: getopt is allowed to have ?
    # without that being in long_options[]
    if '?' in getopt_opts:
        getopt_opts.pop('?')
    if getopt_opts:
        print(f"{filepath}: getopt string has option(s) not in long_options[]: "
              f"{', '.join(getopt_opts)}")
    
    # Overall check 4: similarly for switch block
    if '?' in switch_opts:
        switch_opts.pop('?')
    if switch_opts:
        print(f"{filepath}: switch block has option(s) not in long_options[]: "
              f"{', '.join(switch_opts)}")

if __name__ == "__main__":
    for fname in os.listdir(SUBCOMMAND_DIR):
        if not fname.endswith('_main.cpp') or fname in SKIP_FILES:
            continue
        problems = check_file(os.path.join(SUBCOMMAND_DIR, fname))

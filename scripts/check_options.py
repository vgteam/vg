"""Check that command line options are correctly registered.

Reads the options within the helptext (help_<command>() function),
the long_options[] array, the getopt_long() string, and the switch(c)
block, and checks that they are consistent with each other.

The base of this script was written by ChatGPT.
It has since been developed by @faithokamoto.
"""
import os
import re
from typing import Dict, Optional
from dataclasses import dataclass

SKIP_FILES = {'test_main.cpp', 'help_main.cpp'}
GIRAFFE_EXCEPTIONS = {'max-multimaps', 'batch-size'}

@dataclass
class OptionInfo:
    """Information about a command line option."""

    takes_argument: bool = None
    """Whether the option takes an argument."""
    shortform: Optional[str | int] = None
    """Short option name, or None if not present."""

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

    Looks for `<< "    -<short>, --<long> <arg>  <desc>`
    (well actually it ignores the description)
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
        A dictionary mapping long option names to OptionInfo
    """

    help_opts = dict()
    # Match the line's prefix, which is `<< "    `
    prefix = r'<< "\s+'
    # Match a shortform option: `-<short`
    shortform_patten = r'-.'
    # Match a longform option: `--<long>`
    longform_patten = r'--[a-zA-Z0-9\-]+'
    # Match an optional argument in all-caps
    arg_pattern = r'[A-Z_]+'

    help_pattern = re.compile(
        rf'{prefix}({shortform_patten}),\s({longform_patten})\s({arg_pattern})?'
        )
    # Same pattern but without the shortform option
    long_only_pattern = re.compile(
        rf'{prefix}({longform_patten})\s({arg_pattern})?'
        )

    inside_help = False
    curly_brace_nesting = 0

    for line in text.splitlines():
        stripped = line.strip()
        # Are we inside the helptext printing function?
        if re.search(r'void\shelp_\w+\s*\(', line):
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
            shortform = match.group(1)[1]
            longform = match.group(2)[2:]
            takes_arg = match.group(3) is not None
            help_opts[longform] = OptionInfo(takes_arg, shortform)
            continue

        # Parse out sections of long-only pattern
        match = long_only_pattern.search(stripped)
        if match:
            longform = match.group(1)[2:]
            takes_arg = match.group(2) is not None
            help_opts[longform] = OptionInfo(takes_arg)

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

    Parameters
    ----------
    text : str
        The text of the file to search for long_options.

    Returns
    -------
    Dict[str, OptionInfo]
        A dictionary mapping long option names to OptionInfo
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
                long_name = parts[0].strip().strip('"')
                arg_type = parts[1].strip()
                shortform = parts[3].strip()

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
                if long_name == '0':
                    continue

                takes_arg = (arg_type == 'required_argument')
                options[long_name] = OptionInfo(takes_arg, shortform)
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

def extract_switch_optarg(text: str) -> dict:
    """Compile optarg usage within switch(c) block.

    Looks within the switch block handling options,
    starting with a line with `switch(c)` (plus optional spaces)
    and ending with the outermost closing curly brace `}`.
    (i.e. nested curly braces are handled correctly)

    Ignores lines with comments. Handles fallthroughs,
    breaking cases on either a return or a `break;`.

    Looks for whether each case uses `optarg` or not.

    Parameters
    ----------
    text : str
        The text of the file to search for long_options.

    Returns
    -------
    Dict[str, bool]
        A dictionary mapping shortform names to
        whether they use optarg (i.e. take an argument).
    """

    optarg_usage = {}
    current_cases = []
    has_optarg = False
    inside_switch = False

    def set_optarg_usage(cases: list[str], usage: Optional[bool]):
        """Set optarg usage for a list of cases."""
        for case in cases:
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
        if case_match:
            # Flush old cases to avoid optarg
            # being present after fallthrough
            if has_optarg:
                set_optarg_usage(current_cases, True)
            
            stripped = stripped.split('//')[0].strip()  # Remove comments
            case_value = case_match.group(1)
            if case_value.startswith("'") and case_value.endswith("'"):
                case_key = case_value.strip("'")
            else:
                case_key = case_value  # All-caps variable
            current_cases.append(case_key)
            continue

        # If it's not a new case, it belongs to the current block
        has_optarg = has_optarg or 'optarg' in stripped

        # On crash/deprecation, flush the current case group
        if (('return' in stripped or 'exit(' in stripped 
             or 'deprecated' in stripped) and curly_brace_nesting == 0):
            # None is a marker for a crash
            set_optarg_usage(current_cases, None)

        # On break, flush the current case group
        if stripped == 'break;' and curly_brace_nesting == 0:
            set_optarg_usage(current_cases, has_optarg)

    # Handle trailing block without break (not common but valid)
    if current_cases:
        set_optarg_usage(current_cases, has_optarg)

    return optarg_usage

def check_file(filepath: str):
    """
    Check consistency of options in a single _main.cpp file.
    Returns list of problematic long options.
    """
    with open(filepath) as f:
        text = f.read()

    help_opts = extract_help_options(text)
    long_opts = extract_long_options(text)
    getopt_opts = extract_getopt_string(text)
    switch_opts = extract_switch_optarg(text)

    all_longform = set(help_opts) | set(long_opts)

    # Confirm that the help options are present
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

    # Check longform options between the four sources
    for longform in all_longform:
        cur_help = help_opts.get(longform, OptionInfo())

        # Check 1: all options must appear in long_options[] with shortform
        if longform not in long_opts:
            print(f"{filepath}: --{longform} is not in long_options[]")
            continue
        elif long_opts[longform].shortform is None:
            print(f"{filepath}: --{longform} is in long_options[] but has no "
                  "shortform")
            continue

        cur_long = long_opts[longform]
        should_str = 'should' if cur_long.takes_argument else "shouldn't"

        # Check 1: long_options[] don't use raw numbers
        if isinstance(cur_long.shortform, int):
            print(f"{filepath}: --{longform} has an int ({cur_long.shortform}) "
                  "in long_options[]; use a char or ALL_CAPS variable instead")
            continue

        # Check 2: help vs long_options[]
        if not cur_help.is_unset() and cur_help != cur_long:
            print(f"{filepath}: --{longform} has mismatch between helptext "
                  f"{cur_help} and long_options[] {cur_long}")
            continue

        # Check 3: long_options[] vs getopt string
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

        # Check 4: long_options[] vs switch
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
    
    # getopt is allowed to have ? as a help alias
    # without that being in long_options[]
    if '?' in getopt_opts:
        getopt_opts.pop('?')
    if getopt_opts:
        print(f"{filepath}: getopt string has options not in long_options[]: "
              f"{', '.join(getopt_opts)}")
    
    if '?' in switch_opts:
        switch_opts.pop('?')
    if switch_opts:
        print(f"{filepath}: switch block has options not in long_options[]: "
              f"{', '.join(switch_opts)}")

if __name__ == "__main__":
    subcommand_dir = 'src/subcommand'
    # Some subcommand files have a different format
    ignore_files = SKIP_FILES
    for fname in os.listdir(subcommand_dir):
        if not fname.endswith('_main.cpp') or fname in ignore_files:
            continue
        fpath = os.path.join(subcommand_dir, fname)
        problems = check_file(fpath)

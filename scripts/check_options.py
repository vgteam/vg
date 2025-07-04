"""Check that command line options are correctly registered.

Reads the options within the helptext (help_<command>() function),
the long_options[] array, the getopt_long() string, and the switch(c)
block, and checks that they are consistent with each other.

Mostly written by ChatGPT with small tweaks by @faithokamoto
"""
import os
import re

def extract_help_options(text: str) -> dict:
    """
    Extract options from the help_<command>() function.
    Returns a dict mapping long option names to tuples: (short_opt or None, takes_argument: bool)
    """
    help_opts = {}
    help_pattern = re.compile(r'<< "\s+(?:^|[^\S\r\n])(-\w),?\s+(--[a-zA-Z0-9\-]+)(?:\s[A-Z_]+)?')
    long_only_pattern = re.compile(r'<< "\s+(?:^|[^\S\r\n])(?!-)(--[a-zA-Z0-9\-]+)(?:\s[A-Z_]+)?')
    
    inside_help = False
    for line in text.splitlines():
        if re.search(r'void\s+help_\w+\s*\(', line):
            inside_help = True
        elif inside_help and '}' in line:
            inside_help = False
        if not inside_help or line.strip().startswith('//'):
            continue

        match = help_pattern.search(line)
        if match:
            short_opt = match.group(1)[1]
            long_opt = match.group(2)[2:]
            takes_arg = bool(re.search(rf'{re.escape(match.group(2))}\s[A-Z_]+', line))
            help_opts[long_opt] = (short_opt, takes_arg)
            continue

        match = long_only_pattern.search(line)
        if match:
            long_opt = match.group(1)[2:]
            takes_arg = bool(re.search(rf'{re.escape(match.group(1))}\s[A-Z_]+', line))
            help_opts[long_opt] = (None, takes_arg)

    return help_opts

def extract_long_options(text: str) -> dict:
    """
    Extract options from long_options[] definition.
    Returns a dict mapping long option names to tuples: (short_opt (char or ALL_CAPS) or None, takes_argument: bool)
    """
    options = {}
    inside = False
    for line in text.splitlines():
        if ('struct option long_options' in line
            or 'std::vector<struct option> long_options' in line):
            inside = True
        elif inside and '};' in line:
            inside = False
        if not inside:
            continue
        line = line.strip()
        if line.startswith('{') and not line.startswith('{0'):
            parts = line.split('//')[0].strip('{} \t,\n').split(',')
            if len(parts) >= 4:
                long_name = parts[0].strip().strip('"')
                arg_type = parts[1].strip()
                shortform = parts[3].strip().strip("'")
                if long_name == '0':
                    continue
                takes_arg = (arg_type == 'required_argument')
                short_opt = None if shortform == '0' else shortform
                options[long_name] = (short_opt, takes_arg)
    return options

def extract_getopt_string(text: str) -> set:
    """
    Extract short options from the getopt_long() string (e.g. "hx:rt:")
    Returns a set of short options, with ':' indicating required arguments
    """
    match = re.search(r'getopt_long\s*\([^,]+,[^,]+,\s*"([^"]+)"', text)
    if match:
        opts = match.group(1)

    # Second try to catch vg giraffe/augment weirdness
    if not match:
        match = re.search(r'short_options = "[^"]+"', text)
    if match:
        opts = match.group(0)

    # Give up
    if not match:
        return set()
    
    result = set()
    i = 0
    while i < len(opts):
        if i+1 < len(opts) and opts[i+1] == ':':
            result.add(f'{opts[i]}:')
            i += 2
        else:
            result.add(opts[i])
            i += 1
    return result

def extract_switch_optarg(text: str) -> dict:
    """
    Returns a dict of short_opt -> uses_optarg (True/False), accounting for fallthroughs.
    """
    switch_block = re.search(r'switch\s*\(\s*c\s*\)\s*{(.+?)}', text, re.DOTALL)
    if not switch_block:
        return {}

    body = switch_block.group(1)
    lines = body.splitlines()

    optarg_usage = {}
    current_cases = []
    block_lines = []

    def process_block(cases, block):
        block_text = "\n".join(block)
        uses_optarg = 'optarg' in block_text
        for case in cases:
            optarg_usage[case] = uses_optarg

    for line in lines:
        stripped = line.strip()

        # Detect new case
        case_match = re.match(r'case\s+([A-Za-z0-9_\'"]+)\s*:', stripped)
        if case_match:
            stripped = stripped.split('//')[0].strip()  # Remove comments
            case_value = case_match.group(1)
            if case_value.startswith("'") and case_value.endswith("'"):
                case_key = case_value.strip("'")
            else:
                case_key = case_value  # All-caps variable
            current_cases.append(case_key)
            continue

        # If it's not a new case, it belongs to the current block
        block_lines.append(stripped)

        # On break/return, flush the current case group
        if stripped == 'break;' or stripped.startswith('return'):
            process_block(current_cases, block_lines)
            current_cases = []
            block_lines = []

    # Handle trailing block without break (not common but valid)
    if current_cases or block_lines:
        process_block(current_cases, block_lines)

    return optarg_usage

def check_file(filepath: str) -> list:
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

    all_longopts = set(help_opts) | set(long_opts)

    for long_opt in all_longopts:
        help_short, help_arg = help_opts.get(long_opt, (None, None))
        long_short, long_arg = long_opts.get(long_opt, (None, None))

        # Check 1: long_options[] don't use raw numbers
        if isinstance(long_short, int):
            print(f"{filepath}: --{long_opt} has a raw number ({long_short}) "
                  "in long_options[]; use a character or ALL_CAPS variable instead")
            continue

        # Check 2: help vs long_options[]
        if ((help_short, help_arg) != (None, None) 
            and (help_arg != long_arg if help_short is None else
                 (help_short, help_arg) != (long_short, long_arg))):
            print(f"{filepath}: --{long_opt} has mismatch between helptext "
                  f"(short: -{help_short}, arg: {help_arg}) and long_options[] "
                  f"(-{long_short}, arg: {long_arg})")
            continue

        # Check 3: long_options[] vs getopt string
        if long_short and len(long_short) == 1:
            short_entry = f'{long_short}:' if long_arg else long_short
            if short_entry not in getopt_opts:
                print(f"{filepath}: --{long_opt}'s -{long_short} should be in "
                      f"getopt string as '{short_entry}'")
                continue

            if not long_arg and f'{long_short}:' in getopt_opts:
                print(f"{filepath}: --{long_opt}'s -{long_short} shouldn't have "
                      f"a : after it in getopt string")
                continue

        # Check 4: long_options[] vs switch(optarg)
        if long_short:
            used = switch_opts.get(long_short)

            if used is None:
                print(f"{filepath}: --{long_opt}'s -{long_short} is "
                      "not used in the switch block")
                continue

            if used != long_arg:
                print(f"{filepath}: --{long_opt}'s -{long_short} "
                      f"should {'not ' if not long_arg else ''}use optarg")
                continue
    
if __name__ == "__main__":
    subcommand_dir = 'src/subcommand'
    for fname in os.listdir(subcommand_dir):
        if not fname.endswith('_main.cpp'):
            continue
        fpath = os.path.join(subcommand_dir, fname)
        problems = check_file(fpath)

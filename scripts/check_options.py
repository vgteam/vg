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
    help_pattern = re.compile(r'(?:^|[^\S\r\n])(-\w),?\s+(--[a-zA-Z0-9\-]+)(?:\s[A-Z_]+)?')
    long_only_pattern = re.compile(r'(?:^|[^\S\r\n])(?!-)(--[a-zA-Z0-9\-]+)(?:\s[A-Z_]+)?')
    
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
        if 'struct option long_options' in line:
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
    if not match:
        return set()
    opts = match.group(1)
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
    Returns a dict of short_opt -> uses_optarg (True/False)
    """
    switch_block = re.search(r'switch\s*\(\s*c\s*\)\s*{(.+?)}', text, re.DOTALL)
    if not switch_block:
        return {}
    body = switch_block.group(1)
    cases = re.split(r'case\s+', body)[1:]
    optarg_usage = {}
    for case in cases:
        lines = case.strip().splitlines()
        if not lines:
            continue
        key = lines[0].strip().strip(':')
        uses_optarg = any('optarg' in l for l in lines)
        optarg_usage[key] = uses_optarg
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

    problems = []

    for long_opt in all_longopts:
        help_short, help_arg = help_opts.get(long_opt, (None, None))
        long_short, long_arg = long_opts.get(long_opt, (None, None))

        # Check 1: help vs long_options[]
        if ((help_short, help_arg) != (None, None) 
            and (help_short, help_arg) != (long_short, long_arg)):
            problems.append(long_opt)
            continue

        # Check 2: long_options[] vs getopt string
        if long_short and long_short.islower():
            short_entry = f'{long_short}:' if long_arg else long_short
            if short_entry not in getopt_opts:
                problems.append(long_opt)
                continue

        # Check 3: long_options[] vs switch(optarg)
        if long_short:
            used = switch_opts.get(long_short)
            if used is not None and used != long_arg:
                problems.append(long_opt)

    return problems
    
if __name__ == "__main__":
    subcommand_dir = 'src/subcommand'
    for fname in os.listdir(subcommand_dir):
        if not fname.endswith('_main.cpp'):
            continue
        fpath = os.path.join(subcommand_dir, fname)
        problems = check_file(fpath)
        if problems:
            print(f"{fname}: {', '.join(sorted(problems))}")

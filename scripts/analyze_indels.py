#!/usr/bin/python3

"""
analyze_indels.py: take GAM JSON as input and produce tables of non-softclip
indel positions in the read and lengths
"""

import argparse
import os
import sys
import time
import subprocess
import collections
import json

    
def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("--input", type=argparse.FileType('r'), default=sys.stdin,
                        help="JSON GAM to process")
    parser.add_argument("--positions", required=True, type=argparse.FileType('w'),
                        help="TSV of positions and counts to write")
    parser.add_argument("--lengths", required=True, type=argparse.FileType('w'),
                        help="TSV of lengths and counts to write")
    parser.add_argument("--mapqs", required=True, type=argparse.FileType('w'),
                        help="TSV of indel/gapless categories, MAPQs and counts to write")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # These will track stats over all the indels.
    
    # Where does each indel start in read space?
    indel_starts = collections.Counter()
    # How long is each indel?
    indel_lengths = collections.Counter()
    # For reads with and without indels, what are your MAPQs?
    read_mapqs = {'indel': collections.Counter(), 'gapless': collections.Counter()}
    
    for line in options.input:
        # Load each line
        alignment = json.loads(line)
        
        # Parse the alignment using these counters
        read_offset = 0
        ref_offset = 0
        
        # Determine if the read is 'indel' or 'gapless'
        status = 'gapless'
        
        mappings = alignment.get("path", {}).get("mapping", []) 
        
        for i, mapping in enumerate(mappings):
            edits = mapping.get("edit", [])
            for j, edit in enumerate(edits):
                # Parse each edit
                from_length = edit.get("from_length", 0)
                to_length = edit.get("to_length", 0)
                
                # If we're the first or last edit, we're never an indel; we're a softclip instead
                is_border = (i == 0 and j == 0) or (i == len(mappings) - 1 and j == len(edits) - 1)
                
                # If the lengths differ and we aren't on the end of the read we're an indel.
                is_indel = from_length != to_length and not is_border
                
                if is_indel:
                    # Record this indel
                    indel_starts[read_offset] += 1
                    indel_lengths[to_length - from_length] += 1
                    # Mark that the read has an indel
                    status = 'indel'
                
                # Advance the position counters
                read_offset += to_length
                ref_offset += from_length
                
        # Record the read's MAPQ
        read_mapqs[status][alignment.get("mapping_quality", 0)] += 1
                
    # Now dump tables for plotting
    for value, count in indel_starts.items():
        options.positions.write("{}\t{}\n".format(value, count))
        
    for value, count in indel_lengths.items():
        options.lengths.write("{}\t{}\n".format(value, count))
      
    for status, counts in read_mapqs.items():
        for value, count in counts.items():
            options.mapqs.write("{}\t{}\t{}\n".format(status, value, count))


def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        

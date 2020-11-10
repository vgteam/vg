#!/usr/bin/python3
# filter_unpaired_mappings.py: Filter alignments so that every pair of lines represents a valid read pair
"""

Filter out unpaired mappings so that the results can be used as input for interleaved, paired
alignment.  Ensure that for each even value of i, then lines i and i+1 contain mappings for
the first and second fragment for a given read.  Input must be sorted by name

Examples:

  # Note: GAM sorting with jq takes loads of memory
  vg view -a unsorted.gam | jq -s -c 'sort_by(.name)[]' | filter_unpaired_mappings.py - json-gam | vg view -JaG - > interleaved.gam

  samtools sort -n unsorted.bam | filter_unpaired_mappings.py - sam | samtools view - -O BAM > interleaved.bam

"""
import argparse
import os
import sys
import time
import subprocess
import collections
import json
try:
    import pysam
    have_pysam = True
except:
    have_pysam = False
    
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
    
    parser.add_argument("input",
                        help="JSON GAM, BAM or SAM file (sorted by name) to filter (use - for stdin).")
    parser.add_argument("type", choices = ['json-gam', 'sam', 'bam'],
                        help="Input format")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def filter_sam(options):
    """
    use pysam to iterate over sam/bam input
    """
    
    in_file = pysam.AlignmentFile(options.input, 'r' if options.type == 'sam' else 'rb')
    out_file = pysam.AlignmentFile('-', 'w' if options.type == 'sam' else 'wb', template=in_file)

    prev_name = None
    first_record = None

    num_written = 0
    num_read = 0
    for record in in_file:
        num_read += 1
        if record.is_secondary:
            continue
        if record.is_read1:
            first_record = record
        elif record.is_read2 and first_record and record.query_name == first_record.query_name:
            out_file.write(first_record)
            out_file.write(record)
            num_written += 2

        # samtools sorts names differently then python 
        #if prev_name and record.query_name < prev_name:
        #    raise RuntimeError('Records not sorted by name: {} {}'.format(prev_name, record.query_name))
            
        prev_name = record.query_name

    in_file.close()
    out_file.close()

    sys.stderr.write('Filtered {} / {} lines\n'.format(num_read - num_written, num_read))

def filter_json_gam(options):
    """
    process gam json made with vg view -a my.gam | jq -s -c 'sort_by(.name)[]'
    """

    in_file = sys.stdin if options.input == '-' else open(options.input, 'r')

    prev_name = None
    first_record = None

    num_written = 0
    num_read = 0    
    for line in in_file:
        num_read += 1
        gam = json.loads(line)
        if 'is_secondary' in gam and gam['is_secondary']:
            continue
        if 'fragment_next' in gam:
            first_record = gam
        elif 'fragment_prev' in gam and first_record and gam['fragment_prev']['name'] == first_record['name']:
            assert first_record['fragment_next']['name'] == gam['name']
            print(json.dumps(first_record))
            print(json.dumps(gam))
            num_written += 2
        elif 'fragment_next' not in gam and 'fragment_prev' not in gam:
            raise RuntimeError('fragment_prev/next not set for record: {}'.format(json.dumps(gam)))

        if prev_name and gam['name'] < prev_name:
            raise RuntimeError('Records not sorted by name: {} {}\n'.format(prev_name, gam['name']) + 
                               'Use vg view -a {} | jq -s -c \'sort_by(.name)[]\''.format(options.input))

            
        prev_name = gam['name']

    in_file.close()

    sys.stderr.write('Filtered {} / {} lines\n'.format(num_read - num_written, num_read))
        

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object

    if options.type in ['sam', 'bam']:
        if not have_pysam:
            raise RuntimeError('pysam needed from SAM/BAM support')
        return filter_sam(options)
    elif options.type == 'json-gam':
        return filter_json_gam(options)

def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        

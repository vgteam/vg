#!/usr/bin/env python3



"""
giraffe-facts.py: process a GAM file from the new minimizer-based mapper (vg
gaffe AKA giraffe) and report runtime statistics by pipeline stage.
"""

import argparse
import os
import sys
import time
import subprocess
import collections
import json

# We depend on our local histogram.py
import histogram

# Define a constant list of all the stages, in order.
STAGES = ['minimizer', 'seed', 'cluster', 'extend', 'align', 'winner']
    
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
                        help="line-oriented JSON GAM to process")
    parser.add_argument("outdir",
                        help="directory to place output in")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def accumulate_read_times(reads, stages=STAGES):
    """
    Given an iterator that emits parsed JSON reads dicts, yield one list of
    cumulative runtimes, in stages order, for each read. Each list also ends
    with the total runtime, after all the stages.
    """
    
    for read in reads:
        # For every read
        
        # We make a list of cumulative times
        cumulative_times = []
        # And we track the time we are summing up
        time_counter = 0.0
        
        for stage in stages:
            # For each stage in order
            
            # Grab its runtime from the read
            stage_time = read.get('annotation', {}).get('stage_' + stage + '_seconds', 0.0)
            
            # Sum it in
            time_counter += stage_time
            
            # Record the cumulative sum through this stage
            cumulative_times.append(time_counter)
            
        # Now do the total time
        cumulative_times.append(read.get('annotation', {}).get('map_seconds', 0.0))
        
        yield cumulative_times
            
def read_line_oriented_json(lines):
    """
    For each line in the given stream, yield it as a parsed JSON object.
    """
    
    for line in lines:
        yield json.loads(line)


def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    print("Thank you for subscribing to Giraffe Facts")
    
    # Make the output directory if it doesn't exist
    os.makedirs(options.outdir, exist_ok=True)
    
    # Open a TSV file to draw a histogram from
    tsv_path = os.path.join(options.outdir, 'times.tsv')
    tsv = open(tsv_path, 'w')
    
    # Make a place to total up times by stage and cumulative time
    time_totals = [0.0 for _ in STAGES] + [0.0]
    
    # Count all the reads
    read_count = 0
    
    for times in accumulate_read_times(read_line_oriented_json(options.input)):
        # For the cumulative runtimes for each read, ending with total time
        
        for stage, cumulative_time in zip(STAGES, times):
            # Dump cumulative times to the TSV we will plot a histogram from
            tsv.write('{}\t{}\n'.format(stage, cumulative_time))
        
        # Also include total time
        tsv.write('total\t{}\n'.format(times[-1]))
        
        # Sum up all times in accumulator
        for i, time in enumerate(times):
            time_totals[i] += time
            
        # Count the read
        read_count += 1
    
    # After processing all the reads
    
    # Close the TSV
    tsv.close()
    
    # Plot the histogram
    svg_path = os.path.join(options.outdir, 'times.svg')
    histogram.main(['histogram.py', tsv_path, '--save', svg_path,
        '--title', 'Runtime Histogram',
        '--x_label',  'Time (seconds)',
        '--line',
        '--bins', '100',
        '--log',
        '--cumulative',
        '--y_label', 'Cumulative Count',
        '--legend_overlay', 'lower right',
        '--categories', 'total'] + STAGES)
   
   
    # Now do a table
    
    # How long is the longest stage name
    stage_width = max((len(x) for x in STAGES))
    # Leave room for the header
    stage_header = "Stage"
    stage_width = max(stage_width, len(stage_header))
    # And for the "Overall" entry
    stage_overall = "Overall"
    stage_width = max(stage_width, len(stage_overall))
    
    # How about the reads per second column
    rps_header = "Reads per Second"
    rps_header2 = "(cumulative)"
    rps_width = max(len(rps_header), len(rps_header2))
    
    print("┌{}┬{}┐".format('─' * stage_width, '─' * rps_width))
    print("│{}│{}│".format(stage_header.rjust(stage_width), rps_header.rjust(rps_width)))
    print("│{}│{}│".format(' ' * stage_width, rps_header2.rjust(rps_width)))
    print("├{}┼{}┤".format('─' * stage_width, '─' * rps_width))
    
    for stage, total_cumulative_time in zip(STAGES, time_totals):
        # Compute cumulative reads per second values
        reads_per_second = read_count / total_cumulative_time if total_cumulative_time != 0 else float('NaN')
        
        # Justify right with spaces in a field of the correct width.
        print(("│{}│{:>" + str(rps_width) + ".2f}│").format(stage.rjust(stage_width), reads_per_second))
    
    # And do overall reads per second, after a divider
    print("├{}┼{}┤".format('─' * stage_width, '─' * rps_width))
    reads_per_second_overall = read_count / time_totals[-1] if time_totals[-1] != 0 else float('NaN')
    print(("│{}│{:>" + str(rps_width) + ".2f}│").format(stage_overall.rjust(stage_width), reads_per_second_overall))
    
    print("└{}┴{}┘".format('─' * stage_width, '─' * len(rps_header)))

def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        


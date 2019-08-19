#!/usr/bin/env python3



"""
giraffe-facts.py: process a GAM file from the new minimizer-based mapper (vg
gaffe AKA giraffe) and report runtime statistics by filter.
"""

import argparse
import os
import sys
import time
import subprocess
import collections
import itertools
import json
import random

# We depend on our local histogram.py
import histogram

FACTS = ["Giraffes are the tallest living terrestrial animal.",
         "There are nine subspecies of giraffe, each occupying separate regions of Africa. Some researchers consider some subspecies to be separate species entirely.",
         "Giraffes' horn-like structures are called 'ossicones'. They consist mostly of ossified cartilage.",
         "Male giraffes compete for dominance in fights in which they strike each other with their necks.",
         "There are more than 1600 giraffes in captivity worldwide.",
         "The name 'giraffe' has roots in the Arabic 'zarafah', meaning 'fast-walker'.",
         "Before the name 'giraffe' came into standard use, giraffes were commonly called 'camelopards'.",
         "There are 10 known extinct species of giraffe.",
         "The closest living relative to the giraffe is the okapi, an endangered hoofed mammal from the Congo.",
         "Full grown giraffes are usually between 14 and 18 feet tall.",
         "The tallest recorded giraffe was 19.3 feet tall.",
         "Adult male giraffes weigh an average of 2628 lbs., whereas females weight 1825 lbs.",
         "Giraffes have the ability to close their nostrils to protect against sandstorms and ants.",
         "Giraffes have 18-inch-long prehensile tongues, which they use for grasping foliage and for grooming.",
         "Male giraffes' spots grow darker as they age.",
         "Under their fur coat, giraffes have grey skin.",
         "Female giraffes have hair on their ossicones, whereas males' ossicones are bald.",
         "Giraffes use the weight of their head to maintain their balance when they gallop.",
         "Giraffes can run at 37 miles per hour for short distances, and 31 miles per hour for several miles.",
         "Giraffes sleep for about half an hour a day.",
         "Giraffes have the same number of vertebrae as most mammals. The length of their neck comes from longer vertebrae (over 10 inches each).",
         "Giraffes' neck is fairly short at birth, probably to make birthing easier for mothers.",
         "A giraffe's heart can weigh more than 25 lbs.",
         "Giraffes have structures like check valves in their necks' veins to prevent blood from rushing to their head when they bend down to drink.",
         "Giraffes have a four-chambered stomach similar to cattle.",
         "An adult girafffe can eat 75 lbs. of foliage per day.",
         "While generally herbivorous, giraffes have been observed eating meat and bone from carcasses.",
         "The giraffe's gestation period is 14 months.",
         "Newborn giraffes are about 6 feet tall.",
         "Giraffes are lions' most common prey.",
         "Most of giraffes' mounting behavior is between two males, often after a fight for dominance.",
         "Giraffes allow red-billed oxpeckers (a bird species) to perch on them to feed on ticks.",
         "Egyptian heiroglyphs use the giraffe as a character, pronounced 'sr'.",
         "Designers of suits for fighter pilots studied giraffe skin, since figher pilots are also at risk of passing out when blood rushes to the legs.",
         "The Humr people of Sudan use giraffe liver to create a putatively hallucinogenic drink called 'umm nyolokh'. The drink's psychoactive properties may come from the giraffe's diet of acacia plants.",
         "The giraffe is the national animal of Tanzania.",
         "There are around 100,000 giraffes in the wild as of 2016.",
         "Giraffes only need to drink every few days. Most of their water comes from the vegetation they eat.",
         "Giraffes give birth standing up, so newborn giraffes fall over 5 feet upon being born.",
         "Giraffes usually sleep standing upright.",
         "Male giraffes detect oestrus in females by tasting their urine.",
         "June 21 is World Giraffe Day.",
         "Toys R' Us has used Geoffrey the Giraffe as its mascot since 1965, although earlier advertisements in the 1950's used another giraffe: Dr. G. Raffe.",
         "Giraffe hooves are 1 foot in diameter.",
         "About 50% of giraffe calves die in their first year, mostly due to predation.",
         "Kirahvi sanoo öri öri öri öri öri öri.",
         "The giraffe's average walking speed is 10 miles per hour.",
         "The giraffe's tongue is colored dark blue.",
         "Some of giraffes' vocalizations are too low to be heard by human ears.",
         "Giraffes have never been observed swimming.",
         "Mozambique requires power lines to be 39 feet high so giraffes can safely pass underneath."]

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
    
def make_stats(read):
    """
    Given a read dict parsed from JSON, compute a stats OrderedDict for the read.
    
    Run on an empty dict, makes a zero-value stats dict.
    
    A stats dict maps from filter name to a Counter of filter stats.
    The filter stats include:
        
        - 'passed_count_total' which is the count of results passing the filter.
        - 'failed_count_total' which is the count of results failing the filter.
        - 'passed_count_correct' which is the count of correct results passing the filter.
        - 'failed_count_correct' which is the count of correct results failing the filter.
        
    Additionally, each of these '_count_' stats has a '_size_' version,
    describing the total size of all items meeting the specified criteria (as
    opposed to the number of items).
        
    Filters appear in the OrderedDict in an order corresponding to their filter
    number in the GAM.
    """
    
    # This is the annotation dict from the read
    annot = read.get('annotation', {})
    
    # This will map from filter number int to filter name
    filters_by_index = {}
    
    # This will map from filter name to Counter of filter stats
    filter_stats = collections.defaultdict(collections.Counter)
    
    for annot_name in annot.keys():
        # For each annotation
        if not annot_name.startswith('filter_'):
            # Skip the ones that aren't about filters
            continue
            
        # Names look like 'filter_2_cluster-score-threshold_cluster_passed_size_correct'
            
        # Break into components on underscores
        (_, filter_num, filter_name, filter_stage, filter_status, filter_accounting, filter_metric) = annot_name.split('_')
        
        # Collect integers
        filter_num = int(filter_num)
        filter_stat_value = annot[annot_name]
        
        # Record the filter being at this index if not known already
        filters_by_index[filter_num] = filter_name
        
        # Record the stat value
        filter_stats[filter_name]['{}_{}_{}'.format(filter_status, filter_accounting, filter_metric)] = filter_stat_value
        
    # Now put them all in this OrderedDict in order
    ordered_stats = collections.OrderedDict()
    for filter_index in sorted(filters_by_index.keys()):
        filter_name = filters_by_index[filter_index]
        ordered_stats[filter_name] = filter_stats[filter_name]
        
    return ordered_stats

def add_in_stats(destination, addend):
    """
    Add the addend stats dict into the destinatin stats dict.
    Implements += for stats dicts.
    """
    
    # Internally we're just recursive += on dicts.
    for k, v in addend.items():
        if isinstance(v, dict):
            # Recurse into dict
            add_in_stats(destination[k], v)
        else:
            # Use real += and hope it works
            destination[k] += v

def read_line_oriented_json(lines):
    """
    For each line in the given stream, yield it as a parsed JSON object.
    """
    
    for line in lines:
        yield json.loads(line)
        
class Table(object):
    """
    Format a table of output nicely in fixed-width text.
    """
    
    # Interface
    
    def __init__(self, widths, out=sys.stdout):
        """
        Start a table with the given column widths (a list of integers) in
        characters, printing to the given stream.
        """
        
        # Remember the base widths
        self.widths = widths
        
        # Remember the out stream
        self.out = out
        
        # Remember the previous actual column widths used, if any.
        # None if no wor has been produced.
        self.last_widths = None
        
        # Remember if we need a dividing line
        self.need_line = False
        
    def line(self):
        """
        Say to divide the previous row from the next row.
        """
        
        self.need_line = True
        
    def row(self, values, justify='l', merge=None, line_top=False, line_bottom=False):
        """
        Given a list of values, one per column, for up to the number of columns
        in the table, draw a table row.
        
        Justify can be 'l', 'r', 'c' or a list/string of those per value.
        
        If merge is given, it must be a list of the number of cells to merge
        horizontally for each value.
        
        Different merge values without a line_top separator will look bad.
        If line_top is set, divide from the previous row.
        If line_bottom is set, divide from the next row.
        """
        
        # Compute new merged widths
        merged_widths = self.compute_merges(merge)
        
        # Start or continue the table
        if self.last_widths is None:
            # Start the table
            self.start(merged_widths)
        elif self.need_line or line_top:
            # Divide from the previous row.
            self.sep(self.last_widths, merged_widths)
            
        # Print the actual row
        self.cells(values, justify, merged_widths)
            
        # Remember this row's widths for next time.
        self.last_widths = merged_widths
        
        # Remember if we need a line
        self.need_line = line_bottom
        
    def close(self):
        """
        Close off the table at the bottom.
        """
        
        if self.last_widths is None:
            self.last_widths = self.widths
        
        self.end(self.last_widths)
        
        self.last_widths = None
        
    def inner_width(self):
        """
        Get the total width of the table across all columns, between the outer edges.
        """
        
        return sum(self.widths) + len(self.widths) - 1
        
    # Internal methods
        
    def box(self, part):
        """
        Return the box-drawing character to draw the given part of a box.
        Parts are {(t)op, (m)iddle, (b)ottom} crossed with {(l)eft, (m)iddle,
        (r)ight} as two-character strings, plus (v)ertical and (h)orizontal as one-character strings.
        """
        
        skin = {
            'tl': '┌',
            'tm': '┬',
            'tr': '┐',
            'bl': '└',
            'bm': '┴',
            'br': '┘',
            'ml': '├',
            'mm': '┼',
            'mr': '┤',
            'v': '│',
            'h': '─'
        }
        
        return skin[part]
        
    def horizontal(self, left, junction, right, column, widths=None):
        """
        Print a line across (either top, middle, or bottom).
        
        Takes the leftmost, between-column, rightmost, and in-column characters
        as box() character ID strings.
        
        Can use a specified widths list, usually self.widths.
        """
        
        if widths is None:
            widths = self.widths
        
        # Start edge
        self.out.write(self.box(left))
        
        for i, width in enumerate(widths):
            # For each column
            # Do its top line
            self.out.write(self.box(column) * width)
            if i + 1 != len(widths):
                # Do the separator
                self.out.write(self.box(junction))
                
        # End edge
        self.out.write(self.box(right))
        
        self.out.write('\n')
        
    def start(self, widths_after):
        """
        Print an opening line at the top of the table.
        Needs to know the widths of the cells on the next table line.
        """
        
        self.horizontal('tl', 'tm', 'tr', 'h', widths_after)
        
    def end(self, widths_before):
        """
        Print a closing line at the bottom of the table.
        Needs to know the widths of the cells on the previous table line.
        """
        
        self.horizontal('bl', 'bm', 'br', 'h', widths_before)
        
    def sep(self, widths_before, widths_after):
        """
        Print a middle separator line across the table.
        Needs to know the widths of the cells on the previous and next table lines.
        Both sets of widths must describe a table of the same total width.
        """
        
        # Start edge
        self.out.write(self.box('ml'))
        
        # Compute total width (cells and separators), not counting final border
        total_width = sum(widths_before) + len(widths_before) - 1
        
        # Track what cell we are in on top
        before_cursor = 0
        # And what column its trailing border is at
        before_border = widths_before[before_cursor]
        # Track what cell we are in on the bottom
        after_cursor = 0
        # And what column its trailing border is at
        after_border = widths_after[after_cursor]
        # Track what column of internal table width we are in.
        col = 0
        
        while col < total_width:
            if col == before_border:
                if col == after_border:
                    # Junction on both sides
                    char = self.box('mm')
                    
                    # Advance after
                    after_cursor += 1
                    after_border += widths_after[after_cursor] + 1
                else:
                    # Junction on top only
                    char = self.box('bm')
                    
                # Advance before
                before_cursor += 1
                before_border += widths_before[before_cursor] + 1
            elif col == after_border:
                # Junction on bottom only
                char = self.box('tm')
                
                # Advance after
                after_cursor += 1
                after_border += widths_after[after_cursor] + 1
            else:
                # No junction
                char = self.box('h')
                
            # Print the character
            self.out.write(char)
            
            # Go to the next column
            col += 1
        
        
        # End edge
        self.out.write(self.box('mr'))
        
        self.out.write('\n')
        
    def compute_merges(self, merges=None):
        """
        Given a list of cell counts to merge horizontally, compute new widths from self.widths.
        
        If merges is None, use self.widths.
        """
        
        widths = self.widths
        
        if merges is not None:
            new_widths = []
            width_cursor = 0
            for merge in merges:
                # Compute a new column by merging the given number of old columns.
                merged_width = 0
                for i in range(merge):
                    # Take the widths of all cells
                    merged_width += widths[width_cursor]
                    width_cursor += 1
                # Take the separating columns between cells
                merged_width += merge - 1
                new_widths.append(merged_width)
            while width_cursor < len(widths):
                # Copy any unmerged columns
                new_widths.append(widths[i])
                
            widths = new_widths
            
        return widths
        
    def cells(self, values, justify, widths):
        """
        Given a list of values, one per column, for up to the number of columns
        in the table, draw a table row.
        
        Justify can be 'l', 'r', 'c', or a list/string of those per value.
        
        Column count/widths must be passed.
        """
       
        # Start the row
        self.out.write(self.box('v'))
        
        for i, (value, width) in enumerate(itertools.zip_longest(values, widths)):
            # For each item and its column and width...
            if width is None:
                # Too many items
                raise RuntimeError("Ran out of table width values ({}) for {} columns".format(
                    len(widths), len(values)))
              
            # Compute the item string
            item_string = str(value) if value is not None else ''
              
            # Decide on justification for this item
            if justify == 'l':
                item_just = 'l'
            elif justify == 'r':
                item_just = 'r'
            if justify == 'c':
                item_just = 'c'
            elif i < len(justify):
                item_just = justify[i]
            else:
                item_just = 'l'
                
            # Actually justify it in a field of the necessary width
            if item_just == 'l':
                justified_item = item_string.ljust(width)
            elif item_just == 'r':
                justified_item = item_string.rjust(width)
            elif item_just == 'c':
                justified_item = item_string.center(width)
            else:
                raise RuntimeError('Invalid justification: {}'.format(item_just))
                
            # Output the content
            self.out.write(justified_item)
            
            if (i + 1 != len(widths)):
                # This isn't the last item. Do a separator.
                self.out.write(self.box('v'))
                
                
        # End the row
        # TODO: Same as the separator
        self.out.write(self.box('v'))
        
        self.out.write('\n')
                
def print_table(read_count, stats_total, out=sys.stdout):
    """
    Take the read count, and the accumulated total stats dict.
    
    Print a nicely formatted table to the given stream.
    
    """
    
    if stats_total is None:
        # Handle the empty case
        assert(read_count == 0)
        out.write('No reads.\n')
        return
    
    # Now do a table
    
    # First header line for each column
    headers = []
    # Second header line for wach column
    headers2 = []
    # Column min widths from headers
    header_widths = []
    
    # How long is the longest filter name
    filter_width = max((len(x) for x in stats_total.keys()))
    # Leave room for the header
    filter_header = "Filter"
    filter_width = max(filter_width, len(filter_header))
    # And for the "Overall" entry
    filter_overall = "Overall"
    filter_width = max(filter_width, len(filter_overall))
    
    headers.append(filter_header)
    headers2.append('')
    header_widths.append(filter_width)
    
    # And the passing count columns (average)
    passing_header = "Passing"
    passing_header2 = "(Per Read)"
    passing_width = max(len(passing_header), len(passing_header2))
    
    headers.append(passing_header)
    headers2.append(passing_header2)
    header_widths.append(passing_width)
    
    # And the failing count columns (average)
    failing_header = "Failing"
    failing_header2 = "(Per Read)"
    failing_width = max(len(failing_header), len(failing_header2))
    
    headers.append(failing_header)
    headers2.append(failing_header2)
    header_widths.append(failing_width)
    
    # And the correct result lost count header
    lost_header = "Correct"
    lost_header2 = "Lost"
    # How big a number will we need to hold?
    # Look at the reads lost at all filters
    lost_reads = [stats_total[filter_name]['failed_count_correct'] for filter_name in stats_total.keys()]
    max_filter_stop = max(lost_reads)
    # How many correct reads are lost overall by filters?
    overall_lost = sum(lost_reads)
    lost_width = max(len(lost_header), len(lost_header2), len(str(max_filter_stop)), len(str(overall_lost)))
    
    headers.append(lost_header)
    headers2.append(lost_header2)
    header_widths.append(lost_width)
    
    # And the incorrect result rejected count header
    rejected_header = "Incorrect"
    rejected_header2 = "Removed"
    # How big a number will we need to hold?
    # Look at the reads rejected at all filters
    rejected_reads = [(stats_total[filter_name]['failed_count_total'] - stats_total[filter_name]['failed_count_correct'])
        for filter_name in stats_total.keys()]
    max_filter_stop = max(rejected_reads)
    # How many incorrect reads are rejected overall by filters?
    overall_rejected = sum(rejected_reads)
    rejected_width = max(len(rejected_header), len(rejected_header2), len(str(max_filter_stop)), len(str(overall_rejected)))
    
    headers.append(rejected_header)
    headers2.append(rejected_header2)
    header_widths.append(rejected_width)
    
    # Start the table
    table = Table(header_widths)
    
    table.row(["Giraffe Facts"], 'c', merge=[len(header_widths)])
    table.line()
    table.row(['Reads' + str(read_count).rjust(table.inner_width() - 5)], merge=[len(header_widths)])
    table.line()
    table.row(headers, 'c')
    table.row(headers2, 'c')
    table.line()
    
    
    for filter_name in stats_total.keys():
        # Grab average results passing this filter per read
        total_passing = stats_total[filter_name]['passed_count_total']
        average_passing = total_passing / read_count if read_count != 0 else float('NaN')
        
        # Grab average results failing this filter per read
        total_failing = stats_total[filter_name]['failed_count_total']
        average_failing = total_failing / read_count if read_count != 0 else float('NaN')
        
        # Grab reads that are lost.
        # No reads are lost at the final stage.
        lost = stats_total[filter_name]['failed_count_correct']
        
        # And reads that are (rightly) rejected
        rejected = stats_total[filter_name]['failed_count_total'] - lost
        
        row = [filter_name]
        align = 'c'
        # Add the provenance columns
        row += ['{:.2f}'.format(average_passing), '{:.2f}'.format(average_failing), lost, rejected]
        align += 'rrrr'
        
        # Output the finished row
        table.row(row, align)
        
    table.line()
        
    # Compose the overall row
    row = [filter_overall]
    align = 'c'
    # Add the provenance columns
    row += ['', overall_lost, overall_rejected]
    align += 'rr'
    
    table.row(row, align)
    
    # Close off table
    table.close()


def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    print(random.choice(FACTS), file = sys.stderr)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Make the output directory if it doesn't exist
    os.makedirs(options.outdir, exist_ok=True)
    
    # Make a place to total up all the stats
    stats_total = None 
    
    # Count all the reads
    read_count = 0
    
    
    
    for stats in (make_stats(read) for read in read_line_oriented_json(options.input)):
        # For the stats dict for each read
        
        if stats_total is None:
            stats_total = stats
        else:
            # Sum up all the stats
            add_in_stats(stats_total, stats)
        
        # Count the read
        read_count += 1
    
    # After processing all the reads
    
    # Print the table now in case plotting fails
    print_table(read_count, stats_total)
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        


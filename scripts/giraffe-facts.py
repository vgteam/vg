#!/usr/bin/python3



"""
giraffe-facts.py: process a GAM file from the new minimizer-based mapper (vg giraffe) and report runtime statistics by filter.
"""

import argparse
import os
import sys
import time
import subprocess
import collections
import io
import itertools
import json
import random
import math

# Force output to UTF-8. Eventually we can use reconfigure() if we drop 3.6
# and earlier.
# We need to do this before these streams get picked as any argument default
# values.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf8')

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
    
    parser.add_argument("input", type=str,
                        help="GAM to process")
    parser.add_argument("--vg", type=str, default="vg",
                        help="vg binary to use")
    parser.add_argument("outdir",
                        help="directory to place output in")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def sniff_params(read):
    """
    Given a read dict parsed from JSON, compute a mapping parameter dict for the read.
    
    The read will have param_XXX annotations. Turn those into a dict from XXX to value.
    
    These should be the same for every read.
    """
    
    # This is the annotation dict from the read
    annot = read.get('annotation', {})
    
    # This is the params dict to fill in
    params = {}
    
    for annot_name in annot.keys():
        if annot_name.startswith('param_'):
            # Split the param annotations on underscore
            (_, param_name) = annot_name.split('_')
            
            # Save the values under the name
            params[param_name] = annot[annot_name]
    
    return params
    
# Stats under NO_FILTER are not associated with a filter
NO_FILTER = "__none__"

def make_stats(read):
    """
    Given a read dict parsed from JSON, compute a stats OrderedDict for the read.
    
    Run on an empty dict, makes a zero-value stats dict.
    
    A stats dict maps from filter name to a Counter of filter stats.
    The filter stats include:
        
        - 'passed_count_total' which is the count of results passing the
          filter.
        - 'failed_count_total' which is the count of results failing the
          filter.
        - 'passed_count_correct' which is the count of correct results passing
          the filter.
        - 'failed_count_correct' which is the count of correct results failing
          the filter.
        
    Additionally, each of these '_count_' stats has a '_size_' version,
    describing the total size of all items meeting the specified criteria (as
    opposed to the number of items).
    
    For the 'seed' stage, correctness information is not yet available, so only
    the '_total' values will be defined. '_correct' values will be set to None
    (instead of 0).
    
    The Counter for a filter also has sub-Counters embedded in it for
    expressing distributions of filter statistic values, to assist in filter
    design.
    
        - 'statistic_distribution_correct': statistic value counts for items
          deemed correct
        - 'statistic_distribution_noncorrect': statistic value counts for items
          not deemed correct
        
    NaN values of the statistics are filtered out.
        
    Filters appear in the OrderedDict in an order corresponding to their filter
    number in the GAM.
    
    The stats dict may also have an entry for NO_FILTER, with stats not
    associated with a filter.
    """
    
    # This is the annotation dict from the read
    annot = read.get('annotation', {})
    
    # This will map from filter number int to filter name
    filters_by_index = {}
    
    # This will map from filter name to Counter of filter stats
    filter_stats = collections.defaultdict(collections.Counter)
    
    for annot_name in annot.keys():
        # For each annotation
        if annot_name.startswith('filter_'):
            # If it is an individual filter info item
            
            # Names look like 'filter_2_cluster-score-threshold_cluster_passed_size_correct'
                
            # Break into components on underscores
            (_, filter_num, filter_name, filter_stage, filter_status, filter_accounting, filter_metric) = annot_name.split('_')
            
            # Collect integers
            filter_num = int(filter_num)
            filter_stat_value = annot[annot_name]
            
            # Record the filter being at this index if not known already
            filters_by_index[filter_num] = filter_name
            
            if filter_stage == 'minimizer':
                # Wer are filtering items produced by the minimizer stage.
                # At the minimizer stage, correct and incorrect are not defined yet.
                if filter_metric == 'correct':
                    # Make sure we didn't get any counts
                    assert filter_stat_value == 0
                    # None out the correct stat so we can detect this when making the table
                    filter_stat_value = None
            
            # Record the stat value
            filter_stats[filter_name]['{}_{}_{}'.format(filter_status, filter_accounting, filter_metric)] = filter_stat_value
        
        elif annot_name.startswith('filterstats_'):
            # It is a whole collection of correct or not-necessarily-correct filter statistic distribution values, for plotting.
            
            # Break into components on underscores (correctness will be 'correct' or 'noncorrect'
            (_, filter_num, filter_name, filter_stage, filter_correctness) = annot_name.split('_')
            
            distribution = collections.Counter()
            
            for item in annot[annot_name]:
                # Parse all the statistic vlues
                item = float(item)
                
                if math.isnan(item):
                    # Discard NANs
                    continue
                    
                # Count all instances of the same value
                distribution[item] += 1
                
            # Save the statistic distribution
            filter_stats[filter_name]['statistic_distribution_{}'.format(filter_correctness)] = distribution

        elif annot_name.startswith('last_correct_stage'):
            stage = annot[annot_name]
            if stage == 'none':
                filter_stats['hard-hit-cap']['last_correct_stage'] = 1
            elif stage == 'cluster':
                filter_stats['cluster-coverage']['last_correct_stage'] = 1
            elif stage == 'extend':
                filter_stats['extension-set']['last_correct_stage'] = 1
            elif stage == 'align':
                filter_stats['max-alignments']['last_correct_stage'] = 1
        
    # Now put them all in this OrderedDict in order
    ordered_stats = collections.OrderedDict()
    for filter_index in sorted(filters_by_index.keys()):
        filter_name = filters_by_index[filter_index]
        ordered_stats[filter_name] = filter_stats[filter_name]
        
    # Add in special non-filter stats
    ordered_stats[NO_FILTER] = collections.Counter()
    for k in ['time_used']:
        if k in read:
            ordered_stats[NO_FILTER][k] = read[k]
        
    return ordered_stats

def add_in_stats(destination, addend):
    """
    Add the addend stats dict into the destination stats dict.
    Implements += for stats dicts.
    """
    
    for k, v in addend.items():
        if v is None:
            # None will replace anything and propagate through
            destination[k] = None
        elif isinstance(v, dict):
            # Recurse into dict
            if k in destination:
                add_in_stats(destination[k], v)
        else:
            # Use real += and hope it works
            destination[k] += v

def read_line_oriented_json(lines):
    """
    For each line in the given iterable of lines (such as a stream), yield it as a parsed JSON object.
    """
    
    for line in lines:
        line = line.strip()
        if len(line) > 0:
            yield json.loads(line)


def read_read_views(vg, filename):
    """
    Given a vg binary and a filename, iterate over subsets of the parsed read dicts for each read in the file.

    The subsets will have the annotation and time_used fields.
    """

    # Extract just the annotations and times of reads as JSON, with a # header
    # We don't know all the annotation field names in advance so we have to dump them all.
    filter_process = subprocess.Popen([vg, "filter", "--tsv-out", "annotation;time_used", filename], stdout=subprocess.PIPE)

    lines = iter(filter_process.stdout)
    # Drop header line
    next(lines)

    for line in lines:
        # Parse the TSV and reconstruct a view of the full read dict.
        line = line.decode('utf-8')
        line = line.strip()
        if len(line) == 0:
            continue
        parts = line.split("\t")
        assert len(parts) == 2
        read = {"annotation": json.loads(parts[0]), "time_used": float(parts[1])}

        yield read

    return_code = filter_process.wait()
    assert return_code == 0
        
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
        
    def full_row(self, left_value, right_value):
        """
        Draw a full-width row in the table with a left-justified and a
        right-justified value.
        """
        
        full_value = left_value + right_value.rjust(self.inner_width() - len(left_value))
        self.row([full_value], merge=[len(self.widths)])
        
        
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
                
def print_table(read_count, stats_total, params=None, out=sys.stdout):
    """
    Take the read count, the accumulated total stats dict, and an optional dict
    of mapping parameters corresponding to values for filters.
    
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
    
    # Find all the filters
    filters = [k for k in stats_total.keys() if k != NO_FILTER]
    
    # Compute filter row headings
    filter_headings = list(filters)
    
    if params is not None:
        # Annotate each filter with its parameter value
        annotated_headings = []
        for heading in filter_headings:
            # For each filter
            # It may be a compound thing||thing filter
            parts = heading.split('||')
            
            # We will fill this with all the relevant filter cutoff values
            filter_values = []
            for part in parts:
                if part in params:
                    filter_values.append(params[part])
                    
            if len(filter_values) == 0:
                # No parameters
                annotated_headings.append(heading)
            else:
                # Annotate with the parameters
                annotated_headings.append('{} ({})'.format(heading, ', '.join((str(x) for x in filter_values))))
                
        filter_headings = annotated_headings
    
    # How long is the longest filter name
    filter_width = max(itertools.chain((len(x) for x in filter_headings), [0]))
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
    passing_header2 = "(/Read)"
    passing_width = max(len(passing_header), len(passing_header2))
    
    headers.append(passing_header)
    headers2.append(passing_header2)
    header_widths.append(passing_width)
    
    # And the failing count columns (average)
    failing_header = "Failing"
    failing_header2 = "(/Read)"
    failing_width = max(len(failing_header), len(failing_header2))
    
    headers.append(failing_header)
    headers2.append(failing_header2)
    header_widths.append(failing_width)
    
    # And the number of correct reads lost at each stage
    lost_stage_header = "Lost"
    lost_stage_header2 = "reads"
    lost_stage_reads = [x for x in (stats_total[filter_name].get('last_correct_stage', 0) for filter_name in filters) if x is not None]
    max_stage = max(itertools.chain(lost_stage_reads, [0]))
    overall_lost_stage = sum(lost_stage_reads)
    lost_stage_width = max(len(lost_stage_header), len(lost_stage_header2), len(str(max_stage)), len(str(overall_lost_stage)))
    
    headers.append(lost_stage_header)
    headers2.append(lost_stage_header2)
    header_widths.append(lost_stage_width)

    # And the correct result lost count header
    lost_header = "Lost"
    lost_header2 = ""
    # How big a number will we need to hold?
    # Look at the reads lost at all filters
    # Account for None values for stages that don't have correctness defined yet.
    lost_reads = [x for x in (stats_total[filter_name]['failed_count_correct'] for filter_name in filters) if x is not None]
    max_filter_stop = max(itertools.chain(lost_reads, [0]))
    # How many correct reads are lost overall by filters?
    overall_lost = sum(lost_reads)
    lost_width = max(len(lost_header), len(lost_header2), len(str(max_filter_stop)), len(str(overall_lost)))
    
    headers.append(lost_header)
    headers2.append(lost_header2)
    header_widths.append(lost_width)
    
    # And the total rejected count header
    rejected_header = "Cut"
    rejected_header2 = ""
    # How big a number will we need to hold?
    # Look at the reads rejected at all filters
    rejected_reads = [stats_total[filter_name]['failed_count_total'] for filter_name in filters]
    max_filter_stop = max(itertools.chain(rejected_reads, [0]))
    # How many incorrect reads are rejected overall by filters?
    overall_rejected = sum(rejected_reads)
    rejected_width = max(len(rejected_header), len(rejected_header2), len(str(max_filter_stop)), len(str(overall_rejected)))
    
    headers.append(rejected_header)
    headers2.append(rejected_header2)
    header_widths.append(rejected_width)
    
    # Now do precision and recall
    # How should we format them?
    pr_format = '{:.2f}'
    precision_header = "P"
    precision_header2 = ""
    precision_width = max(len(precision_header), len(precision_header2), len(pr_format.format(1.0)), len('N/A'))
    headers.append(precision_header)
    headers2.append(precision_header2)
    header_widths.append(precision_width)
    recall_header = "R"
    recall_header2 = ""
    recall_width = max(len(recall_header), len(recall_header2), len(pr_format.format(1.0)), len('N/A'))
    headers.append(recall_header)
    headers2.append(recall_header2)
    header_widths.append(recall_width)
    
    
    # Start the table
    table = Table(header_widths)
    
    table.row(["Giraffe Facts"], 'c', merge=[len(header_widths)])
    table.line()
    table.full_row('Reads', str(read_count))
    if 'time_used' in stats_total[NO_FILTER] and stats_total[NO_FILTER]['time_used'] != 0:
        table.full_row('Mapping speed', '{:0.2f} RPS'.format(read_count / stats_total[NO_FILTER]['time_used']))
    table.line()
    table.row(headers, 'c')
    table.row(headers2, 'c')
    table.line()
    
    
    for i, filter_name in enumerate(filters):
        # Grab average results passing this filter per read
        total_passing = stats_total[filter_name]['passed_count_total']
        average_passing = total_passing / read_count if read_count != 0 else float('NaN')
        
        # Grab average results failing this filter per read
        total_failing = stats_total[filter_name]['failed_count_total']
        average_failing = total_failing / read_count if read_count != 0 else float('NaN')
        
        # Grab reads that are lost.
        # No reads are lost at the final stage.
        lost = stats_total[filter_name]['failed_count_correct']
        
        lost_stage = stats_total[filter_name]['last_correct_stage']

        # And reads that are rejected at all
        rejected = stats_total[filter_name]['failed_count_total']
        
        if lost is None:
            # Correctness is not defined yet.
            # TODO: have a way to see if the correct mapping never shows up.
            lost = 'N/A'
            
        # Compute precision
        try:
            precision = pr_format.format(stats_total[filter_name]['passed_count_correct'] /
                stats_total[filter_name]['passed_count_total'])
        except:
            precision = 'N/A'
        
        # Compute recall
        try:
            recall = pr_format.format(stats_total[filter_name]['passed_count_correct'] / 
                (stats_total[filter_name]['passed_count_correct'] +
                stats_total[filter_name]['failed_count_correct']))
        except:
            recall = 'N/A'
        
        row = [filter_headings[i]]
        align = 'c'
        # Add the provenance columns
        row += ['{:.2f}'.format(average_passing), '{:.2f}'.format(average_failing), lost_stage, lost, rejected,
            precision, recall]
        align += 'rrrrrr'
        
        # Output the finished row
        table.row(row, align)
        
    table.line()
        
    # Compose the overall row
    row = [filter_overall]
    align = 'c'
    # Add the provenance columns
    row += ['', '', overall_lost_stage, overall_lost, overall_rejected, '', '']
    align += 'rr'
    
    table.row(row, align)
    
    # Close off table
    table.close()
    
def plot_filter_statistic_histograms(out_dir, stats_total):
    """
    For each filter in the stats dict, see if it has nonempty
    'statistic_distribution_correct' and/or 'statistic_distribution_noncorrect'
    Counters. Then if so, plot a histogram comparing correct and noncorrect
    distributions, or just the noncorrect distribution if that is the only one
    available (because correctness isn't known).
    
    Store histograms in out_dir.
    """
    
    for filter_name in stats_total.keys():
        correct_counter = stats_total[filter_name]['statistic_distribution_correct']
        noncorrect_counter = stats_total[filter_name]['statistic_distribution_noncorrect']
        
        if not ((isinstance(correct_counter, dict) and len(correct_counter) > 0) or
            (isinstance(noncorrect_counter, dict) and len(noncorrect_counter) > 0)):
            
            # No stats to plot
            continue
        
        # Open a TSV file to draw a histogram from
        tsv_path = os.path.join(out_dir, 'stat_{}.tsv'.format(filter_name))
        tsv = open(tsv_path, 'w')
        
        # Some stages don't have correctness annotation. So we track if we saw
        # correct and noncorrect things to identify them.
        have_correct = False
        have_noncorrect = False
        
        if isinstance(correct_counter, dict) and len(correct_counter) > 0:
            # We have correct item stats.
            have_correct = True
            for value, count in correct_counter.items():
                # Output format: label, value, repeats
                tsv.write('correct\t{}\t{}\n'.format(value, count))
                
        if isinstance(noncorrect_counter, dict) and len(noncorrect_counter) > 0:
            # We have noncorrect item stats.
            have_noncorrect = True
            for value, count in noncorrect_counter.items():
                # Output format: label, value, repeats
                tsv.write('noncorrect\t{}\t{}\n'.format(value, count))
              
        tsv.close()
        
        # Now make the plot
        svg_path = os.path.join(out_dir, 'stat_{}.svg'.format(filter_name))
        
        args = ['histogram.py', tsv_path, '--save', svg_path,
            '--title', '{} Statistic Histogram'.format(filter_name),
            '--x_label',  'Statistic Value',
            '--bins', '20',
            '--y_label', 'Frequency']
        if have_correct and have_noncorrect:
            args.append('--legend_overlay')
            args.append('best')
            args.append('--categories')
            args.append('correct')
            args.append('noncorrect')
        histogram.main(args)
                
            


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
    
    # Record mapping parameters from special magic GAM chunk, if any, or a read
    params = None
   
    # Get the params from a magic chunk.
    # TODO: This is a whole pass through a possibly big file!
    params_json = subprocess.check_output([options.vg, "view", "--extract-tag", "PARAMS_JSON", options.input]).decode('utf-8')
    lines = params_json.split("\n")
    for parsed_params in read_line_oriented_json(lines):
        if params is None:
            params = parsed_params

    for read in read_read_views(options.vg, options.input):
        # For the data we need on each read
        
        if params is None:
            # Go get the mapping parameters
            params = sniff_params(read)
    
        # For the stats dict for each read
        stats = make_stats(read)
        if stats_total is None:
            stats_total = stats
        else:
            # Sum up all the stats
            add_in_stats(stats_total, stats)
        
        # Count the read
        read_count += 1
    
    # After processing all the reads
    
    # Print the table now in case plotting fails
    print_table(read_count, stats_total, params)
    
    # Make filter statistic histograms
    plot_filter_statistic_histograms(options.outdir, stats_total)
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        


#!/usr/bin/python3



"""
giraffe-facts.py: process a GAM file from the new minimizer-based mapper (vg giraffe) and report runtime statistics by filter.
"""

import argparse
import concurrent.futures
import functools
import os
import sys
import time
import subprocess
import collections
import io
import itertools
import json
import random
import re
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
    parser.add_argument("--vg", "-v", type=str, default="vg",
                        help="vg binary to use")
    parser.add_argument("--threads", "-t", type=int, default=16,
                        help="total threads to use")
    parser.add_argument("--stages", "-s", action="store_true", default=False,
                        help="show stage associated with each filter")
    parser.add_argument("--filter-help", "-f", action="store_true", default=False,
                        help="show documentation for each filter")
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

class LossyCounter(collections.Counter):
    """
    A class for representing a distribution in manageable memory.

    If it tries to get too many distinct values, they start to get assigned to
    ranges instead.

    When you put an item in or get something by the item it gets coerced to a
    representative value for the range. When you ask about what is in the
    collection (or try to delete) you only see the representative items.

    Only comes up with a binning scheme at one particular point in time, so you
    can trick it into picking too-fine bins and then growing unboundedly.

    If made using an iterable, only condenses *after* iterating the iterable.
    """

    # All the internal Counter methods (in Python 3.10 at least) go through
    # items()/keys()/values() and []. But also the addition/math methods call
    # Counter() so we need to override them.

    def __init__(self, *args, **kwargs):
        """
        Make a new LossyCounter.
        """
        super().__init__(*args, **kwargs)
        # How big can we get before we start binning floats?
        self._max_len = 10000
        # What bin size are we using (or None if no bins in use yet)
        self._bin_size = None
        self._condense_if_needed()

    def _to_key(self, elem):
        """
        Map from unbinned key to the key it should actually use.
        """
        if self._bin_size is None or not isinstance(elem, float) or not math.isfinite(elem):
            # We shouldn't or can't bin this
            return elem

        # Round floats to the nearest bin using math.remainder. See
        # <https://stackoverflow.com/a/70210770>
        return elem - math.remainder(elem, self._bin_size)

    def _condense_if_needed(self):
        """
        If we are now too long, determine bin size and rehash everything in bins.
        """
        if len(self) > 10 * self._max_len and self._bin_size is not None:
            sys.stderr.write(f"Re-binning because length has reached {len(self)}\n")
            self._bin_size = None
        if self._bin_size is None and len(self) > self._max_len:
            old_len = len(self)
            old_total = sum(self.values())
            # We need to pick a bin size form our float items
            min_value = min((x for x in self.keys() if isinstance(x, float) and math.isfinite(x)), default=None)
            max_value = max((x for x in self.keys() if isinstance(x, float) and math.isfinite(x)), default=None)

            # Copy all our keys and values to temp storage, before we do
            # anything that could upset __getitem__
            all_items = list(self.items())
            # Delete them all
            self.clear()

            if min_value is not None and max_value is not None:
                # Set for this many bins
                self._bin_size = (max_value - min_value) / self._max_len
                if self._bin_size == 0:
                    # We divided too hard; the values are too close together
                    self._bin_size = min_value
                if self._bin_size == 0:
                    # Still not going to work
                    self._bin_size = max_value
            if self._bin_size is None or self._bin_size == 0:
                # If we have like a jillion strings and 0, bin by integers.
                self._bin_size = 1

            # Now that _bin_size is not None, we can't recurse into here.

            for k, v in all_items:
                # Re-insert all the items by binning. Since we're a counter we can add to aggregate.
                self[k] += v

            # Make sure nothing went wrong
            assert len(self) <= old_len
            assert sum(self.values()) == old_total


    def __getitem__(self, elem):
        """
        Get an item from the collection by key.
        """
        return super().__getitem__(self._to_key(elem))


    def __setitem__(self, elem, value):
        """
        Set an item in the collection to a value.
        """
        assigned_key = self._to_key(elem)
        super().__setitem__(assigned_key, value)
        # Now make sure we don't get too big
        self._condense_if_needed()

    def __add__(self, other):
        return LossyCounter(super().__add__(other))

    def __sub__(self, other):
        return LossyCounter(super().__sub__(other))

    def __or__(self, other):
        return LossyCounter(super().__or__(other))

    def __and__(self, other):
        return LossyCounter(super().__and__(other))

    def __pos__(self):
        return LossyCounter(super().__pos__())

    def __neg__(self):
        return LossyCounter(super().__neg__())

def flat_items(struct):
    """
    Yield dotted key strings and values for all leaf values in the given dict.
    """

    for k, v in struct.items():
        if isinstance(v, dict):
            for subkey, subvalue in flat_items(v):
                yield f"{k}.{subkey}", subvalue 
        else:
            yield k, v

def get_flat_key(struct, key, default=None):
    """
    Get a value from a nested dict by dotted key string.
    """

    if "." in key:
        first, rest = key.split(".", 1)
        return get_flat_key(struct.get(first, {}), rest, default=default)
    else:
        return struct.get(key, default)
    
# Stats under NO_FILTER are not associated with a filter
NO_FILTER = (-1, "__none__", "")

def make_stats(read):
    """
    Given a read dict parsed from JSON, compute a stats OrderedDict for the read.
    
    Run on an empty dict, makes a zero-value stats dict.
    
    A stats dict maps from filter number, name pair to a LossyCounter of filter stats.
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
    
    The LossyCounter for a filter also has sub-Counters embedded in it for
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
    associated with a filter, such as count or time used.
    """
    
    # This is the annotation dict from the read
    annot = read.get('annotation', {})

    # This will map from filter number int to filter name
    filters_by_index = {}

    # And this from filter number int to stage name
    stages_by_index = {}

    # This will map from filter name to LossyCounter of filter stats
    filter_stats = collections.defaultdict(LossyCounter)

    filter_dict = annot.get("filter", {})
    filterstats_dict = annot.get("filterstats", {})
    for filter_num_str, filter_data in filter_dict.items():
        filter_num = int(filter_num_str)
        filter_name = filter_data["name"]
        filter_stage = filter_data["stage"]

        # Record the filter being at this index if not known already
        filters_by_index[filter_num] = filter_name
        # And also the stage
        stages_by_index[filter_num] = filter_stage

        for filter_status in ("passed", "failed"):
            for filter_accounting in ("count", "size"):
                for filter_metric in ("total", "correct"):
                    # Fetch the stat value
                    filter_stat_value = filter_data.get(filter_status, {}).get(f"{filter_accounting}_{filter_metric}", None)
                    # Record the stat value
                    filter_stats[filter_name]['{}_{}_{}'.format(filter_status, filter_accounting, filter_metric)] = filter_stat_value

        # Get the stats that go with this filter (if any)
        filterstats_data = filterstats_dict.get(filter_num_str, {})

        for filter_correctness in ("correct", "noncorrect"):
            # Get the subset dict for this correctness
            filterstats_subset_data = filterstats_data.get(filter_correctness, {})
            # Get all the values and weights
            value_list = filterstats_subset_data.get("values", [])
            weight_list = filterstats_subset_data.get("weights", None)
            if weight_list is None:
                weight_list = [1] * len(value_list)

            
            distribution = LossyCounter()
            
            for item, weight in zip(value_list, weight_list):
                # Parse all the statistic vlues
                item = float(item)
                
                if math.isnan(item):
                    # Discard NANs
                    continue
                    
                # Weight the values 
                distribution[item] += weight
                
            # Save the statistic distribution
            filter_stats[filter_name]['statistic_distribution_{}'.format(filter_correctness)] = distribution

    # Now we need the filter number of the first filter in the given stage name, including "none".
    # TODO: This is all duplicate work over all reads!
    first_filter_number_in = {}
    for index, stage in stages_by_index.items():
        first_filter_number_in[stage] = min(first_filter_number_in.get(stage, float("inf")), index)

    if "last_correct_stage" in annot:
        stage = annot["last_correct_stage"]
        if stage in first_filter_number_in:
            # Assign a last correct stage point to the first filter in the named stage, which maybe lost the item.
            filter_stats[filters_by_index[first_filter_number_in[stage]]]['last_correct_stage'] = 1

    # Now put them all in this dict by number and then name, and then stage to smuggle the stage out
    ordered_stats = collections.OrderedDict()
    for filter_index in sorted(filters_by_index.keys()):
        filter_name = filters_by_index[filter_index]
        ordered_stats[(filter_index, filter_name, stages_by_index[filter_index])] = filter_stats[filter_name]
        
    # Add in special non-filter stats (not distributions)
    ordered_stats[NO_FILTER] = collections.OrderedDict()
    for k in ['time_used']:
        if k in read:
            ordered_stats[NO_FILTER][k] = read[k]
    
    # And the count
    ordered_stats[NO_FILTER]["count"] = 1 if read else 0

    return ordered_stats

def add_in_stats(destination, addend):
    """
    Add the addend stats dict into the destination stats dict.
    Implements += for stats dicts.
    """

    def add_recursive(dest, add):
        for k, v in add.items():
            if v is None:
                # None will replace anything and propagate through
                # TODO: Why???
                dest[k] = None
            elif isinstance(v, dict):
                # Recurse into dict
                if k in dest:
                    add_recursive(dest[k], v)
                else:
                    dest[k] = v
            elif k in dest:
                # Use real += and hope it works
                dest[k] += v
            else:
                # Add missing keys
                dest[k] = v

    add_recursive(destination, addend)
    return destination

def read_line_oriented_json(lines):
    """
    For each line in the given iterable of lines (such as a stream), yield it as a parsed JSON object.
    """
    
    for line in lines:
        line = line.strip()
        if len(line) > 0:
            yield json.loads(line)


def read_read_lines(vg, filename, threads=16):
    """
    Given a vg binary and a filename, iterate over nonempty TSV lines for each read in the file.

    The lines will have the annotation and time_used fields. They will be
    enumerated as bytes objects, with trailing whitespace, and some may be empty.
    """

    # Extract just the relevant annotations and times of reads as JSON, with a # header
    # Use several threads because we spend a lot of time crunching JSON, and use small batches to limit memory use.
    filter_process = subprocess.Popen([vg, "filter", f"-t{threads}", "-B10", "--tsv-out", "annotation.filter;annotation.filterstats;annotation.last_correct_stage;time_used", filename], stdout=subprocess.PIPE)

    lines = iter(filter_process.stdout)
    # Drop header line
    next(lines)
    
    # Yield the rest.
    yield from lines

    return_code = filter_process.wait()
    assert return_code == 0

def read_line_to_read_dict(line):
    """
    Given a bytes object that may be a blank line from read_read_lines(),
    return a read dict of the subset of the read information we loaded.

    The result will have the annotation.filterm annotation.filterstats,
    annotation.last_correct_stage, and time_used fields, or be empty and
    represent a blank line.
    """

    # Parse the TSV and reconstruct a view of the full read dict.
    line = line.decode('utf-8')
    line = line.strip()
    if len(line) == 0:
        return {}
    parts = line.split("\t")
    assert len(parts) == 4
    filter_json, filterstats_json, last_correct_stage, time_used_string = parts
    read = {
        "annotation": {
            "filter": json.loads(filter_json) or {},
            "filterstats": json.loads(filterstats_json) or {},
            "last_correct_stage": last_correct_stage
        },
        "time_used": float(time_used_string) if time_used_string != "null" else float("nan")
    }

    return read

def read_line_to_stats(line):
    """
    Map from a read TSV line (possibly empty) to a stats dict.
    """

    return make_stats(read_line_to_read_dict(line))

def batched(input_iterator, batch_size):
    """
    Yield batch lists of up to batch_size items from the given iterator.
    """
    
    if hasattr(itertools, "batched"):
        # Itertools implements this
        return itertools.batched(input_iterator, batch_size)
    else:
        # Implement it ourselves.
        # This is not the same as any other object. It will fill the empty spots in the batches.
        sentinel = object()
        # We can zip an iterator with itself to take batches, like in grouper()
        # at
        # <https://docs.python.org/3/library/itertools.html#itertools-recipes>.
        # Then we just have to trim down the extended batches.
        return ([x for x in batch if x is not sentinel] for batch in itertools.zip_longest(*([iter(input_iterator)] * batch_size), fillvalue=sentinel))

def map_reduce(input_iterable, map_function, reduce_function, zero_function, threads=16):
    """
    Do a parallel map-reduce operation on the input iterable.

    Return the reduction result of all the map_function results reduced by the reduce_function.

    The reduce_function may modify either argument in place but must return the reduction result.

    If there are no items, returns the result of the zero_function.
    """

    if threads <= 1:
        # Do it all in one thread
        reduced = zero_function()
        for item in input_iterable:
            reduced = reduce_function(reduced, map_function(item))
        return reduced
   
    MAP_CHUNK_SIZE = 100
    REDUCE_CHUNK_SIZE = 1000
    max_tasks = threads * 2
    
    executor = concurrent.futures.ProcessPoolExecutor(threads)

    # Make sure we have an iterator and not like a list
    input_iterator = iter(input_iterable)
    # Then bvatch it
    input_batches = batched(input_iterator, MAP_CHUNK_SIZE)

    # We can't actually use the executor's map() since it loads the whole input
    # iterable into memory to start. Instead we do some futures wizardry. See
    # <https://alexwlchan.net/2019/adventures-with-concurrent-futures/>.

    # Futures for map tasks in progress
    map_in_flight = {executor.submit(map_function, item) for item in itertools.islice(input_iterator, max_tasks)}
    # Futures for reduce tasks in progress
    reduce_in_flight = set()

    # Reduce task currently being buffered
    reduce_buffer = []

    while map_in_flight or reduce_in_flight:
        # Until all the work is done

        if reduce_in_flight:
            # Only wait for reduces to complete if no maps are running.
            # Otherwise just poll.
            timeout = 0 if map_in_flight else None
            # Get any completed reduces.
            reduced, reduce_in_flight = concurrent.futures.wait(reduce_in_flight, timeout=0, return_when=concurrent.futures.FIRST_COMPLETED)
            for f in reduced:
                # Put their results in the buffer to be reduced
                reduce_buffer.append(f.result())
                if len(reduce_buffer) >= REDUCE_CHUNK_SIZE:
                    # Kick off a reduction.
                    # We won't make more than one new job per finished job here.
                    reduce_in_flight.add(executor.submit(functools.reduce, reduce_function, reduce_buffer))
                    reduce_buffer = []
                    

        if map_in_flight:
            # Get at least one done mapping result
            mapped, map_in_flight = concurrent.futures.wait(map_in_flight, return_when=concurrent.futures.FIRST_COMPLETED)
            
            for f in mapped:
                # Put their results in the buffer to be reduced
                reduce_buffer.append(f.result())
                if len(reduce_buffer) >= REDUCE_CHUNK_SIZE:
                    # Kick off a reduction.
                    # We won't make more than one new job per finished job here.
                    reduce_in_flight.add(executor.submit(functools.reduce, reduce_function, reduce_buffer))
                    reduce_buffer = []
            
        new_tasks = max_tasks - len(map_in_flight) - len(reduce_in_flight)
        if new_tasks > 0:
            # Top up mapping to have at most the set number of jobs in flight
            map_in_flight |= {executor.submit(map_function, item) for item in itertools.islice(input_iterator, new_tasks)}

    # When everything is done, finish reduction in this thread
    result = functools.reduce(reduce_function, reduce_buffer, zero_function())

    return result
        
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
                
def print_table(stats_total, params=None, out=sys.stdout, show_stages=False):
    """
    Take the accumulated total stats dict, and an optional dict
    of mapping parameters corresponding to values for filters.
    
    Print a nicely formatted table to the given stream.
    
    """
    
    if not stats_total:
        # Handle the empty case
        out.write('No reads.\n')
        return
    
    read_count = stats_total.get(NO_FILTER, {}).get("count", 0)

    # Now do a table
    
    # First header line for each column
    headers = []
    # Second header line for wach column
    headers2 = []
    # Column min widths from headers
    header_widths = []
    
    # Find all the filter value, name pairs
    filters = sorted([k for k in stats_total.keys() if k != NO_FILTER])
    
    # Compute filter row headings
    filter_headings = [k[1] for k in filters]
    # And their stages, which we might show
    filter_stages = [k[2] for k in filters]

    if show_stages:
        # How long is the longest stage name
        stage_width = max(itertools.chain((len(x) for x in filter_stages), [0]))
        # Leave room for the header.
        # Filters are given the stage they filter the items of, so label the
        # stage names accordingly.
        stage_header = "Item"
        stage_width = max(stage_width, len(stage_header))
        
        headers.append(stage_header)
        headers2.append('')
        header_widths.append(stage_width)

    
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
    lost_stage_reads = [x for x in (stats_total[filter_key].get('last_correct_stage', 0) for filter_key in filters) if x is not None]
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
    lost_reads = [x for x in (stats_total[filter_key]['failed_count_correct'] for filter_key in filters) if x is not None]
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
    rejected_reads = [stats_total[filter_key]['failed_count_total'] for filter_key in filters]
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

    # Track the previous stage so we can mark off when stages change.
    prev_stage = None
    
    for i, filter_key in enumerate(filters):
        # Grab average results passing this filter per read
        total_passing = stats_total[filter_key]['passed_count_total']
        average_passing = total_passing / read_count if read_count != 0 else float('NaN')
        
        # Grab average results failing this filter per read
        total_failing = stats_total[filter_key]['failed_count_total']
        average_failing = total_failing / read_count if read_count != 0 else float('NaN')
        
        # Grab reads that are lost.
        # No reads are lost at the final stage.
        lost = stats_total[filter_key]['failed_count_correct']
        
        lost_stage = stats_total[filter_key]['last_correct_stage']

        # And reads that are rejected at all
        rejected = stats_total[filter_key]['failed_count_total']
        
        if lost is None:
            # Correctness is not defined yet.
            # TODO: have a way to see if the correct mapping never shows up.
            lost = 'N/A'
            
        # Compute precision
        try:
            precision = pr_format.format(stats_total[filter_key]['passed_count_correct'] /
                stats_total[filter_key]['passed_count_total'])
        except:
            precision = 'N/A'
        
        # Compute recall
        try:
            recall = pr_format.format(stats_total[filter_key]['passed_count_correct'] / 
                (stats_total[filter_key]['passed_count_correct'] +
                stats_total[filter_key]['failed_count_correct']))
        except:
            recall = 'N/A'

        row = []
        align = ''

        if show_stages:
            if prev_stage is None or prev_stage != filter_stages[i]:
                # Write the stage name once per stage
                row.append(filter_stages[i])
                if prev_stage is not None:
                    # When the stage changes, put a line
                    table.line()
                prev_stage = filter_stages[i]
            else:
                # Leave a blank in the stage column
                row.append('')
            align += 'l'
        
        row.append(filter_headings[i])
        align += 'c'
        # Add the provenance columns
        row += ['{:.2f}'.format(average_passing), '{:.2f}'.format(average_failing), lost_stage, lost, rejected,
            precision, recall]
        align += 'rrrrrr'
        
        # Output the finished row
        table.row(row, align)
        
    table.line()
        
    # Compose the overall row
    row = []
    align = ''

    if show_stages:
        row.append('')
        align += 'l'

    row.append(filter_overall)
    align += 'c'
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
    
    for filter_key in (k for k in sorted(stats_total.keys()) if k != NO_FILTER):
        correct_counter = stats_total[filter_key]['statistic_distribution_correct']
        noncorrect_counter = stats_total[filter_key]['statistic_distribution_noncorrect']
        
        if not ((hasattr(correct_counter, "items") and len(correct_counter) > 0) or
            (hasattr(noncorrect_counter, "items") and len(noncorrect_counter) > 0)):
            
            # No stats to plot
            continue
        
        # Open a TSV file to draw a histogram from
        tsv_path = os.path.join(out_dir, 'stat_{}.tsv'.format(filter_key[1]))
        tsv = open(tsv_path, 'w')
        
        # Some stages don't have correctness annotation. So we track if we saw
        # correct and noncorrect things to identify them.
        have_correct = False
        have_noncorrect = False
        
        if hasattr(correct_counter, "items") and len(correct_counter) > 0:
            # We have correct item stats.
            have_correct = True
            for value, count in correct_counter.items():
                # Output format: label, value, repeats
                tsv.write('correct\t{}\t{}\n'.format(value, count))
                
        if hasattr(noncorrect_counter, "items") and len(noncorrect_counter) > 0:
            # We have noncorrect item stats.
            have_noncorrect = True
            for value, count in noncorrect_counter.items():
                # Output format: label, value, repeats
                tsv.write('noncorrect\t{}\t{}\n'.format(value, count))
              
        tsv.close()
        
        # Now make the plot
        svg_path = os.path.join(out_dir, 'stat_{}.svg'.format(filter_key[1]))
        
        args = ['histogram.py', tsv_path, '--save', svg_path,
            '--title', '{} Statistic Histogram'.format(filter_key[1]),
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

def capitalize_first(string):
    """
    Capitalize the first letter of a string while elaving the other letters alone.
    """

    return string[0].upper() + string[1:]

def explain_filters(filters_and_stages, vg):
    """
    Go get the vg giraffe help for all the filters in filter_names, in order, and print it.
    """

    # Get all the lines of help text
    help_lines = subprocess.run([vg, "giraffe", "--help"], stderr=subprocess.PIPE, check=False).stderr.decode('utf-8').split('\n')

    # This will hold help line lists by filter name
    filter_helps = collections.defaultdict(list)

    # Add the magic built-in filters with no help
    filter_helps["window-downsampling"].append("Downsample minimizers to one per read-length-based sliding window")
    filter_helps["any-hits"].append("Drop minimizers without any hits in the graph")
    filter_helps["no-chain-overlap"].append("Drop chains that share read-graph matchings with previously aligned chains")


    # Extract the names without the stages
    filter_names = [f[0] for f in filters_and_stages]

    # These are the filters we are still looking for.
    unseen_filters = set(filter_names)
    
    # And this is the filter whole help we are currently reading
    current_filter = None
    for line in help_lines:
        # Do a state machine down the help lines.
        if current_filter is not None:
            # We're reading a help section
            if not line.startswith(' ') or line.lstrip().startswith('-'):
                # Starting a new header or argument.
                current_filter = None
            else:
                # Continuing the filter help for the previous option.
                # Add on the line
                filter_helps[current_filter].append(line.strip())
        if current_filter is None:
            # Look for all the filters we haven't started yet on this line,
            # longest first to handle filters which are substrings of each
            # other
            sorted_filters = sorted(unseen_filters, key=len, reverse=True)
            for filter_to_try in sorted_filters:
                if filter_to_try in line:
                    # This line starts the section for this filter
                    current_filter = filter_to_try
                    filter_helps[current_filter].append(line.rstrip())
                    unseen_filters.remove(filter_to_try)
                    break


    for filter_name, help_lines in filter_helps.items():
        # Rewrite the filter option specs out of the first lines and capitalize it.
        help_lines[0] = capitalize_first(re.sub("^ +(-[a-zA-Z0-9], )?--[-a-zA-Z0-9]+( [A-Z]+)? +", "", help_lines[0]))
        # And the trailing defaults
        help_lines[-1] = re.sub("\[[0-9.]+\]$", "", help_lines[-1])

    # Now print the results by stage
    filters_by_stage = collections.defaultdict(list)
    for filter_name, stage in filters_and_stages:
        filters_by_stage[stage].append(filter_name)
    
    for stage in filters_by_stage.keys():
        print(f"Filters on {stage} items:")
        for filter_name in filters_by_stage[stage]:
            filter_help_string = ' '.join(filter_helps[filter_name]) if filter_name in filter_helps else "!!!UNKNOWN!!!"
            print(f"    * {filter_name}: {filter_help_string}")

def explain_filters_in(stats_total, vg):
    """
    Print an explanation of all the filters in the given stats dict, in order.
    """

    filters_and_stages = []
    for k in stats_total.keys():
        if k != NO_FILTER:
            for piece in k[1].split('||'):
                # Make a list of all filter parameter anmes and their stages
                filters_and_stages.append((piece, k[2]))

    explain_filters(filters_and_stages, vg)

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
    params_json = subprocess.check_output([options.vg, "view", "--extract-tag", "PARAMS_JSON", "--first", options.input]).decode('utf-8')
    lines = params_json.split("\n")
    for parsed_params in read_line_oriented_json(lines):
        if params is None:
            params = parsed_params

    # Get the stream of TSV lines for reads
    read_lines = read_read_lines(options.vg, options.input, threads=max(1, options.threads // 2))
   
    # Map it to stats and reduce it to total stats
    stats_total = map_reduce(read_lines, read_line_to_stats, add_in_stats, collections.OrderedDict, threads=max(1, options.threads // 2))

    # After processing all the reads
    
    # Print the table now in case plotting fails
    print_table(stats_total, params, show_stages=options.stages)

    if options.filter_help:
        # Explain all the filters
        explain_filters_in(stats_total, options.vg)
    
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
        


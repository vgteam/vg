#!/usr/bin/python3
# -*- coding: utf-8 -*-
# filter_variants_on_repeats.py: Drop low-frequency variants on repeats
"""

Low-frequency variants on repeats might be counterproductive in a graph.

They can be more likely to be matched spuriously due to an error than due to a
genuine instance of the variant.

This script removes such variants from VCFs.

Repeat count information is pulled from the hgsql command, which must be
available.

"""

import argparse
import os
import sys
import time
import subprocess
import collections

import vcf
import tsv
from intervaltree import IntervalTree
import itertools

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
    
    parser.add_argument("vcf", type=argparse.FileType('r'), default=sys.stdin,
        help="VCF file to filter")
    parser.add_argument("--out_vcf", type=argparse.FileType('w'),
        default=sys.stdout,
        help="VCF file to write passing entries to")
    parser.add_argument("--assembly", default="hg38",
        help="assembly database name to query")
    parser.add_argument("--contig", default="chr17",
        help="contig to pull variants from")
    parser.add_argument("--start", type=int, default=43044294,
        help="start position on contig")
    parser.add_argument("--end", type=int, default=43125484,
        help="end position on contig")
    parser.add_argument("--error_rate", type=float, default=0.01,
        help="error rate to assume for sequencing reads")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

class RepeatDb(object):

    def __init__(self, assembly, contig, start, end):
        """
        Given a range on a contig, get all the repeats overlapping that range.
        
        Keeps an IntervalTree of element names, and a Counter from element
        name to number of that element in the range.
        
        No protection against SQL injection.
        
        """
        
        # Make the interval tree
        self.tree = IntervalTree()
        
        # Make a counter for repeats with a certain name
        self.counts = collections.Counter()
        
        command = ["hgsql", "-e", "select repName, genoName, genoStart, genoEnd "
            "from {}.rmsk where genoName = '{}' and genoStart > '{}' "
            "and genoEnd < '{}';".format(assembly, contig, start, end)]
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        
        for parts in itertools.islice(tsv.TsvReader(process.stdout), 1, None):
            # For each line except the first, broken into fields
            
            # Add the item to the tree covering its range. Store the repeat type
            # name as the interval's data.
            self.tree.addi(int(parts[2]), int(parts[3]), parts[0])
            
            # Count it
            self.counts[parts[0]] += 1

    def get_copies(self, contig, pos):
        """
        Given a contig name and a position, estimate the copy number of that
        position in the genome.
        
        Return the number of instances expected (1 for non-repetitive sequence).
        """
        
        # TODO: use contig
        
        # Get the set of overlapping things
        overlaps = self.tree[pos]
        
        # Keep track of the number of copies of the most numerous repeat
        # observed.
        max_copies = 1
        
        for interval in overlaps:
            # For each repeat we are in
            
            # Max in how many copies of it exist
            max_copies = max(max_copies, self.counts[interval.data])
            
        return max_copies
            
   
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Make a repeat database
    db = RepeatDb(options.assembly, options.contig, options.start, options.end)
    
    # Strip "chr" from the contig.
    # TODO: figure out if we need to
    short_contig = (options.contig[3:] if len(options.contig) > 3
        else options.contig)
    
    # Read the input VCF
    reader = vcf.Reader(options.vcf)
    
    # Make a writer for the passing records
    writer = vcf.Writer(options.out_vcf, reader)
    
    # Track statistics
    total_variants = 0
    kept_variants = 0
    non_snp_variants = 0
    
    for variant in reader.fetch(short_contig, options.start, options.end):
        # For each variant in our region
        
        # Count it
        total_variants += 1
        
        # See how numerous it should be
        count = db.get_copies("chr" + variant.CHROM, variant.POS)
        
        # And how common it is
        frequency = variant.INFO["AF"][0]
        
        if variant.is_snp:
            # SNPs are modeled
        
            # What's the minimum frequency to get more true than fake hits?
            min_frequency = (count * options.error_rate /
                (3 - 2 * options.error_rate))
                
            # Should we keep it?
            keep = (frequency >= min_frequency)
            
            sys.stderr.write(
                "{} {}:{} has {} copies, frequency {:.4f} vs. {:.4f}\n".format(
                "☑" if keep else "☐", variant.CHROM, variant.POS, count,
                frequency, min_frequency))
        else:
            # We have to keep everything that's not a SNP because we have no
            # error model
            
            keep = True
            sys.stderr.write(
                "{} {}:{} has {} copies, frequency {:.4f} (non-SNP)\n".format(
                "☑" if keep else "☐", variant.CHROM, variant.POS, count,
                frequency))
            non_snp_variants += 1
            
        if keep:
            # Pass the variants we're keeping
            writer.write_record(variant)
            # Count it
            kept_variants += 1
            
    sys.stderr.write("Finished! Kept {} ({} non-SNP) / {} variants\n".format(
        kept_variants, non_snp_variants, total_variants))
        
    
    return 0
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        
        
        
        
        
        
        
        
        
        


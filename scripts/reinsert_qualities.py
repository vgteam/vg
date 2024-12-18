# Stephen Hwang's FASTQ quality inserter into SAM files.
# Adds qualities from a FASTQ back into a SAM that is missing them.
# License: "I can put it online this afternoon or go ahead" - Stephen Hwang
# https://ucsc-gi.slack.com/archives/D02GGLLQXUM/p1673976340012069

import re
import sys
from math import log
from statistics import stdev


class FastAreader:
    """
    Class to contain the necessary methods to parse out fasta files. Reads fasta files either from filenames passed into the class, or from STDIN.

    Author: David Bernick
    Initialized: filename that is either passed in to the class or an empty string
    Methods: doOpen(): either reads in STDIN or opens the file to read its lines, readFasta(): parses the fasta file, separates the actual sequence from the header, removes the newline characters, and yields a generator
    """
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        """ Return input from either  STDIN or filename """
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):
        """ Return generator after filtering out header, cleaning newlines, and whitespace from sequence """
        header = ''
        sequence = ''
        # open the file to read its lines
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            # if the line doesn't start with > it is a sequence
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                # join together sequences under the same header
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header,sequence


class FastQreader :
    """
    Class to contain the necessary methods to parse out fasta files. Reads fasta files either from filenames passed into the class, or from STDIN.

    Author: David Bernick
    Initialized: filename that is either passed in to the class or an empty string
    Methods: doOpen(): either reads in STDIN or opens the file to read its lines, readFasta(): parses the fasta file, separates the actual sequence from the header, removes the newline characters, and yields a generator
    """
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        """ Return input from either  STDIN or filename """
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFastq (self):
        """ Return generator after filtering out header, cleaning newlines, and whitespace from sequence """
        header = ''
        sequence = ''
        # open the file to read its lines

        # print('starting reading')
        read_num = 1

        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            score = ''
            on_sequence = True
            line = fileH.readline().strip()

            # skip to first fasta header
            while not line.startswith('@'):
                line = fileH.readline()
            header = line[1:].rstrip()

            all_header = header.split('_')[0] + '_'
            # print('all_header', all_header)
            # print ('on reads')

            for line in fileH:
                # if the line doesn't start with @ it is a sequence or score
                # print(line)
                # print('@' + all_header + str(read_num))


                # if line.startswith('@' + all_header + str(read_num)):         # @S#_
                if line.startswith('@' + all_header):         # @S#_
                    # print(header, read_num)

                # if line.startswith ('@S'):         # @S#_
                # if re.match(r'^@S\d_\d]', line):
                    # print('match')
                    yield header, sequence, score
                    read_num += 1
                    header = line[1:].rstrip()
                    sequence = ''
                    score = ''
                    on_sequence = True
                # join together sequences under the same header
                else:
                    # print('no match')
                    # if line.strip() != '+':
                    # if not line.strip().startswith('+S'):
                    if not line.strip().startswith('+' + all_header):
                    # if not re.match(r'^\+S\d_\d]', line):
                        if on_sequence:
                            sequence += ''.join(line.rstrip().split()).upper()
                        else:
                            score += ''.join(line.rstrip().split()).upper()
                    elif on_sequence:
                        on_sequence = False

        yield header,sequence,score



class SAMreader:        # assumes everything on a single line
    """
    Class to contain the necessary methods to parse out fasta files. Reads fasta files either from filenames passed into the class, or from STDIN.

    Author: David Bernick
    Initialized: filename that is either passed in to the class or an empty string
    Methods: doOpen(): either reads in STDIN or opens the file to read its lines, readFasta(): parses the fasta file, separates the actual sequence from the header, removes the newline characters, and yields a generator
    """
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        """ Return input from either  STDIN or filename """
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFile(self):
        ''' Read file line-by-line. '''
        for line in self.doOpen():
            yield line.strip()

    def readSAM(self):
        ''' Parse HMM rosalind file into x, alphabet, path, and states. '''
        global_headers = []

        lines = self.readFile()
        next_line = next(lines)

        while next_line.startswith('@'):
            global_headers.append(next_line)
            next_line = next(lines)
        print('\n'.join(global_headers))        # print SAM header lines

        # now on sequence: then continue to end
        yield next_line
        for line in lines:
            yield line


def reverseComplement(seq):
    ''' Return reverse complement of a sequence. '''
    complement = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G'}
    return ''.join([complement.get(base, 'N') for base in seq.upper()[::-1]])





class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''
    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] >output')

        self.parser.add_argument('-s', '--sam', action='store', nargs='?',
                                 required=True, help='fastq (not compressed)')
        self.parser.add_argument('-f', '--fastq', action='store', nargs='?',
                                 required=True, help='maf file')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)








################################################################################

def main():
    '''

python sam_reader.py -s /public/groups/vg/sjhwang/vg_scripts/bin/reads/sim_HiFi_other_tools/sim_pbsim2/sim_NA19239/sim_NA19239/sam/tmp/head.sam \
                     -f /public/groups/vg/sjhwang/vg_scripts/bin/reads/sim_HiFi_other_tools/sim_pbsim2/sim_NA19239/sim_NA19239/sam/tmp/head.fastq \
                > sam_with_quality.sam
    '''

    # sam_file_path = '/public/groups/vg/sjhwang/vg_scripts/bin/reads/sim_HiFi_other_tools/sim_pbsim2/sim_NA19239/sim_NA19239/sam/tmp/head.sam'
    # fastq_file_path = '/public/groups/vg/sjhwang/vg_scripts/bin/reads/sim_HiFi_other_tools/sim_pbsim2/sim_NA19239/sim_NA19239/sam/tmp/head.fastq'
    thisCommandLine = CommandLine()
    sam_file_path = thisCommandLine.args.sam
    fastq_file_path = thisCommandLine.args.fastq

    sam_obj = SAMreader(sam_file_path)
    fastq_obj = FastQreader(fastq_file_path)

    for fastq_line, sam_line in zip(fastq_obj.readFastq(), sam_obj.readSAM()):
        fastq_header, fastq_sequence, fastq_score = fastq_line

        # print(sam_line.split('\t'))
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, other = sam_line.split('\t')

        # make sure header and sequence is the same and the len of score is the length of sequence
        if fastq_header != qname:
            print('oh no: header', fastq_header)
            break

        if fastq_sequence.upper() != seq.upper():
            if reverseComplement(fastq_sequence.upper()) == seq.upper():
                fastq_score = fastq_score[::-1]
            else:
                print('oh no: sequence', fastq_header)
                print(fastq_sequence)
                print(seq)
                break

        if len(seq) != len(fastq_score):
            print('oh no: length', fastq_header)
            print(len(seq), len(fastq_score))
            # sam_line = [qname, fastq_score]
            # print('\t'.join(sam_line))
            break

        # print(header)
        # print(sequence)
        # print('score', score)
        # print(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, other)
        sam_line = [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, fastq_score, other]
        print('\t'.join(sam_line))




main()

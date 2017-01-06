#!/usr/bin/python

import sys

threshold = 100

for line in sys.stdin:
    fields = line.split(' ')
    aln_name = fields[0]
    true_chr = fields[1]
    true_pos = int(fields[2])
    aln_chr = fields[3]
    aln_pos = int(fields[4])
    aln_mapq = int(fields[5])
    aln_correct = 1 if aln_chr == true_chr and abs(true_pos - aln_pos) < threshold else 0
    print aln_name, aln_correct, aln_mapq

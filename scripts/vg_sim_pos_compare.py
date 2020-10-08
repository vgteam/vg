#!/usr/bin/python3

### positional comparison script useful for tagging alignments with their true position

# length.bp       unaligned.bp    known.nodes     known.bp        novel.nodes     novel.bp

import sys

threshold = int(sys.argv[1])

for line in sys.stdin:
    fields = line.split(' ')
    # every input has a true position, and if it has less than the expected number of fields we assume alignment failed
    aln_name = fields[0]
    if len(fields) != 13:
        print(aln_name, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        continue
    aln_chr = fields[1]
    aln_pos = int(fields[2])
    aln_mapq = int(fields[3])
    aln_score = 0
    try: aln_score = int(fields[4])
    except: pass
    true_chr = fields[5]
    true_pos = int(fields[6])
    length = int(fields[7])
    unaligned = int(fields[8])
    known_nodes = int(fields[9])
    known_bp = int(fields[10])
    novel_nodes = int(fields[11])
    novel_bp = int(fields[12])
    aln_correct = 1 if aln_chr == true_chr and abs(true_pos - aln_pos) < threshold else 0
    print(aln_name, aln_correct, aln_mapq, aln_score, length, unaligned, known_nodes, known_bp, novel_nodes, novel_bp)

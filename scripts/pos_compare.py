#!/usr/bin/python

### positional comparison script useful for tagging alignments with their true position
## for use with results of this form
# vg map -iG x.paired.sim -x x.xg -g x.gcsa \
#     | vg annotate -x x.xg -p -a - \
#     | vg view -a - \
#     | jq -c  '[.name, .refpos[0].name, .refpos[0].offset, if .mapping_quality == null then 0 else .mapping_quality end ] | @sh'  \
#     | sed 's/"//g' | sed "s/'//g" | sed s/null/0/g \
#     | sort >x.pos
## then merged and passed through this script
# join x.true_pos x.pos | ./pos_compare.py
## the output is a table of aln_name, correct (0/1), aln_mapq
## the same result can be examined for a SAM-capable aligner like bwa mem:
# bwa mem x.fasta x.fq.gz | cut -f 1,3,4,5 | grep -v ^@ | sort >x.bwa_pos

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

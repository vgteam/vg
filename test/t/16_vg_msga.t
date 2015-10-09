#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg


plan tests 8

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -t 1 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 549b183e90036f767acdd10e7d5ba125 "MSGA produces the expected graph for GRCh38 HLA-V"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -t 4 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) 549b183e90036f767acdd10e7d5ba125 "graph for GRCh38 HLA-V is unaffected by the number of alignment threads"

is $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 10000 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ ) $(vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -k 16 -n 'gi|568815592:29791752-29792749' -n 'gi|568815454:1057585-1058559' -b 'gi|568815454:1057585-1058559' -B 300 | vg view - | grep ^S | cut -f 3 | sort | md5sum | cut -f 1 -d\ )  "varying alignment bandwidth does not affect output graph"

is $(vg msga -f msgas/s.fa -k 16 -b s1 -B 20 | vg view - | md5sum | cut -f 1 -d\ ) ddd2fc488ac58a95e932d0eabf1be800 "banded alignment can detect and include large deletions in the graph"

is $(vg msga -f msgas/s.fa -k 16 -b s2 -B 20 | vg view - | md5sum | cut -f 1 -d\ ) ddd2fc488ac58a95e932d0eabf1be800 "assembly around indels works regardless of the base sequence that is used"

vg construct -v tiny/tiny.vcf.gz -r tiny/tiny.fa >t.vg
is $(vg msga -g t.vg -s CAAATTTTCTGGAGTTCTAT | vg stats -s - | wc -l) 1 "soft clips at node boundaries (start) are included correctly"
is $(vg msga -g t.vg -s TTCTATAATATG | vg stats -s - | wc -l) 1 "soft clips at node boundaries (end) are included correctly"
rm t.vg

#is $(vg msga -f msgas/s.fa -k 16 -b s1 -B 20 | vg view -j - | jq -M --sort-keys .) '{"edge":[{"from":1,"to":2},{"from":1,"to":5},{"from":2,"to":3},{"from":3,"to":4},{"from":4,"to":6},{"from":5,"to":6}],"node":[{"id":1,"sequence":"TCAGATTCTCATCCCTCCTCAAGGG"},{"id":2,"sequence":"CTTCTAACTACTCCACATCAAAGCT"},{"id":3,"sequence":"ACCCAGGCCATTTTAAGTTTCCTGT"},{"id":4,"sequence":"GGACT"},{"id":5,"sequence":"CTTCTGGACT"},{"id":6,"sequence":"AAGGACAAAGGTGCGGGGAG"}],"path":[{"mapping":[{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":1}},{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":2}},{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":3}},{"edit":[{"from_length":5,"to_length":5}],"position":{"node_id":4}},{"edit":[{"from_length":20,"to_length":20}],"position":{"node_id":6}}],"name":"s1"},{"mapping":[{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":1}},{"edit":[{"from_length":10,"to_length":10}],"position":{"node_id":5}},{"edit":[{"from_length":20,"to_length":20}],"position":{"node_id":6}}],"name":"s2"}]}' "banded alignment can detect and include large deletions in the graph"

#is $(vg msga -f msgas/s.fa -k 16 -b s2 -B 20 | vg view -j - | jq -M -c --sort-keys .) '{"edge":[{"from":1,"to":2},{"from":1,"to":4},{"from":2,"to":3},{"from":3,"to":6},{"from":4,"to":5},{"from":5,"to":6},{"from":6,"to":7}],"node":[{"id":1,"sequence":"TCAGATTCTCATCCCTCCTCAAGGG"},{"id":2,"sequence":"CT"},{"id":3,"sequence":"TCT"},{"id":4,"sequence":"CTTCTAACTACTCCACATCAAAGCT"},{"id":5,"sequence":"ACCCAGGCCATTTTAAGTTTCCTGT"},{"id":6,"sequence":"GGACTAAGGACAAAGGTGCGGGGA"},{"id":7,"sequence":"G"}],"path":[{"mapping":[{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":1}},{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":4}},{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":5}},{"edit":[{"from_length":24,"to_length":24}],"position":{"node_id":6}},{"edit":[{"from_length":1,"to_length":1}],"position":{"node_id":7}}],"name":"s1"},{"mapping":[{"edit":[{"from_length":25,"to_length":25}],"position":{"node_id":1}},{"edit":[{"from_length":2,"to_length":2}],"position":{"node_id":2}},{"edit":[{"from_length":3,"to_length":3}],"position":{"node_id":3}},{"edit":[{"from_length":24,"to_length":24}],"position":{"node_id":6}},{"edit":[{"from_length":1,"to_length":1}],"position":{"node_id":7}}],"name":"s2"}]}' "assembly around indels works regardless of the base sequence that is used"

is $(vg msga -f msgas/s.fa -k 16 -b s1 -B 20 | vg view -j - | jq -M -c --sort-keys '{"node": .node, "edge": .edge}') $(vg msga -f msgas/s.fa -f msgas/s-rev.fa -k 16 -b s1 -B 20 | vg view -j - | jq -M -c --sort-keys '{"node": .node, "edge": .edge}') "adding in existing sequences in reverse doesn't change graph"

# TODO: you should get the same graph no matter the orientation of the input sequences, but you don't. This is probably due to the orientation-dependence of the kmer search.

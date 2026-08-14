The graphs here were produced by `vg msga`, a subcommand which was deprecated
in a 2020 commit (#3136) and has now (vg 1.76.0) been entirely removed.

To update these graphs you'd need to get an old enough version of vg. You would
also need the input `.fa` files which used to be in here. The following commands
were used to create the graphs which exist here now:

- `vg msga -f cyclic/cycle.fa -b s1 -w 64 -t 1 > msgas/c1.vg`
- `vg msga -g <(vg msga -f cyclic/cycle.fa -b s1 -w 19 -O 18 -k 4 -t 1 | vg paths -d -v - | vg mod -U 10 -) -f cyclic/cycle.fa -t 1 > msgas/c2.vg`
- `cat tiny/tiny.fa tiny/flat-s69-n1-l50-e0.05.fa tiny/flat-s77-n1-l50-e0.05.fa > flats.fa` and `vg msga -f flats.fa -b x > msgas/flat.vg`
- `vg msga -f GRCh38_alts/FASTA/HLA/B-3106.fa -D > msgas/hla_b.vg`
- `vg msga -f GRCh38_alts/FASTA/HLA/K-3138.fa -B 256 -k 22 -K 11 -X 1 -E 4 -Q 22 -D > msgas/hla_k.vg`
- `vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 > msgas/hla_v.vg`
- `vg msga -f msgas/l.fa -b a1 -w 16 > msgas/l.vg`
- `vg msga -w 20 -f msgas/q.fa > msgas/q.vg`
- `vg msga -w 20 -f msgas/s.fa > msgas/s.vg`
- `vg msga -g msgas/s.vg -s TCAGATTCTCATCCCTCCTCAAGGGCTTCTGTAGCTTTGATGTGGAGTAGTTCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAG -w 16 -N > msgas/s2.vg`
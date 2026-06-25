The graphs here were produced by `vg msga`, a subcommand which was deprecated
in a 2020 commit (#3136) and has now (vg 1.76.0) been entirely removed.

To update these graphs you'd need to get an old enough version of vg. You would
also need the input `.fa` files which used to be in here. The following commands
were used to create the graphs which exist here now:

- `vg msga -f msgas/l.fa -b a1 -w 16 > msgas/l.vg`
- `vg msga -w 20 -f msgas/q.fa > msgas/q.vg`
- `vg msga -w 20 -f msgas/s.fa > msgas/s.vg`
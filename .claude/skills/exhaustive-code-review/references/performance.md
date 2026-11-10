# Performance lens

Ask only: does this make vg slower or use more memory in a way that matters?

vg runs on whole human pangenomes: graphs with hundreds of millions of nodes,
tens of millions of reads. A cost per node, per edge, per read or per snarl is
multiplied by those numbers.

- **Complexity.** A search, `find`, or path walk inside a loop over nodes,
  steps or reads. Walking a whole path to find one step. Building a
  set or map per item when one built once would do.
- **Allocation in hot loops.** A `std::string`, `std::vector` or `std::function`
  made per base, per step, or per read in a mapping inner loop. `get_sequence`
  copies; prefer `get_base`/`get_subsequence` or `for_each_base` where only part
  is needed.
- **Copies.** Passing large objects by value, copying a graph, or returning a
  container that is then copied again.
- **Memory.** A per-node table on a whole-genome graph: check its size at
  100 million nodes. Holding all reads in memory when they could stream.
- **I/O.** Reading a file twice, unbuffered writes, flushing (`endl`) per
  record in output loops.
- **Parallelism.** Work moved out of, or a lock added inside, a parallel
  region. A critical section around the work itself rather than its output.

Say how big the cost is (per what, times how many) before ranking it. A cost
that only runs once per program is not a finding. Avoidable work done per read
or per anchor in mapping or surjection is `important`, not `minor`: it is
multiplied by every read of every run.

# Concurrency lens

Ask only: what happens when two threads run this at once?

vg uses OpenMP. Parallel work comes from `#pragma omp parallel` blocks,
`#pragma omp task`, `for_each_handle(..., true)`, `for_each_edge(..., true)`,
`for_each_path_handle` variants with a parallel flag, and the
`vg::io::for_each_parallel` family that reads records on many threads.

## Shared writes

- A shared container, counter or stream written from a parallel loop without
  `#pragma omp critical`, an `atomic`, or a per-thread buffer indexed by
  `omp_get_thread_num()`.
- A per-thread buffer sized by `omp_get_max_threads()` outside the parallel
  region but indexed inside a nested region, where the count can differ.
- Output emitted from many threads into one `vg::io::ProtobufEmitter` or one
  stream: check it goes through a lock or an emitter meant for it.

## Shared reads

- Reading a graph from many threads is safe only if nothing mutates it. Check
  for lazy caches inside "const" methods (a `mutable` member filled on first
  use).
- A static local or global initialized on first use from several threads.

## Ordering

- Output whose order depends on thread scheduling, where tests or users expect
  a stable order. Say whether the order matters.
- A `critical` section whose name is shared with unrelated code, so the two
  serialize each other for no reason (a performance finding, not a race).

## Deadlock

- Nested `critical` sections with the same name, locks taken in different
  orders, or a task that waits on a task that cannot start because all threads
  are waiting.

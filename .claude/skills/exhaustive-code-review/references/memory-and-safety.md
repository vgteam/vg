# Memory and safety lens

Ask only: can this code read or write memory it should not, or be made to by
its input?

## Lifetimes

- A reference or pointer into a `std::vector`, `std::string` or hash map that
  outlives an insertion into that container.
- A `string_view`, `const char*` or reference bound to a temporary.
- A lambda that captures by reference and runs later, on another thread or
  after the enclosing function returns.
- A `handle_t`, `path_handle_t` or `step_handle_t` used after the graph that
  made it was changed in a way that invalidates it (for example after
  `destroy_handle`, `apply_ordering` or a compaction). Each graph implementation
  documents its own rules; check the one in use.
- Ownership: `unique_ptr` for owned objects, no `new` without a matching owner.
  A raw pointer member should say who owns it.

## Bounds

- Indexing with `[]` by a value that came from a file, a command-line option,
  or arithmetic that can underflow.
- `substr`, `memcpy` and offset arithmetic on node sequences, with offsets in
  the wrong orientation.
- Reading past the end of a record when a file line has fewer fields than
  expected.

## Input from outside

vg reads GFA, GAM, GAF, VCF, FASTA and FASTQ files and index files that users
make or download. Treat each as possibly malformed.

- A parser that assumes a field exists, is numeric, or fits in its type.
- A count or length read from a file and used to size an allocation.
- A file name or path name passed to a shell (`system`, `popen`) without
  quoting.

## Undefined behavior

- Shifts by the type's width or more, signed overflow, reading an inactive
  union member, `reinterpret_cast` between unrelated types, misaligned access.
- A function that should return a value but can fall off the end (vg builds
  with `-Werror=return-type`, so this usually fails to compile; check
  `[[noreturn]]` helpers).

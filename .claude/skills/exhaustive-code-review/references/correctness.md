# Correctness lens

Ask only: does this code compute the right thing for every input it can get?

## Intent

- Does the change do what its commit messages or PR description say? Name any
  case the description promises that the code does not handle.
- Do the doc comments match the code? A stated precondition nothing enforces is
  a finding (see the `vg-code-review` skill).

## Edge cases to try in your head

- Empty input: an empty graph, a path with no steps, a node with an empty
  sequence, an empty read, a file with only a header.
- One of something: a single-node graph, a trivial chain, a snarl with only the
  edge between its bounds.
- Orientation: each handle forward and reverse; a path that visits a node in
  reverse; a child reversed in its parent; a traversal read from the end.
- Cycles and self-loops: an edge from a node to itself, a looping chain, a
  circular path, a reversing edge (a node's right side to its own right side).
- Boundaries: offset 0 and offset equal to the node length; the first and last
  step of a path; the largest `size_t` as "no distance".
- IDs: node IDs that do not start at 1, have gaps, or exceed 32 bits.

## Arithmetic and types

- Signed and unsigned mixing, especially `size_t` minus something that can be
  larger, and comparisons of `int` with `size_t`.
- Off-by-one in half-open ranges; vg path subranges are 0-based and
  end-exclusive, VCF and GFF are 1-based.
- Integer overflow in sums of lengths or scores.

## Error handling

- Is every failure reported, with the subcommand's logger, and does it stop the
  program or get handled? A function that returns a sentinel needs every caller
  to check it.
- Is the error message true and does it name the bad input?

## State

- Is every variable set on every path before it is read?
- Does a loop that removes or adds while iterating keep its iterator valid?
- Does code that edits a graph keep its invariants: each edge has both ends,
  paths still spell valid walks, IDs are unique?

Report each finding as the input that triggers it and the wrong result.

# VG Project Notes

## Building
- New `.cpp` files auto-discovered
- Build with `make -j8` or `make obj/whatever.o` to build just one .o.
- Changing a dependency's source under `deps/` is slow: the Makefile cleans and rebuilds the dependency, copies its headers into `include/` again, and every object depends on the libraries, so all of vg recompiles. Incremental builds across projects are not possible yet. TODO: Make them possible, for example by copying headers with `rsync` so unchanged ones keep their times, and by making objects depend on the headers they include rather than on the libraries.
  - Until then, to try a change to a library such as libbdsg: run `make` in the dependency's own directory, copy its `.a` into `lib/`, delete `bin/vg`, and run `make -o lib/libbdsg.a` (naming the library you changed), which relinks without remaking the library or recompiling vg.
- You may be getting errors from `clangd`. If these errors seem spurious, stop and demand a `clangd` that works properly.

## Testing

### Running Bash-TAP Tests
Use `prove -v` (not `bash`) to execute Bash-TAP tests. This provides proper test harness output and better error reporting.

**Important**: Run `prove` from the `test/` directory:
```bash
cd test
prove -v t/26_deconstruct.t
```

### Running Unit Tests
To run all unit tests:
```bash
./bin/vg test
```
- `./bin/vg test "[tag]"` runs tests matching a tag

#### Writing Unit Tests
- Framework: Catch v2 (header-only)
- Include: `#include "catch.hpp"` (in `src/unittest/catch.hpp`)
- Macros: `TEST_CASE("name", "[tags]")`, `SECTION("name")`, `REQUIRE(cond)`
- Namespace: `vg::unittest`
- Directory: `src/unittest/`

### Running All Tests
```bash
make test
```

## Writing Code

### HandleGraph API
The interfaces in libhandlegraph model a bidirected sequence graph (where nodes have DNA sequences and edges can connect to either the start or end of each involved node).

#### Core types
- `handle_t` - opaque 64-bit value
- `nid_t` - node ID type
- `edge_t` = `pair<handle_t, handle_t>`

#### Key HandleGraph methods
- `get_handle(nid_t, bool is_reverse=false)` → `handle_t`
- `get_id(handle_t)` → `nid_t`
- `get_is_reverse(handle_t)` → `bool`
- `flip(handle_t)` → `handle_t` (toggle orientation)
- `get_sequence(handle_t)` → `string` (in handle's orientation)
- `follow_edges(handle_t, bool go_left, iteratee)` - iterate neighbors
- `for_each_handle(iteratee, bool parallel=false)` - iterate all nodes
- `for_each_edge(iteratee, bool parallel=false)` - iterate all edges
- `has_edge(handle_t left, handle_t right)` → `bool`

#### MutableHandleGraph additions
- `create_handle(string seq)` / `create_handle(string seq, nid_t id)` → `handle_t`
- `create_edge(handle_t left, handle_t right)`
- `destroy_handle(handle_t)` / `destroy_edge(handle_t, handle_t)`

#### HandleGraph algorithms
- Things like `topological_sort.hpp` and copy_graph.hpp` are in `deps/libhandlegraph/src/include/handlegraph/algorithms`.

#### bdsg::HashGraph
- Header: `deps/libbdsg/bdsg/include/bdsg/hash_graph.hpp`
- Implements MutablePathMutableHandleGraph
- Go-to handlegraph implementation to use
- In libbdsg

### Utilities
- `reverse_complement(string)` → `string` in src/utility.hpp

## Terminology
Use these words, which the code already uses, in comments and names. A synonym makes a reader look for a second idea.
- A *handle* (`handle_t`) is a node in one *orientation*: forward or reverse. A handle reads the node's sequence in its orientation, and leaves the node through its right side. Say a handle is *oriented* forward or reverse; `flip()` gives the other orientation.
- An *edge* joins a side of one node to a side of another. `follow_edges()` gives the handles an edge leads to.
- A *Visit* is one step of a walk through the graph: a node, or a child snarl, in an orientation. A *SnarlTraversal* or *path* *visits* nodes and snarls; it does not "pass" or "reach" them.
- A *snarl* is a part of the graph between two *bounds*. Its `start()` Visit is oriented into the snarl and its `end()` Visit out of it. Snarls *contain* chains.
- A *chain* is a run of snarls that share bounds, each snarl's end being the next one's start. The nodes in a chain are oriented *forward along the chain* or *in reverse along the chain*.
- In a `SnarlDecomposition`, a *net handle* (`net_handle_t`) is one *traversal* of a snarl, chain or node: from its start, end or a *tip* to its start, end or a tip. A *sentinel* is a snarl's bound, facing into or out of the snarl. A traversal is *realizable* when a walk in the graph can make it.
- A *dead end* is a side of a node with no edges. A *tip* is a dead end inside a snarl or chain.

## Writing Code and Comments
- When writing or reviewing doc comments, follow `.claude/skills/doc-comments/SKILL.md`.
- When writing or reviewing comments, identifiers or test names, follow `.claude/skills/concept-phrasing/SKILL.md`.
- When adding, splitting or reviewing classes and structs, follow `.claude/skills/class-identity/SKILL.md`.
- When writing, restructuring or reviewing unit tests, follow `.claude/skills/test-design/SKILL.md`.
- When adding, changing or reviewing a skill in `.claude/skills/`, follow `.claude/skills/writing-skills/SKILL.md`.

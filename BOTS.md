# VG Project Notes

## Building
- New `.cpp` files auto-discovered
- Build with `make -j8` or `make obj/whatever.o` to build just one .o.
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


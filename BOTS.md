# VG Project Notes

## Building
- New `.cpp` files auto-discovered
- Build with `make -j8` or `make obj/whatever.o` to build just one .o.
- Changing a dependency's source under `deps/` is slow: the Makefile cleans and rebuilds the dependency, copies its headers into `include/` again, and every object depends on the libraries, so all of vg recompiles. Incremental builds across projects are not possible yet. TODO: Make them possible, for example by copying headers with `rsync` so unchanged ones keep their times, and by making objects depend on the headers they include rather than on the libraries.
  - Until then, to try a change to a library such as libbdsg: run `make` in the dependency's own directory, copy its `.a` into `lib/`, delete `bin/vg`, and run `make -o lib/libbdsg.a` (naming the library you changed), which relinks without remaking the library or recompiling vg.
- You may be getting errors from `clangd`. If these errors seem spurious, stop and demand a `clangd` that works properly.

## Waiting for Background Work
- Never wait for a process by polling `pgrep -f` or `pkill -f` with a pattern. The shell running the check has the pattern in its own command line, so the check finds itself and never ends. Wait on the job itself instead: run it in the background and wait to be told it finished, `wait` on its PID, or poll for a file the job writes as its last step.

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
- Unit tests exercise library code, so test names and comments name the library setting (a setter or field) or describe the behaviour, never a command-line option.

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

### Snarls
- Find snarls with `IntegratedSnarlFinder`. Only the older `CactusSnarlFinder` can produce unary snarls, and it is to be removed; new code may require a `SnarlManager` with no unary snarls, as a stated precondition.
- Give the orientation of a child relative to its chain as an `is_reverse` flag, as `is_reversed_in_parent()` does.

## Terminology
Use these words, which the code already uses, in comments and names, without defining them again. A synonym makes a reader look for a second idea.

A term here names its idea in sentences about that idea. The same word can mean something else in another context: its ordinary English sense, or a meaning that a module, class or function defines in its own comments. Leave a word alone where it is clear in its context, even if this list gives it another meaning.
- A *handle* is a node in one *orientation*: forward or reverse. A handle reads the node's sequence in its orientation, and leaves the node through its right side. Say a handle is *oriented* forward or reverse. Flipping a handle gives the other orientation.
- An *edge* joins a side of one node to a side of another.
- A *visit* is one step of a walk through the graph: a node, or a child snarl, in an orientation. A traversal of a snarl, or a *path*, *visits* nodes and snarls; they does not "pass" or "reach" them. The capital Visit suggests the Protobuf `Visit` type, which can only be to nodes; the lower-case "visit" is preferred for the general concept. Generally paths are thought of as only visiting nodes, but in some contexts we might think of a path as visiting snarls or chains.
- A *snarl* is a part of the graph between two *bounds*. Its start is oriented into the snarl and its end out of it. Snarls *contain* chains.
- A *chain* is a run of snarls that share bounds, each snarl's end being the next one's start. Each node or snarl in a chain is forward or *reversed in parent* relative to the chain.
- In a snarl decomposition, a *net handle* is one *traversal* of a snarl, chain or node: from its start, end or a *tip* to its start, end or a tip. A *sentinel* is a snarl's bound, facing into or out of the snarl. A traversal is *realizable* when a walk in the graph can make it.
- A *dead end* is a side of a node with no edges. A *tip* is a dead end inside a snarl or chain. The outer sides of the boundary nodes of a snarl or chain are outside it, so a dead end there is not one of its tips, though it can be a tip of the parent.
- *variation graph*: A graph whose nodes carry DNA sequence and whose paths spell genomes, so that variants appear as alternative routes between shared sequence.
- *bidirected graph*: A graph in which each edge joins a side of one node to a side of another, so that a walk can read each node forward or reverse.
- *node*: A piece of sequence in the graph, with an ID and two sides.
- *node ID*: The integer that names a node.
- *side*: One of the two ends of a node, left (its start) or right (its end), where edges attach.
- *graph position*: A place in the graph: a node ID, an orientation and an offset into the node in that orientation. Call it a position for short where nothing else has one; other things, such as a node in a chain, have positions of their own.
- *connected component*: A set of nodes that edges connect to each other and to nothing outside the set.
- *path*: A named walk through the graph that the graph stores, such as a reference chromosome or a haplotype.
- *step*: One visit of a path to a node, in an orientation.
- *circular path*: A path whose last step leads back to its first, as for a circular chromosome.
- *path sense*: The kind of a path: reference, for a linear reference such as GRCh38; haplotype, for one haplotype of a sample; or generic, for anything else, such as a gene or a read.
- *sample*: The individual a haplotype path comes from, or the assembly a reference path belongs to.
- *haplotype*: One copy of a genome region that an individual carries; for paths, the haplotype number says which copy (0 or 1 in a diploid).
- *locus*: The chromosome or contig that a path belongs to, or for a generic path the name of the thing it represents.
- *phase block*: A stretch of a contig over which a haplotype is phased, numbered so that a haplotype cut where phasing is unknown becomes several paths.
- *subrange*: The 0-based, end-exclusive coordinates that a reference or generic path covers within a longer path that the graph does not hold in full.
- *PanSN*: The naming scheme `sample#haplotype#contig` for paths in a pangenome.
- *bubble*: The older, narrower term in the literature for something like a snarl: alternative routes between two nodes.
- *ultrabubble*: A snarl in which every walk from the start reaches the end without cycles or turning around.
- *boundary node*: One of the two nodes whose sides delimit a snarl or chain; a snarl's start is oriented into it and its end out of it.
- *trivial chain*: A chain that holds just one node, which is both of its bounds.
- *trivial snarl*: A snarl with nothing inside it but an edge between its boundary nodes, implied by two consecutive nodes of a chain and often left out of the snarl tree.
- *unary snarl*: A snarl whose start and end are the same node in opposite orientations, so it hangs off one side of that node. Eventually it was figured out that it is always possible to create a decomposition of a graph that did not use any unary snarls, and this is now the preferred approach.
- *looping chain*: A chain whose first and last node are the same node, as on a circular chromosome.
- *start-end connected*: Said of a snarl or chain that some walk inside it goes from its start to its end; the snarl decomposition can contain snarls and chains that are not.
- *nesting*: The containment of snarls in chains and chains in snarls; a snarl contains a chain when it contains all the chain's nodes and snarls.
- *snarl tree*: The tree of nested snarls and chains, with every node, snarl and chain of the graph as a vertex, alternating snarl levels and chain levels; also called the snarl decomposition.
- *child*: A snarl or node that is a component of a chain, or a chain that a snarl contains directly; the snarl or chain is its parent.
- *top-level chain*: A chain with no snarl around it, at the top of the snarl tree.
- *root*: The top of a SnarlDecomposition. A snarl with undefined bounds that holds the top-level chains.
- *netgraph*: The view of a snarl in which each child chain is seen as a single element, since a chain, like a node, joins the rest of the graph at two sides.
- *net graph element*: A snarl, chain or node in a snarl decomposition. A thing that net handles are traversals of. An *element* for short.
- *site*: A snarl. A place where haplotypes differ. Used in code that finds walks through it or calls variants in it.
- *endpoint*: One end of a traversal: the start, the end, or a tip.
- *canonical traversal*: The one realizable traversal that every net handle to the same snarl, chain or node maps to, for comparing net handles.
- *distance index*: libbdsg's index of the snarl tree that also answers minimum-distance queries between positions.
- *minimum distance*: The length of the shortest walk between two positions, counting one of them; the largest size_t value when no walk exists.
- *depth*: The level of a chain in the snarl tree, with the root at 0 and top-level chains at 1; nodes and snarls take the depth of their parent chain.
- *reversed in (its) parent (snarl/chain)*: Said of a node or snarl that a start-to-end traversal of its chain crosses end to start, or of a chain that is similarly situated relative to its parent snarl (when its parent snarl's structure suggests a coherent orientation for its contents).


### Module terms
Some words mean something narrower, or something else, inside one module. In that module's code and comments, use the meaning given here.

#### Snarl distance index (`deps/libbdsg/bdsg/include/bdsg/snarl_distance_index.hpp`)
- *snarl*: as above, but without its boundary nodes, which the snarl holds only as sentinels. A snarl may hold no nodes, if edges other than the one between its boundary nodes run through it.
- *sentinel*: the inner side of a boundary node of a snarl, as a net handle; not the node itself.
- *root snarl*: a snarl of the root with no start or end node, grouping root-level chains that have connectivity between them.
- *record*: the stored entry for one snarl, chain or node; its record type says what it stores.
- *simple snarl*: a bubble whose children are nodes that each reach only the boundary nodes, each side of a child reaching a different boundary node, with at most an edge between the boundary nodes besides.
- *regular snarl*: a snarl like a simple snarl, except that its children may be chains that allow no reversals; simple and trivial snarls are regular too.
- *oversized snarl*: a snarl with more children than the size limit, which keeps hub labeling data instead of a table of distances.
- *multicomponent chain*: a chain holding snarls that are not start-end connected, split into components between them, with each node tagged with its component.
- *length*: of a node, its sequence length; of a chain, the minimum length including its boundary nodes, counting a looping chain's shared node once; of a snarl, the minimum length without its boundary nodes.
- *rank in parent*: the index of a child within its parent, used to find distances between children of a snarl.
- *connected component number*: the number of a connected component of the graph; two root-level chains connected in the root share one.
- *externally connected*: when two things in a snarl or chain are connected by a walk that goes out into the containing parent chain/snarl/root and then comes back.
- *temporary index*: the mutable description of a snarl tree from which the index is built.
- *prefix sum*: When computing distances, for a node in a chain, the distance from the start of the chain to the start of the node.

## Machine-Written Files
- New source files, unit tests and write-ups that you write go in `artifacts/`, `artifacts/unittest/` and `artifacts/doc/`, not in `src/` or `doc/`. A file moves to `src/` or `doc/` only when a person has read it and cleared it for reading. See `artifacts/README.md`.
- A file in `artifacts/` other than a `README.md` must not share a name with a file in `src/` or `doc/`, since objects, unit test binaries and includes are found by name. Check this yourself when you add a file.

## Writing Classes and Functions
- When adding, splitting or reviewing classes and structs, use the `class-identity` skill.
- Before finishing a task that adds or redesigns a class or adds a function, lint them with the `lint-new-code` skill. Keep each class and function it makes you split coherent by the `class-identity` skill.

## Writing Comments
- Before finishing a task in which you wrote or changed code comments, use the `intelligible-comments` skill to rewrite them so that a programmer new to the code understands each one on one reading.
- Name a setting by what it means, or by its name in the library code, never by a command-line option. This holds in library comments and in unit tests: their comments, `TEST_CASE` and `SECTION` names, and messages. Only `src/subcommand/` names options. The `lint_prose.py` script in the `intelligible-comments` skill checks this, but finds only long options such as `--min-mapq`, so check for short ones such as `-L` yourself.
- Do not write "gate", "gated" or "gating" for the condition under which something happens. Say what must be true for it to happen ("the depth discount applies where the alleles differ in length by 50 bp or more"). The same script checks this.
- Say what a thing is. Do not define it by likeness or by reference to something else ("oriented like the start and end of a Snarl", "in the order `find_traversals()` uses"). If the next sentence explains the thing, delete the indirect one.

## Documentation
- When writing or reviewing doc comments or design docs, follow `.claude/skills/explanatory-comments/SKILL.md`. (Yes, there are two different skills for working on comments. This one is more about what ought to be explained, and how to organize ideas, rather than how to write them out in words.) 

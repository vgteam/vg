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

### Snarls
- Find snarls with `IntegratedSnarlFinder`. Only the older `CactusSnarlFinder` can produce unary snarls, and it is to be removed; new code may require a `SnarlManager` with no unary snarls, as a stated precondition.
- Give the orientation of a child relative to its chain as an `is_reverse` flag, as `is_reversed_in_parent()` does.

## Terminology
Use these words, which the code already uses, in comments and names, without defining them again. A synonym makes a reader look for a second idea.
- A *handle* (`handle_t`) is a node in one *orientation*: forward or reverse. A handle reads the node's sequence in its orientation, and leaves the node through its right side. Say a handle is *oriented* forward or reverse; `flip()` gives the other orientation.
- An *edge* joins a side of one node to a side of another. `follow_edges()` gives the handles an edge leads to.
- A *Visit* is one step of a walk through the graph: a node, or a child snarl, in an orientation. A *SnarlTraversal* or *path* *visits* nodes and snarls; it does not "pass" or "reach" them.
- A *snarl* is a part of the graph between two *bounds*. Its `start()` Visit is oriented into the snarl and its `end()` Visit out of it. Snarls *contain* chains.
- A *chain* is a run of snarls that share bounds, each snarl's end being the next one's start. Each node or snarl in a chain is forward or *reversed in parent* relative to the chain (`is_reversed_in_parent()`).
- In a `SnarlDecomposition`, a *net handle* (`net_handle_t`) is one *traversal* of a snarl, chain or node: from its start, end or a *tip* to its start, end or a tip. A *sentinel* is a snarl's bound, facing into or out of the snarl. A traversal is *realizable* when a walk in the graph can make it.
- A *dead end* is a side of a node with no edges. A *tip* is a dead end inside a snarl or chain.
- *variation graph* (`HandleGraph`): A graph whose nodes carry DNA sequence and whose paths spell genomes, so that variants appear as alternative routes between shared sequence.
- *bidirected graph* (`HandleGraph`): A graph in which each edge joins a side of one node to a side of another, so that a walk can read each node forward or reverse.
- *node* (`nid_t`, `has_node()`, `get_length()`, `get_sequence()`): A piece of sequence in the graph, with an ID and two sides.
- *node ID* (`nid_t`, `get_id()`, `min_node_id()`, `max_node_id()`): The integer that names a node.
- *side* (`go_left`): One of the two ends of a node, left (its start) or right (its end), where edges attach.
- *position* (`pos_t`, `offset_t`): A place in the graph: a node ID, an orientation and an offset into the node in that orientation.
- *connected component* (`get_connected_component_number()`): A set of nodes that edges connect to each other and to nothing outside the set.
- *path* (`path_handle_t`, `PathHandleGraph`, `get_path_name()`): A named walk through the graph that the graph stores, such as a reference chromosome or a haplotype.
- *step* (`step_handle_t`, `get_next_step()`, `get_previous_step()`, `path_begin()`, `path_end()`): One visit of a path to a node, in an orientation.
- *circular path* (`get_is_circular()`): A path whose last step leads back to its first, as for a circular chromosome.
- *path sense* (`PathSense`, `PathSense::REFERENCE`, `PathSense::HAPLOTYPE`, `PathSense::GENERIC`, `get_sense()`): The kind of a path: reference, for a linear reference such as GRCh38; haplotype, for one haplotype of a sample; or generic, for anything else, such as a gene or a read.
- *sample* (`get_sample_name()`, `PathMetadata::NO_SAMPLE_NAME`): The individual a haplotype path comes from, or the assembly a reference path belongs to.
- *haplotype* (`get_haplotype()`, `PathMetadata::NO_HAPLOTYPE`): One copy of a genome region that an individual carries; for paths, the haplotype number says which copy (0 or 1 in a diploid).
- *locus* (`get_locus_name()`): The chromosome or contig that a path belongs to, or for a generic path the name of the thing it represents.
- *phase block* (`get_phase_block()`, `PathMetadata::NO_PHASE_BLOCK`): A stretch of a contig over which a haplotype is phased, numbered so that a haplotype cut where phasing is unknown becomes several paths.
- *subrange* (`get_subrange()`, `subrange_t`, `PathMetadata::NO_SUBRANGE`): The 0-based, end-exclusive coordinates that a reference or generic path covers within a longer path that the graph does not hold in full.
- *PanSN* (`PathMetadata::create_path_name()`, `parse_path_name()`): The naming scheme `sample#haplotype#contig` for paths in a pangenome.
- *bubble*: The older, narrower term in the literature for a snarl: alternative routes between two nodes.
- *ultrabubble* (`ULTRABUBBLE`): A snarl in which every walk from the start reaches the end without cycles or turning around.
- *boundary node* (`Snarl::start()`, `Snarl::end()`, `get_bound()`, `get_start_bound()`, `get_end_bound()`): One of the two nodes whose sides delimit a snarl or chain; a snarl's start is oriented into it and its end out of it.
- *bound* (`get_bound()`, `get_start_bound()`, `get_end_bound()`): The start or the end of a snarl or chain, as a boundary node traversal or, in a SnarlDecomposition, as a sentinel.
- *trivial chain* (`is_trivial_chain()`, `NODE` (record type)): A chain that holds just one node, which is both of its bounds.
- *trivial snarl* (`TRIVIAL_SNARL`, `is_trivial_snarl()`): A snarl with nothing inside it but an edge between its boundary nodes, implied by two consecutive nodes of a chain and usually left out of the snarl tree.
- *unary snarl* (`UNARY`, `add_unary_child()`): A snarl whose start and end are the same node in opposite orientations, so it hangs off one side of that node.
- *looping chain* (`is_looping_chain()`): A chain whose first and last node are the same node, as on a circular chromosome.
- *start-end connected* (`is_start_end_connected()`, `start_end_reachable`): Said of a snarl or chain that some walk inside it goes from its start to its end; the snarl decomposition can contain snarls and chains that are not.
- *nesting* (`get_parent()`, `for_each_child()`): The containment of snarls in chains and chains in snarls; a snarl contains a chain when it contains all the chain's nodes and snarls.
- *snarl tree* (`SnarlDecomposition`, `SnarlManager`, `SnarlDistanceIndex`): The tree of nested snarls and chains, with every node, snarl and chain of the graph as a vertex, alternating snarl levels and chain levels; also called the snarl decomposition.
- *child* (`for_each_child()`, `get_parent()`): A snarl or node that is a component of a chain, or a chain that a snarl contains directly; the snarl or chain is its parent.
- *top-level chain*: A chain with no snarl around it, at the top of the snarl tree.
- *root* (`get_root()`, `is_root()`): The top of a SnarlDecomposition: a snarl with undefined bounds that holds the top-level chains.
- *netgraph* (`directed_acyclic_net_graph`): The view of a snarl in which each child chain is replaced by a single node, since a chain, like a node, joins the rest of the graph at two sides.
- *SnarlTraversal* (`SnarlTraversal`): A walk through a snarl as a list of Visits, from its start node to its end node, with no edits.
- *SnarlManager* (`SnarlManager`): vg's in-memory store of the snarls and chains found in a graph.
- *endpoint* (`endpoint_t`, `START`, `END`, `TIP`, `starts_at()`, `ends_at()`): One end of a traversal: the start, the end, or a tip.
- *canonical traversal* (`canonical()`): The one realizable traversal that every net handle to the same snarl, chain or node maps to, for comparing net handles.
- *distance index* (`SnarlDistanceIndex`, `minimum_distance()`): libbdsg's index of the snarl tree that also answers minimum-distance queries between positions.
- *minimum distance* (`minimum_distance()`, `distance_in_parent()`): The length of the shortest walk between two positions, counting one of them; the largest size_t value when no walk exists.
- *depth* (`get_depth()`, `get_max_tree_depth()`): The level of a chain in the snarl tree, with the root at 0 and top-level chains at 1; nodes and snarls take the depth of their parent chain.
- *reversed in parent* (`is_reversed_in_parent()`): Said of a node or snarl that a start-to-end traversal of its chain crosses end to start.
- *prefix sum* (`get_prefix_sum_value()`): For a node in a chain, the distance from the start of the chain to the start of the node.

### Module terms
Some words mean something narrower, or something else, inside one module. In that module's code and comments, use the meaning given here.

#### Snarl distance index (`deps/libbdsg/bdsg/include/bdsg/snarl_distance_index.hpp`)
- *snarl*: as above, but without its boundary nodes, which the snarl holds only as sentinels. A snarl may hold no nodes, if edges other than the one between its boundary nodes run through it.
- *sentinel*: the inner side of a boundary node of a snarl, as a net handle; not the node itself (`SENTINEL_HANDLE`, `get_node_from_sentinel()`).
- *root snarl*: a snarl of the root with no start or end node, grouping root-level chains that have connectivity between them (`ROOT_SNARL`, `is_root_snarl()`).
- *record*: the stored entry for one snarl, chain or node; its record type says what it stores (`record_t`, `SnarlTreeRecord`).
- *simple snarl*: a bubble whose children are nodes that each reach only the boundary nodes, each side of a child reaching a different boundary node, with at most an edge between the boundary nodes besides (`SIMPLE_SNARL`, `is_simple_snarl()`).
- *regular snarl*: a snarl like a simple snarl, except that its children may be chains that allow no reversals; simple and trivial snarls are regular too (`REGULAR_SNARL`, `is_regular_snarl()`).
- *oversized snarl*: a snarl with more children than the size limit, which keeps hub labeling data instead of a table of distances (`OVERSIZED_SNARL`, `is_oversized_snarl()`).
- *multicomponent chain*: a chain holding snarls that are not start-end connected, split into components between them, with each node tagged with its component (`MULTICOMPONENT_CHAIN`).
- *length*: of a node, its sequence length; of a chain, the minimum length including its boundary nodes, counting a looping chain's shared node once; of a snarl, the minimum length without its boundary nodes (`minimum_length()`).
- *rank in parent*: the index of a child within its parent, used to find distances between children of a snarl (`get_rank_in_parent()`).
- *connected component number*: the number of a connected component of the graph; two root-level chains connected in the root share one (`get_connected_component_number()`).
- *externally connected*: said of the endpoints of a child of the root that a walk outside the child joins (`is_externally_start_end_connected()`).
- *temporary index*: the description of one connected component's snarl tree from which `get_snarl_tree_records()` builds the index (`TemporaryDistanceIndex`).

## Writing Classes and Functions
- When adding, splitting or reviewing classes and structs, use the `class-identity` skill.
- Before finishing a task that adds or redesigns a class or adds a function, lint them with the `lint-new-code` skill. Keep each class and function it makes you split coherent by the `class-identity` skill.

## Writing Comments
- Before finishing a task in which you wrote or changed code comments, use the `intelligible-comments` skill to rewrite them so that a programmer new to the code understands each one on one reading.

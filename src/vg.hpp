#ifndef VG_VG_H

#define VG_VG_H

#include <vector>
#include <set>
#include <string>
#include <deque>
#include <list>
#include <array>
#include <omp.h>
#include <unistd.h>
#include <limits.h>
#include <algorithm>
#include <random>

#include "gssw.h"
#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "gssw_aligner.hpp"
#include "ssw_aligner.hpp"
#include "region.hpp"
#include "path.hpp"
#include "utility.hpp"
#include "alignment.hpp"

#include "vg.pb.h"
#include "hash_map.hpp"

#include "progressive.hpp"
#include "lru_cache.h"

#include "Variant.h"
#include "Fasta.h"

#include "swap_remove.hpp"

#include "pictographs.hpp"
#include "colors.hpp"

#include "types.hpp"
#include "gfakluge.hpp"

#include "globalDefs.hpp"
#include "Graph.hpp"
#include "helperDefs.hpp"

#include "bubbles.hpp"

#include "nodetraversal.hpp"
#include "nodeside.hpp"

// uncomment to enable verbose debugging to stderr
//#define debug

namespace vg {

/**
 * We create a struct that represents each kmer record we want to send to gcsa2
 */
struct KmerPosition {
    string kmer;
    string pos;
    set<char> prev_chars;
    set<char> next_chars;
    set<string> next_positions;
};

class Aligner; // forward declarations
class QualAdjAligner;

}

namespace vg {

/**
 * Represents a variation graph. Graphs consist of nodes, connected by edges.
 * Graphs are bidirected and may be cyclic. Nodes carry forward-oriented
 * sequences. Edges are directed, with a "from" and to" node, and are generally
 * used to connect the end of the "from" node to the start of the "to" node.
 * However, edges can connect to either the start or end of either node.
 *
 */
class VG : public Progressive {

public:

    /// Protobuf-based representation.
    // NB: we can't subclass this safely, so it's best as a member
    Graph graph;

    /// Manages paths of the graph.
    /// Initialized by setting paths._paths = graph.paths.
    Paths paths;

    /// Name of the graph.
    string name;

    /// Current id for Node to be added next.
    id_t current_id;
    // todo
    //id_t min_id;
    //id_t max_id;

    /// `Node`s by id.
    hash_map<id_t, Node*> node_by_id;

    /// `Edge`s by sides of `Node`s they connect.
    /// Since duplicate edges are not permitted, two edges cannot connect the same pair of node sides.
    /// Each edge is indexed here with the smaller NodeSide first. The actual node order is recorded in the Edge object.
    pair_hash_map<pair<NodeSide, NodeSide>, Edge*> edge_by_sides;

    /// nodes by position in nodes repeated field.
    /// this is critical to allow fast deletion of nodes
    hash_map<Node*, int> node_index;

    // edges by position in edges repeated field.
    // same as for nodes, this allows fast deletion.
    hash_map<Edge*, int> edge_index;

    // edges indexed by nodes they connect.
    
    /// Stores the destinations and backward flags for edges attached to the starts of nodes (whether that node is "from" or "to").
    hash_map<id_t, vector<pair<id_t, bool>>> edges_on_start;
    /// Stores the destinations and backward flags for edges attached to the ends of nodes (whether that node is "from" or "to").
    hash_map<id_t, vector<pair<id_t, bool>>> edges_on_end;

    /// Set the edge indexes through this function. Picks up the sides being
    /// connected by the edge automatically, and silently drops the edge if they
    /// are already connected.
    void set_edge(Edge*);
    void print_edges(void);

    // access the edge indexes through these functions
    
    /// Get nodes and backward flags following edges that attach to this node's start.
    vector<pair<id_t, bool>>& edges_start(Node* node);
    /// Get nodes and backward flags following edges that attach to this node's start.
    vector<pair<id_t, bool>>& edges_start(id_t id);
    /// Get nodes and backward flags following edges that attach to this node's end.
    vector<pair<id_t, bool>>& edges_end(Node* node);
    /// Get nodes and backward flags following edges that attach to this node's end.
    vector<pair<id_t, bool>>& edges_end(id_t id);
    
    // properties of the graph
    size_t size(void); ///< Number of nodes
    size_t length(void); ///< Total sequence length

    // Clear everything
    //void clear(void);

    // constructors

    /// Default constructor.
    VG(void);

    /// Construct from protobufs.
    VG(istream& in, bool showp = false);

    /// Construct from an arbitrary source of Graph protobuf messages (which
    /// populates the given Graph and returns a flag for whether it's valid).
    VG(function<bool(Graph&)>& get_next_graph, bool showp = false);

    /// Construct from sets of nodes and edges. For example, from a subgraph of
    /// another graph.
    VG(set<Node*>& nodes, set<Edge*>& edges);

    /// Takes in a VCF file
    /// and returns a map [node] = vcflib::variant.
    /// Unfortunately this is specific to a given graph
    /// and VCF.
    ///
    /// It will need to throw warnings if the node or variant
    /// is not in the graph.
    ///
    /// This is useful for VCF masking:
    /// 
    ///     if map.find(node) then mask variant
    ///
    /// It's also useful for calling known variants
    /// 
    ///     for m in alignment.mappings:
    ///        node = m.Pos.nodeID
    ///        if node in node_to_vcf:
    ///            return (alignment supports variant)
    ///
    /// It would be nice if this also supported edges (e.g.
    /// for inversions/transversions/breakpoints?).
    // TODO: map<edge_id, variant> or map<pair<NodeID, NodeID>, variant>
    map<id_t, vcflib::Variant> get_node_id_to_variant(vcflib::VariantCallFile vfile);
                       
                       
    /// Chop up the nodes.
    void dice_nodes(int max_node_size);
    /// Does the reverse --- combines nodes by removing edges where doing so has no effect on the graph labels.
    void unchop(void);
    /// Get the set of components that could be merged into single nodes without
    /// changing the path space of the graph. Emits oriented traversals of
    /// nodes, in the order and orientation in which they are to be merged.
    set<list<NodeTraversal>> simple_components(int min_size = 1);
    /// Get the simple components of multiple nodes.
    set<list<NodeTraversal>> simple_multinode_components(void);
    /// Get the strongly connected components of the graph.
    set<set<id_t> > strongly_connected_components(void);
    /// Get only multi-node strongly connected components.
    set<set<id_t> > multinode_strongly_connected_components(void);
    /// Returns true if the graph does not contain cycles.
    bool is_acyclic(void);
    /// Remove all elements which are not in a strongly connected component.
    void keep_multinode_strongly_connected_components(void);
    /// Does the specified node have any self-loops?
    bool is_self_looping(Node* node);
    /// Get simple cycles following Johnson's elementary cycles algorithm.
    set<list<NodeTraversal> > elementary_cycles(void);
    /// Concatenates the nodes into a new node with the same external linkage as
    /// the provided component. After calling this, paths will be invalid until
    /// Paths::compact_ranks() is called.
    Node* concat_nodes(const list<NodeTraversal>& nodes);
    /// Merge the nodes into a single node, preserving external linkages.
    /// Use the orientation of the first node as the basis.
    Node* merge_nodes(const list<Node*>& nodes);
    /// Use unchop and sibling merging to simplify the graph into a normalized form.
    void normalize(int max_iter = 1);
    /// Remove redundant overlaps.
    void bluntify(void);
    /// Turn the graph into a dag by copying strongly connected components expand_scc_steps times
    /// and translating the edges in the component to flow through the copies in one direction.
    VG dagify(uint32_t expand_scc_steps,
              map<id_t, pair<id_t, bool> >& node_translation,
              size_t target_min_walk_length = 0,
              size_t component_length_max = 0);
    /// Generate a new graph that unrolls the current one using backtracking. Caution: exponential in branching.
    VG backtracking_unroll(uint32_t max_length, uint32_t max_depth,
                           map<id_t, pair<id_t, bool> >& node_translation);
    /// Represent the whole graph up to max_length across an inversion on the forward strand.
    VG unfold(uint32_t max_length,
              map<id_t, pair<id_t, bool> >& node_translation);
    /// Assume two node translations, the over is based on the under; merge them.
    map<id_t, pair<id_t, bool> > overlay_node_translations(const map<id_t, pair<id_t, bool> >& over,
                                                           const map<id_t, pair<id_t, bool> >& under);
    /// Use our topological sort to quickly break cycles in the graph, return the edges which are removed.
    /// Very non-optimal, but fast.
    vector<Edge> break_cycles(void);
    /// Remove pieces of the graph which are not part of any path.
    void remove_non_path(void);
    /// Convert edges that are both from_start and to_end to "regular" ones from end to start.
    void flip_doubly_reversed_edges(void);

    /// Build a graph from a GFA stream.
    void from_gfa(istream& in, bool showp = false);
    /// Build a graph from a Turtle stream.
    void from_turtle(string filename, string baseuri, bool showp = false);

    /// Destructor.
    ~VG(void);

    /// Copy constructor.
    VG(const VG& other) {
        init();
        if (this != &other) {
            // cleanup
            clear_indexes();
            // assign
            graph = other.graph;
            paths = other.paths;
            // re-index
            rebuild_indexes();
        }
    }

    /// Move constructor.
    VG(VG&& other) noexcept {
        init();
        graph = other.graph;
        paths = other.paths;
        other.graph.Clear();
        rebuild_indexes();
        // should copy over indexes
    }

    /// Copy assignment operator.
    VG& operator=(const VG& other) {
        VG tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    /// Move assignment operator.
    VG& operator=(VG&& other) noexcept {
        std::swap(graph, other.graph);
        rebuild_indexes();
        return *this;
    }

    // TODO: document all these

    void build_indexes(void);
    void build_node_indexes(void);
    void build_edge_indexes(void);
    void index_paths(void);
    void clear_node_indexes(void);
    void clear_node_indexes_no_resize(void);
    void clear_edge_indexes(void);
    void clear_edge_indexes_no_resize(void);
    void clear_indexes(void);
    void clear_indexes_no_resize(void);
    void resize_indexes(void);
    void rebuild_indexes(void);
    void rebuild_edge_indexes(void);

    /// Literally merge protobufs.
    void merge(Graph& g);
    /// Literally merge protobufs.
    void merge(VG& g);

    /// Clear the paths object (which indexes the graph.paths) and the graph paths themselves.
    void clear_paths(void);
    /// Synchronize in-memory indexes and protobuf graph.
    void sync_paths(void);

    /// Merge protobufs after removing overlaps.
    /// Good when there aren't many overlaps.
    void merge_union(VG& g);
    /// Helper to merge_union.
    void remove_duplicated_in(VG& g);
    /// Remove duplicated nodes and edges.
    void remove_duplicates(void);

    /// Limit the local complexity of the graph, connecting pruned components to a head and tail node
    /// depending on the direction which we come into the node when the edge_max is passed.
    void prune_complex_paths(int length, int edge_max, Node* head_node, Node* tail_node);
    void prune_short_subgraphs(size_t min_size);

    /// Write to a stream in chunked graphs.
    void serialize_to_ostream(ostream& out, id_t chunk_size = 1000);
    void serialize_to_file(const string& file_name, id_t chunk_size = 1000);

    // can we handle this with merge?
    //void concatenate(VG& g);

    /// Get the maximum node ID in the graph.
    id_t max_node_id(void);
    /// Get the minimum node ID in the graph.
    id_t min_node_id(void);
    /// Squish the node IDs down into as small a space as possible. Fixes up paths itself.
    void compact_ids(void);
    /// Add the given value to all node IDs. Preserves the paths.
    void increment_node_ids(id_t increment);
    /// Subtract the given value from all the node IDs. Must not create a node with 0 or negative IDs. Invalidates the paths.
    void decrement_node_ids(id_t decrement);
    /// Change the ID of the node with the first id to the second, new ID not
    /// used by any node. Invalidates any paths containing the node, since they
    /// are not updated.
    void swap_node_id(id_t node_id, id_t new_id);
    /// Change the ID of the given node to the second, new ID not used by any
    /// node. Invalidates the paths. Invalidates any paths containing the node,
    /// since they are not updated.
    void swap_node_id(Node* node, id_t new_id);

    /// Iteratively add when nodes and edges are novel. Good when there are very
    /// many overlaps. TODO: If you are using this with warn on duplicates on,
    /// and you know there shouldn't be any duplicates, maybe you should use
    /// merge instead.
    /// This version sorts paths on rank after adding in the path mappings from
    /// the other graph.
    void extend(VG& g, bool warn_on_duplicates = false);
    /// This version does not sort path mappings by rank. In order to preserve
    /// paths, call Paths::sort_by_mapping_rank() and
    /// Paths::rebuild_mapping_aux() after you are done adding in graphs to this
    /// graph.
    void extend(Graph& graph, bool warn_on_duplicates = false);
    // TODO: Do a member group for these overloads

    /// Add another graph into this graph, attaching tails to heads.
    /// Modify ids of the second graph to ensure we don't have conflicts.
    /// Then attach tails of this graph to the heads of the other, and extend(g).
    void append(VG& g);

    /// Add another graph into this graph.
    /// Don't append or join the nodes in the graphs;
    /// just ensure that ids are unique, then apply extend.
    void combine(VG& g);

    /// %Edit the graph to include the path.
    void include(const Path& path);

    /// %Edit the graph to include all the sequence and edges added by the given
    /// paths. Can handle paths that visit nodes in any orientation. Returns a
    /// vector of Translations, one per node existing after the edit, describing
    /// how each new or conserved node is embedded in the old graph. Note that
    /// this method sorts the graph and rebuilds the path index, so it should
    /// not be called in a loop.
    vector<Translation> edit(const vector<Path>& paths);
    
    /// %Edit the graph to include all the sequences and edges added by the
    /// given path. Returns a vector of Translations, one per original-node
    /// fragment. Completely novel nodes are not mentioned, and nodes with no
    /// Translations are assumed to be carried through unchanged. Invalidates
    /// the rank-based Paths index. Does not sort the graph. Suitable for
    /// calling in a loop.
    vector<Translation> edit_fast(const Path& path);

    /// Find all the points at which a Path enters or leaves nodes in the graph. Adds
    /// them to the given map by node ID of sets of bases in the node that will need
    /// to become the starts of new nodes.
    void find_breakpoints(const Path& path, map<id_t, set<pos_t>>& breakpoints);
    /// Take a map from node ID to a set of offsets at which new nodes should
    /// start (which may include 0 and 1-past-the-end, which should be ignored),
    /// break the specified nodes at those positions. Returns a map from old
    /// node start position to new node pointer in the graph. Note that the
    /// caller will have to crear and rebuild path rank data.
    ///
    /// Returns a map from old node start position to new node. This map
    /// contains some entries pointing to null, for positions past the ends of
    /// original nodes. It also maps from positions on either strand of the old
    /// node to the same new node pointer; the new node's forward strand is
    /// always the same as the old node's forward strand.
    map<pos_t, Node*> ensure_breakpoints(const map<id_t, set<pos_t>>& breakpoints);

    /// Flips the breakpoints onto the forward strand.
    map<id_t, set<pos_t>> forwardize_breakpoints(const map<id_t, set<pos_t>>& breakpoints);

    /// Given a path on nodes that may or may not exist, and a map from start
    /// position in the old graph to a node in the current graph, add all the
    /// new sequence and edges required by the path. The given path must not
    /// contain adjacent perfect match edits in the same mapping, or any
    /// deletions on the start or end of mappings (the removal of which can be
    /// accomplished with the Path::simplify() function).
    ///
    /// Outputs (and caches for subsequent calls) novel nodes in added_seqs, and
    /// Paths describing where novel nodes translate back to in the original
    /// graph in added_nodes. Also needs a map of the original sizes of nodes
    /// deleted from the original graph, for reverse complementing.
    void add_nodes_and_edges(const Path& path,
                             const map<pos_t, Node*>& node_translation,
                             map<pair<pos_t, string>, Node*>& added_seqs,
                             map<Node*, Path>& added_nodes,
                             const map<id_t, size_t>& orig_node_sizes);

    /// Produce a graph Translation object from information about the editing process.
    vector<Translation> make_translation(const map<pos_t, Node*>& node_translation,
                                         const map<Node*, Path>& added_nodes,
                                         const map<id_t, size_t>& orig_node_sizes);

    /// Add in the given node, by value.
    void add_node(const Node& node);
    /// Add in the given nodes, by value.
    void add_nodes(const vector<Node>& nodes);
    /// Add in the given edge, by value.
    void add_edge(const Edge& edge);
    /// Add in the given edges, by value.
    void add_edges(const vector<Edge>& edges);
    /// Add in the given edges, by value.
    void add_edges(const vector<Edge*>& edges);
    /// Add in the given nodes, by value.
    void add_nodes(const set<Node*>& nodes);
    /// Add in the given edges, by value.
    void add_edges(const set<Edge*>& edges);

    /// Count the number of nodes in the graph.
    id_t node_count(void);
    /// Count the number of edges in the graph.
    id_t edge_count(void);
    /// Get the total sequence length of nodes in the graph.
    /// TODO: redundant with length().
    id_t total_length_of_nodes(void);
    /// Get the rank of the node in the protobuf array that backs the graph.
    int node_rank(Node* node);
    /// Get the rank of the node in the protobuf array that backs the graph.
    int node_rank(id_t id);
    /// Get the number of edges attached to the start of a node.
    int start_degree(Node* node);
    /// Get the number of edges attached to the end of a node.
    int end_degree(Node* node);
    /// Get the number of edges attached to the left side of a NodeTraversal.
    int left_degree(NodeTraversal node);
    /// Get the number of edges attached to the right side of a NodeTraversal.
    int right_degree(NodeTraversal node);
    /// Get the edges of the specified node, and add them to the given vector.
    /// Guaranteed to add each edge only once per call.
    void edges_of_node(Node* node, vector<Edge*>& edges);
    /// Get the edges of the specified node.
    vector<Edge*> edges_of(Node* node);
    /// Get the edges from the specified node.
    vector<Edge*> edges_from(Node* node);
    /// Get the edges to the specified node
    vector<Edge*> edges_to(Node* node);
    /// Get the edges of the specified set of nodes, and add them to the given set of edge pointers.
    void edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges);

    /// Get the sides on the other side of edges to this side of the node.
    set<NodeSide> sides_to(NodeSide side);
    /// Get the sides on the other side of edges from this side of the node.
    set<NodeSide> sides_from(NodeSide side);
    /// Get the sides from both sides of the node.
    // TODO: what does this even mean?
    set<NodeSide> sides_from(id_t id);
    /// Get the sides to both sides of the node.
    set<NodeSide> sides_to(id_t id);
    /// Union of sides_to and sides_from.
    set<NodeSide> sides_of(NodeSide side);
    /// Get all sides connecting to this node.
    set<pair<NodeSide, bool>> sides_context(id_t node_id);
    /// Use sides_from an sides_to to determine if both nodes have the same context.
    bool same_context(id_t id1, id_t id2);
    /// Determine if the node is a prev ancestor of this one.
    bool is_ancestor_prev(id_t node_id, id_t candidate_id);
    /// Determine if the node is a prev ancestor of this one by trying to find it in a given number of steps.
    bool is_ancestor_prev(id_t node_id, id_t candidate_id, set<id_t>& seen, size_t steps = 64);
    /// Determine if the node is a next ancestor of this one.
    bool is_ancestor_next(id_t node_id, id_t candidate_id);
    /// Determine if the node is a next ancestor of this one by trying to find it in a given number of steps.
    bool is_ancestor_next(id_t node_id, id_t candidate_id, set<id_t>& seen, size_t steps = 64);
    /// Try to find a common ancestor by walking back up to steps from the first node.
    id_t common_ancestor_prev(id_t id1, id_t id2, size_t steps = 64);
    /// Try to find a common ancestor by walking forward up to steps from the first node
    id_t common_ancestor_next(id_t id1, id_t id2, size_t steps = 64);
    /// To-siblings are nodes which also have edges to them from the same nodes as this one.
    set<NodeTraversal> siblings_to(const NodeTraversal& traversal);
    /// From-siblings are nodes which also have edges to them from the same nodes as this one.
    set<NodeTraversal> siblings_from(const NodeTraversal& traversal);
    /// Full to-siblings are nodes traversals which share exactly the same upstream `NodeSide`s.
    set<NodeTraversal> full_siblings_to(const NodeTraversal& trav);
    /// Full from-siblings are nodes traversals which share exactly the same downstream `NodeSide`s.
    set<NodeTraversal> full_siblings_from(const NodeTraversal& trav);
    /// Get general siblings of a node.
    set<Node*> siblings_of(Node* node);
    /// Remove easily-resolvable redundancy in the graph.
    void simplify_siblings(void);
    /// Remove easily-resolvable redundancy in the graph for all provided to-sibling sets.
    void simplify_to_siblings(const set<set<NodeTraversal>>& to_sibs);
    /// Remove easily-resolvable redundancy in the graph for all provided from-sibling sets.
    void simplify_from_siblings(const set<set<NodeTraversal>>& from_sibs);
    /// Remove intransitive sibling sets, such as where (A, B, C) = S1 but C ∊ S2.
    set<set<NodeTraversal>> transitive_sibling_sets(const set<set<NodeTraversal>>& sibs);
    /// Remove sibling sets which don't have identical orientation.
    set<set<NodeTraversal>> identically_oriented_sibling_sets(const set<set<NodeTraversal>>& sibs);
    /// Determine if pos1 occurs directly before pos2.
    bool adjacent(const Position& pos1, const Position& pos2);

    /// Create a node. Use the VG class to generate ids.
    Node* create_node(const string& seq, id_t id = 0);
    /// Find a particular node.
    Node* get_node(id_t id);
    /// Get the subgraph of a node and all the edges it is responsible for
    /// (where it has the minimal ID) and add it into the given VG.
    void nonoverlapping_node_context_without_paths(Node* node, VG& g);
    /// Expand the context of what's already in the given graph by the given
    /// distance, either in nodes or in bases. Pulls material from this graph.
    void expand_context(VG& g, size_t distance, bool add_paths = true, bool use_steps = true);
    /// Expand the context of the given graph by the given number of steps. 
    void expand_context_by_steps(VG& g, size_t steps, bool add_paths = true);
    /// Expand the context of the given graph by the given number of bases. If
    /// reflect is true, bounce off the ends of nodes to get siblings of nodes
    /// you came from. Can take a set of NodeSides not to look out from, that
    /// act as barriers to context expansion. These barriers will have no edges
    /// attached to them in the final graph.
    void expand_context_by_length(VG& g, size_t length, bool add_paths = true,
        bool reflect = false, const set<NodeSide>& barriers = set<NodeSide>());
    /// Destroy the node at the given pointer. This pointer must point to a Node owned by the graph.
    void destroy_node(Node* node);
    /// Destroy the node with the given ID.
    void destroy_node(id_t id);
    /// Determine if the graph has a node with the given ID.
    bool has_node(id_t id);
    /// Determine if the graph contains the given node.
    bool has_node(Node* node);
    /// Determine if the graph contains the given node.
    bool has_node(const Node& node);
    /// Find a node with the given name, or create a new one if none is found.
    Node* find_node_by_name_or_add_new(string name);
    /// Run the given function on every node.
    void for_each_node(function<void(Node*)> lambda);
    /// Run the given function on every node in parallel.
    void for_each_node_parallel(function<void(Node*)> lambda);
    /// Go through all the nodes in the same connected component as the given node. Ignores relative orientation.
    void for_each_connected_node(Node* node, function<void(Node*)> lambda);
    
    /// Do a DFS search of the bidirected graph. A bidirected DFS starts at some
    /// root node, and traverses first all the nodes found reading out the right
    /// of that node in their appropriate relative orientations (including the
    /// root), and then all the nodes found reading left out of that node in
    /// their appropriate orientations (including the root). If any unvisited
    /// nodes are left in other connected components, the process will repeat
    /// from one such node, until all nodes have been visited in each
    /// orientation.
    void dfs(
        /// Called when node orientattion is first encountered.
        const function<void(NodeTraversal)>& node_begin_fn,
        /// Called when node orientation goes out of scope.
        const function<void(NodeTraversal)>& node_end_fn,
        /// Called to check if we should stop the DFS.
        const function<bool(void)>& break_fn,
        /// Called when an edge is encountered.
        const function<void(Edge*)>& edge_fn,
        /// Called when an edge forms part of the DFS spanning tree.
        const function<void(Edge*)>& tree_fn,
        /// Called when we meet an edge in the current tree component.
        const function<void(Edge*)>& edge_curr_fn,
        /// Called when we meet an edge in an already-traversed tree component.
        const function<void(Edge*)>& edge_cross_fn,
        /// Start only at these node traversals.
        const vector<NodeTraversal>* sources,
        /// When hitting a sink, don't keep walking.
        const set<NodeTraversal>* sinks);         

    /// Specialization of dfs for only handling nodes.
    void dfs(const function<void(NodeTraversal)>& node_begin_fn,
             const function<void(NodeTraversal)>& node_end_fn,
             const vector<NodeTraversal>* sources = NULL,
             const set<NodeTraversal>* sinks = NULL);         

    /// Specialization of dfs for only handling nodes + break function.
    void dfs(const function<void(NodeTraversal)>& node_begin_fn,
             const function<void(NodeTraversal)>& node_end_fn,
             const function<bool(void)>& break_fn);

    /// Is the graph empty?
    bool empty(void);

    /// Generate a digest of the serialized graph.
    const string hash(void);

    /// Remove nodes with no sequence.
    /// These are created in some cases during the process of graph construction.
    void remove_null_nodes(void);
    /// Remove a node but connect all of its predecessor and successor nodes with new edges.
    void remove_node_forwarding_edges(Node* node);
    /// Remove null nodes but connect predecessors and successors, preserving structure.
    void remove_null_nodes_forwarding_edges(void);
    /// Remove edges for which one of the nodes is not present.
    void remove_orphan_edges(void);
    /// Remove edges representing an inversion and edges on the reverse complement.
    void remove_inverting_edges(void);
    /// Determine if the graph has inversions.
    bool has_inverting_edges(void);

    /// Keep paths in the given set of path names. Populates kept_names with the names of the paths it actually found to keep.
    /// The paths specified may not overlap. Removes all nodes and edges not used by one of the specified paths.
    void keep_paths(set<string>& path_names, set<string>& kept_names);
    void keep_path(string& path_name);

    /// Path stats.
    /// Starting from offset in the first node, how many edges do we cross?
    /// path must be nonempty and longer than the given length. offset is
    /// interpreted as relative to the first node in its on-path
    /// orientation, and is inclusive.
    int path_edge_count(list<NodeTraversal>& path, int32_t offset, int path_length);
    /// Determine the offset in its last node at which the path starting at this offset in its first node ends.
    /// path must be nonempty and longer than the given length. offset is
    /// interpreted as relative to the first node in its on-path
    /// orientation, and is inclusive. Returned offset is remaining unused length
    /// in the last node touched.
    int path_end_node_offset(list<NodeTraversal>& path, int32_t offset, int path_length);
    /// Convert the stored paths in this graph to alignments.
    const vector<Alignment> paths_as_alignments(void);
    /// Return sequence string of path.
    const string path_sequence(const Path& path);
    /// Return percent identity between two paths (# matches / (#matches + #mismatches)).
    /// Note: uses ssw aligner, so will only work on small paths.
    double path_identity(const Path& path1, const Path& path2);
    /// Get the sequence of a NodeTraversal.
    string trav_sequence(const NodeTraversal& trav);

    /// Convert a VG graph to superbubble algorithm input format.
    SB_Input vg_to_sb_input();
    /// Find the superbubbles in the given input graph.
    vector<pair<id_t, id_t> > get_superbubbles(SB_Input sbi);
    /// Find the superbubbles in this graph.
    vector<pair<id_t, id_t> > get_superbubbles();

    //map<pair<id_t, id_t>, vector<id_t> > superbubbles(void);

    /// Takes in a pathname and the nucleotide position (like from a vcf) and
    /// returns the node id which contains that position.
    id_t get_node_at_nucleotide(string pathname, int nuc);

    // edges
    /// Create an edge.
    /// If the given edge cannot be created, returns null.
    /// If the given edge already exists, returns the existing edge.
    Edge* create_edge(Node* from, Node* to, bool from_start = false, bool to_end = false);
    /// Create an edge.
    /// If the given edge cannot be created, returns null.
    /// If the given edge already exists, returns the existing edge.
    Edge* create_edge(id_t from, id_t to, bool from_start = false, bool to_end = false);
    /// Make a left-to-right edge from the left NodeTraversal to the right one, respecting orientations.
    /// If the given edge cannot be created, returns null.
    /// If the given edge already exists, returns the existing edge.
    Edge* create_edge(NodeTraversal left, NodeTraversal right);
    /// Make an edge connecting the given sides of nodes.
    /// If the given edge cannot be created, returns null.
    /// If the given edge already exists, returns the existing edge.
    Edge* create_edge(NodeSide side1, NodeSide side2);

    /// Get a pointer to the specified edge.
    /// This can take sides in any order.
    Edge* get_edge(const NodeSide& side1, const NodeSide& side2);
    /// Get a pointer to the specified edge.
    /// This can take sides in any order.
    Edge* get_edge(const pair<NodeSide, NodeSide>& sides);
    /// Get the edge connecting the given oriented nodes in the given order.
    Edge* get_edge(const NodeTraversal& left, const NodeTraversal& right);
    /// Festroy the edge at the given pointer. This pointer must point to an edge owned by the graph.
    void destroy_edge(Edge* edge);
    /// Destroy the edge between the given sides of nodes. These can be in either order.
    void destroy_edge(const NodeSide& side1, const NodeSide& side2);
    /// Destroy the edge between the given sides of nodes. This can take sides in any order
    void destroy_edge(const pair<NodeSide, NodeSide>& sides);
    /// Remove an edge from the node side indexes, so it doesn't show up when you
    /// ask for the edges connected to the side of a node. Makes the edge
    /// untraversable until the indexes are rebuilt.
    void unindex_edge_by_node_sides(const NodeSide& side1, const NodeSide& side2);
    /// Remove an edge from the node side indexes, so it doesn't show up when you
    /// ask for the edges connected to the side of a node. Makes the edge
    /// untraversable until the indexes are rebuilt.
    void unindex_edge_by_node_sides(Edge* edge);
    /// Add an edge to the node side indexes. Doesn't touch the index of edges by
    /// node pairs or the graph; those must be updated seperately.
    void index_edge_by_node_sides(Edge* edge);
    /// Get the edge between the given node sides, which can be in either order.
    bool has_edge(const NodeSide& side1, const NodeSide& side2);
    /// Determine if the graph has an edge. This can take sides in any order.
    bool has_edge(const pair<NodeSide, NodeSide>& sides);
    /// Determine if the graph has an edge. This can take sides in any order.
    bool has_edge(Edge* edge);
    /// Determine if the graph has an edge. This can take sides in any order.
    bool has_edge(const Edge& edge);
    /// Determine if the graph has an inverting edge on the given node.
    bool has_inverting_edge(Node* n);
    /// Determine if the graph has an inverting edge from the given node.
    bool has_inverting_edge_from(Node* n);
    /// Determine if the graph has an inverting edge to the given node.
    bool has_inverting_edge_to(Node* n);
    /// Run the given function for each edge.
    void for_each_edge(function<void(Edge*)> lambda);
    /// Run the given function for each edge, in parallel.
    void for_each_edge_parallel(function<void(Edge*)> lambda);


    /// Circularize a subgraph / path using the head / tail nodes.
    void circularize(id_t head, id_t tail);
    void circularize(vector<string> pathnames);
    /// Connect node -> nodes.
    /// Connects from the right side of the first to the left side of the second.
    void connect_node_to_nodes(NodeTraversal node, vector<NodeTraversal>& nodes);
    /// Connect node -> nodes.
    /// You can optionally use the start of the first node instead of the end.
    void connect_node_to_nodes(Node* node, vector<Node*>& nodes, bool from_start = false);
    /// connect nodes -> node.
    /// Connects from the right side of the first to the left side of the second.
    void connect_nodes_to_node(vector<NodeTraversal>& nodes, NodeTraversal node);
    /// connect nodes -> node.
    // You can optionally use the end of the second node instead of the start.
    void connect_nodes_to_node(vector<Node*>& nodes, Node* node, bool to_end = false);

    // utilities
    // These only work on forward nodes.

    /// Divide a node at a given internal position. Inserts the new nodes in the
    /// correct paths, but can't update the ranks, so they need to be cleared and
    /// re-calculated by the caller.
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    /// Divide a node at a given internal position. This version works on a collection of internal positions, in linear time.
    void divide_node(Node* node, vector<int> positions, vector<Node*>& parts);
    /// Divide a path at a position. Also invalidates stored rank information.
    void divide_path(map<long, id_t>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    /// Convert the graph to Dot format.
    void to_dot(ostream& out,
                vector<Alignment> alignments = {},
                vector<Locus> loci = {},
                bool show_paths = false,
                bool walk_paths = false,
                bool annotate_paths = false,
                bool show_mappings = false,
                bool simple_mode = false,
                bool invert_edge_ports = false,
                bool color_variants = false,
                bool superbubble_ranking = false,
                bool superbubble_labeling = false,
                bool ultrabubble_labeling = false,
                bool skip_missing_nodes = false,
                int random_seed = 0);

    /// Convert the graph to Dot format.
    void to_dot(ostream& out, vector<Alignment> alignments = {}, bool show_paths = false, bool walk_paths = false,
                bool annotate_paths = false, bool show_mappings = false, bool invert_edge_ports = false, int random_seed = 0,
                bool color_variants = false);
    /// Convert the graph to GFA format.
    void to_gfa(ostream& out);
    /// Convert the graph to Turtle format.
    void to_turtle(ostream& out, const string& rdf_base_uri, bool precompress);
    /// Determine if the graph is valid or not, according to the specified criteria.
    bool is_valid(bool check_nodes = true,
                  bool check_edges = true,
                  bool check_paths = true,
                  bool check_orphans = true);

    /// Topologically order nodes.
    /// Makes sure that Nodes appear in the Protobuf Graph object in their topological sort order.
    void sort(void);
    /// Topological sort helper function, not really meant for external use.
    void topological_sort(deque<NodeTraversal>& l);
    /// Swap the given nodes. TODO: what does that mean?
    void swap_nodes(Node* a, Node* b);

    /// Use a topological sort to order and orient the nodes, and then flip some
    /// nodes around so that they are oriented the way they are in the sort.
    /// Populates nodes_flipped with the ids of the nodes that have had their
    /// orientations changed. TODO: update the paths that touch nodes that
    /// flipped around
    void orient_nodes_forward(set<id_t>& nodes_flipped);

    /// For each path, assign edits that describe a total match of the mapping to the node.
    void force_path_match(void);
    /// For each path, if a mapping has no edits then make it a perfect match against a node.
    /// This is the same as force_path_match, but only for empty mappings.
    void fill_empty_path_mappings(void);

    
    /// Align without base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const string& sequence,
                    Aligner* aligner,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    int8_t full_length_bonus = 0,
                    bool banded_global = false,
                    size_t max_span = 0,
                    bool print_score_matrices = false);
    /// Align without base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const Alignment& alignment,
                    Aligner* aligner,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    int8_t full_length_bonus = 0,
                    bool banded_global = false,
                    size_t max_span = 0,
                    bool print_score_matrices = false);
    
    /// Align with default Aligner.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const Alignment& alignment,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    int8_t full_length_bonus = 0,
                    bool banded_global = false,
                    size_t max_span = 0,
                    bool print_score_matrices = false);
    /// Align with default Aligner.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const string& sequence,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    int8_t full_length_bonus = 0,
                    bool banded_global = false,
                    size_t max_span = 0,
                    bool print_score_matrices = false);
    
    /// Align with base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align_qual_adjusted(const Alignment& alignment,
                                  QualAdjAligner* qual_adj_aligner,
                                  size_t max_query_graph_ratio = 0,
                                  bool pinned_alignment = false,
                                  bool pin_left = false,
                                  int8_t full_length_bonus = 0,
                                  bool banded_global = false,
                                  size_t max_span = 0,
                                  bool print_score_matrices = false);
    /// Align with base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align_qual_adjusted(const string& sequence,
                                  QualAdjAligner* qual_adj_aligner,
                                  size_t max_query_graph_ratio = 0,
                                  bool pinned_alignment = false,
                                  bool pin_left = false,
                                  int8_t full_length_bonus = 0,
                                  bool banded_global = false,
                                  size_t max_span = 0,
                                  bool print_score_matrices = false);
    
    


    /// Calls a function on all node-crossing paths with up to length across node boundaries.
    /// Considers each node in forward orientation to produce the kpaths around it.
    void for_each_kpath(int k, bool path_only, int edge_max,
                        function<void(NodeTraversal)> handle_prev_maxed,
                        function<void(NodeTraversal)> handle_next_maxed,
                        function<void(list<NodeTraversal>::iterator, list<NodeTraversal>&)> lambda);
    /// Calls a function on all kpaths of the given node.
    void for_each_kpath_parallel(int k, bool path_only, int edge_max,
                                 function<void(NodeTraversal)> handle_prev_maxed,
                                 function<void(NodeTraversal)> handle_next_maxed,
                                 function<void(list<NodeTraversal>::iterator, list<NodeTraversal>&)> lambda);
    /// Calls a function on all node-crossing paths with up to length across node boundaries.
    /// Considers each node in forward orientation to produce the kpaths around it.
    void for_each_kpath(int k, bool path_only, int edge_max,
                        function<void(NodeTraversal)> handle_prev_maxed,
                        function<void(NodeTraversal)> handle_next_maxed,
                        function<void(size_t,Path&)> lambda);
    /// Calls a function on all kpaths of the given node.
    void for_each_kpath_parallel(int k, bool path_only, int edge_max,
                                 function<void(NodeTraversal)> handle_prev_maxed,
                                 function<void(NodeTraversal)> handle_next_maxed,
                                 function<void(size_t,Path&)> lambda);
    /// Calls a function on all kpaths of the given node.
    void for_each_kpath_of_node(Node* node, int k, bool path_only, int edge_max,
                                function<void(NodeTraversal)> handle_prev_maxed,
                                function<void(NodeTraversal)> handle_next_maxed,
                                function<void(list<NodeTraversal>::iterator, list<NodeTraversal>&)> lambda);
    /// Calls a function on all kpaths of the given node.
    void for_each_kpath_of_node(Node* n, int k, bool path_only, int edge_max,
                                function<void(NodeTraversal)> handle_prev_maxed,
                                function<void(NodeTraversal)> handle_next_maxed,
                                function<void(size_t,Path&)> lambda);

    /// Get kpaths. TODO: what is this for?
    void kpaths(set<list<NodeTraversal> >& paths, int length, bool path_only, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    /// Get kpaths. TODO: what is this for?
    void kpaths(vector<Path>& paths, int length, bool path_only, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);

    /// Get kpaths on a particular node. TODO: what is this for?
    void kpaths_of_node(Node* node, set<list<NodeTraversal> >& paths,
                        int length, bool path_only, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    /// Get kpaths on a particular node. TODO: what is this for?
    void kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, bool path_only, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    /// Get kpaths on a particular node. TODO: what is this for?
    void kpaths_of_node(id_t node_id, vector<Path>& paths, int length, bool path_only, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    /// Given an oriented start node, a length in bp, a maximum number of edges
    /// to cross, and a stack of nodes visited so far, fill in the set of paths
    /// with all the paths starting at the oriented start node and going left off
    /// its end no longer than the specified length, calling maxed_nodes on nodes
    /// which can't be visited due to the edge-crossing limit. Produces paths
    /// ending with the specified node. TODO: postfix should not be (potentially)
    /// copied on every call.
    void prev_kpaths_from_node(NodeTraversal node, int length, bool path_only, int edge_max, bool edge_bounding,
                               list<NodeTraversal> postfix, set<list<NodeTraversal> >& walked_paths,
                               const vector<string>& followed_paths,
                               function<void(NodeTraversal)>& maxed_nodes);
    /// Do the same as prec_kpaths_from_node, except going right, producing a path starting with the specified node.
    void next_kpaths_from_node(NodeTraversal node, int length, bool path_only, int edge_max, bool edge_bounding,
                               list<NodeTraversal> prefix, set<list<NodeTraversal> >& walked_paths,
                               const vector<string>& followed_paths,
                               function<void(NodeTraversal)>& maxed_nodes);

    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(id_t from, id_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    // traversal
    /// Get the nodes attached to the left side of the given NodeTraversal, in their proper orientations.
    void nodes_prev(NodeTraversal n, vector<NodeTraversal>& nodes);
    /// Get the nodes attached to the left side of the given NodeTraversal, in their proper orientations.
    vector<NodeTraversal> nodes_prev(NodeTraversal n);
    /// Get traversals before this node on the same strand. Same as nodes_prev but using set.
    set<NodeTraversal> travs_to(NodeTraversal node);

    /// Get the nodes attached to the right side of the given NodeTraversal, in their proper orientations.
    void nodes_next(NodeTraversal n, vector<NodeTraversal>& nodes);
    /// Get the nodes attached to the right side of the given NodeTraversal, in their proper orientations.
    vector<NodeTraversal> nodes_next(NodeTraversal n);
    /// Get traversals after this node on the same strand. Same as nodes_next but using set.
    set<NodeTraversal> travs_from(NodeTraversal node);

    /// Get traversals either before or after this node on the same strand.
    set<NodeTraversal> travs_of(NodeTraversal node);

    /// Count the nodes attached to the left side of the given NodeTraversal.
    int node_count_prev(NodeTraversal n);
    /// Count the nodes attached to the right side of the given NodeTraversal.
    int node_count_next(NodeTraversal n);

    // paths
    /// Create a path.
    Path create_path(const list<NodeTraversal>& nodes);
    /// Create a path.
    Path create_path(const vector<NodeTraversal>& nodes);
    /// Get the string sequence for all the NodeTraversals on the given path.
    string path_string(const list<NodeTraversal>& nodes);
    /// Get the string sequence for traversing the given path.
    /// Assumes the path covers the entirety of any nodes visited. Handles backward nodes.
    string path_string(const Path& path);
    /// Expand a path. TODO: what does that mean?
    void expand_path(const list<NodeTraversal>& path, vector<NodeTraversal>& expanded);
    /// Fill in the node_start map with the first index along the path at which each node appears.
    /// Caller is responsible for dealing with orientations.
    void node_starts_in_path(const list<NodeTraversal>& path,
                             map<Node*, int>& node_start);
    /// Return true if nodes share all paths and the mappings they share in these paths
    /// are adjacent, in the specified relative order and orientation.
    bool nodes_are_perfect_path_neighbors(NodeTraversal left, NodeTraversal right);
    /// Return true if the mapping completely covers the node it maps to and is a perfect match.
    bool mapping_is_total_match(const Mapping& m);
    /// Concatenate the mappings for a pair of nodes; handles multiple mappings per path.
    map<string, vector<Mapping>> concat_mappings_for_node_pair(id_t id1, id_t id2);
    /// Concatenate mappings for a list of nodes that we want to concatenate.
    /// Returns, for each path name, a vector of merged mappings, once per path
    /// traversal of the run of nodes. Those merged mappings are in the
    /// orientation of the merged node (so mappings to nodes that are traversed
    /// in reverse will have their flags toggled). We assume that all mappings on
    /// the given nodes are full-length perfect matches, and that all the nodes
    /// are perfect path neighbors.
    map<string, vector<Mapping>> concat_mappings_for_nodes(const list<NodeTraversal>& nodes);

    /// Expand a path. TODO: what does that mean?
    /// These versions handle paths in which nodes can be traversed multiple
    /// times. Unfortunately since we're throwing non-const iterators around, we
    /// can't take the input path as const.
    void expand_path(list<NodeTraversal>& path, vector<list<NodeTraversal>::iterator>& expanded);
    /// Find node starts in a path. TODO: what does that mean?
    /// To get the starts out of the map this produces, you need to dereference
    /// the iterator and then get the address of the NodeTraversal (stored in the
    /// list) that you are talking about.
    void node_starts_in_path(list<NodeTraversal>& path,
                             map<NodeTraversal*, int>& node_start);

    // kmers
    /// Call a function for each kmer in the graph, in parallel.
    void for_each_kmer_parallel(int kmer_size,
                                bool path_only,
                                int edge_max,
                                function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                                int stride = 1,
                                bool allow_dups = false,
                                bool allow_negatives = false);
    /// Call a function for each kmer in the graph.
    void for_each_kmer(int kmer_size,
                       bool path_only,
                       int edge_max,
                       function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                       int stride = 1,
                       bool allow_dups = false,
                       bool allow_negatives = false);
    /// Call a function for each kmer on a node.
    void for_each_kmer_of_node(Node* node,
                               int kmer_size,
                               bool path_only,
                               int edge_max,
                               function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                               int stride = 1,
                               bool allow_dups = false,
                               bool allow_negatives = false);

    /// For the given kmer of the given length starting at the given
    /// offset into the given Node along the given path, fill in end_node and
    /// end_offset with where the end of the kmer falls (counting from the right
    /// side of the NodeTraversal), prev_chars with the characters that preceed
    /// it, next_chars with the characters that follow it, prev_ and
    /// next_positions with the ((node ID, orientation), offset) pairs of the
    /// places you can come from/go next (from the right end of the kmer).
    /// Refuses to follow more than edge_max edges. Offsets are in the path
    /// orientation. Meant for gcsa2.
    void kmer_context(string& kmer,
                      int kmer_size,
                      bool path_only,
                      int edge_max,
                      bool forward_only,
                      list<NodeTraversal>& path,
                      list<NodeTraversal>::iterator start_node,
                      int32_t start_offset,
                      list<NodeTraversal>::iterator& end_node,
                      int32_t& end_offset,
                      set<tuple<char, id_t, bool, int32_t>>& prev_positions,
                      set<tuple<char, id_t, bool, int32_t>>& next_positions);

    /// Do the GCSA2 kmers for a node. head_node and tail_node must both be non-
    /// null, but only one of those nodes actually needs to be in the graph. They
    /// will be examined directly to get their representative characters. They
    /// also don't need to be actually owned by the graph; they can be copies.
    void gcsa_handle_node_in_graph(Node* node, int kmer_size, bool path_only,
                                   int edge_max, int stride,
                                   bool forward_only,
                                   Node* head_node, Node* tail_node,
                                   function<void(KmerPosition&)> lambda);

    /// Call a function for each GCSA2 kemr position in parallel.
    /// GCSA kmers are the kmers in the graph with each node
    /// existing in both its forward and reverse-complement orientation. Node IDs
    /// in the GCSA graph are 2 * original node ID, +1 if the GCSA node
    /// represents the reverse complement, and +0 if it does not. Non-reversing
    /// edges link the forward copy of the from node to the forward copy of the
    /// to node, and similarly for the reverse complement copies, while reversing
    /// edges link the forward copy of the from node to the *reverse complement*
    /// copy of the to node, and visa versa. This allows us to index both the
    /// forward and reverse strands of every node, and to deal with GCSA's lack
    /// of support for reversing edges, with the same trick. Note that
    /// start_tail_id, if zero, will be replaced with the ID actually used for the
    /// start/end node before lambda is ever called.
    void for_each_gcsa_kmer_position_parallel(int kmer_size, bool path_only,
                                              int edge_max, int stride,
                                              bool forward_only,
                                              id_t& head_id, id_t& tail_id,
                                              function<void(KmerPosition&)> lambda);
    
    /// Get the GCSA2 kmers in the graph.
    void get_gcsa_kmers(int kmer_size, bool path_only,
                        int edge_max, int stride,
                        bool forward_only,
                        const function<void(vector<gcsa::KMer>&, bool)>& handle_kmers,
                        id_t& head_id, id_t& tail_id);

    /// Writhe the GCSA2 kmer file for the graph to the goven stream.
    void write_gcsa_kmers(int kmer_size, bool path_only,
                          int edge_max, int stride,
                          bool forward_only,
                          ostream& out,
                          id_t& head_id, id_t& tail_id);

    /// Write the GCSA2 kmers to a temp file with the given base. Return the name of the file.
    string write_gcsa_kmers_to_tmpfile(int kmer_size,
                                       bool paths_only,
                                       bool forward_only,
                                       id_t& head_id, id_t& tail_id,
                                       size_t doubling_steps = 2,
                                       size_t size_limit = 200,
                                       const string& base_file_name = ".vg-kmers-tmp-");

    /// Construct the GCSA2 index for this graph.
    void build_gcsa_lcp(gcsa::GCSA*& gcsa,
                        gcsa::LCPArray*& lcp,
                        int kmer_size,
                        bool paths_only,
                        bool forward_only,
                        size_t doubling_steps = 2,
                        size_t size_limit = 200,
                        const string& base_file_name = ".vg-kmers-tmp-");

    /// Take all nodes that would introduce paths of > edge_max edge crossings, remove them, and link their neighbors to
    /// head_node or tail_node depending on which direction the path extension was stopped.
    /// For pruning graph prior to indexing with gcsa2.
    void prune_complex(int path_length, int edge_max, Node* head_node, Node* tail_node);
    /// Wrap the graph with heads and tails before doing the prune.
    /// Utility function for preparing for indexing.
    void prune_complex_with_head_tail(int path_length, int edge_max);

private:
    /// Call the given function on each kmer. If parallel is specified, goes
    /// through nodes one per thread. If node is not null, looks only at kmers of
    /// that specific node.
    void _for_each_kmer(int kmer_size,
                        bool path_only,
                        int edge_max,
                        function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                        bool parallel,
                        int stride,
                        bool allow_dups,
                        bool allow_negatives,
                        Node* node = nullptr);
    
    /// Private method to funnel other align functions into. max_span specifies
    /// the min distance to unfold the graph to, and is meant to be the longest
    /// path that the specified sequence could cover, accounting for deletions.
    /// If it's less than the sequence's length, the sequence's length is used.
    Alignment align(const Alignment& alignment,
                    Aligner* aligner,
                    QualAdjAligner* qual_adj_aligner,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    int8_t full_length_bonus = 0,
                    bool banded_global = false,
                    size_t max_span = 0,
                    bool print_score_matrices = false);


public:

    /// Generate random reads.
    /// Note that even if either_strand is false, having backward nodes in the
    /// graph will result in some reads from the global reverse strand.
    Alignment random_read(size_t read_len, mt19937& rng, id_t min_id, id_t max_id, bool either_strand);

    /// Find subgraphs.
    void disjoint_subgraphs(list<VG>& subgraphs);
    /// Get the head nodes (nodes with edges only to their right sides). These are required to be oriented forward.
    void head_nodes(vector<Node*>& nodes);
    /// Get the head nodes (nodes with edges only to their right sides). These are required to be oriented forward.
    vector<Node*> head_nodes(void);
    /// Determine if a node is a head node.
    bool is_head_node(id_t id);
    /// Determine if a node is a head node.
    bool is_head_node(Node* node);
    /// Get the distance from head of node to beginning of graph, or -1 if limit exceeded.
    int32_t distance_to_head(NodeTraversal node, int32_t limit = 1000);
    /// Get the distance from head of node to beginning of graph, or -1 if limit exceeded.
    int32_t distance_to_head(NodeTraversal node, int32_t limit,
                             int32_t dist, set<NodeTraversal>& seen);
    /// Get the tail nodes (nodes with edges only to their left sides). These are required to be oriented forward.
    vector<Node*> tail_nodes(void);
    /// Get the tail nodes (nodes with edges only to their left sides). These are required to be oriented forward.
    void tail_nodes(vector<Node*>& nodes);
    /// Determine if a node is a tail node.
    bool is_tail_node(id_t id);
    /// Determine if a node is a tail node.
    bool is_tail_node(Node* node);
    /// Get the distance from tail of node to end of graph, or -1 if limit exceeded.
    int32_t distance_to_tail(NodeTraversal node, int32_t limit = 1000);
    /// Get the distance from tail of node to end of graph, or -1 if limit exceeded.
    int32_t distance_to_tail(NodeTraversal node, int32_t limit,
                             int32_t dist, set<NodeTraversal>& seen);
    /// Get the distance from tail of node to end of graph, or -1 if limit exceeded.
    int32_t distance_to_tail(id_t id, int32_t limit = 1000);
    /// Collect the subgraph of a Node. TODO: what does that mean?
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    /// Join head nodes of graph to common null node, creating a new single head.
    Node* join_heads(void);
    /// Join head nodes of graph to specified node. Optionally from the start/to the end of the new node.
    void join_heads(Node* node, bool from_start = false);
    /// Join tail nodes of graph to specified node. Optionally from the start/to the end of the new node.
    void join_tails(Node* node, bool to_end = false);

    /// Add singular head and tail null nodes to graph.
    void wrap_with_null_nodes(void);

    /// Add a start node and an end node, where all existing heads in the graph
    /// are connected to the start node, and all existing tails in the graph are
    /// connected to the end node. Any connected components in the graph which do
    /// not have either are connected to the start at an arbitrary point, and the
    /// end node from nodes going to that arbitrary point. If start_node or
    /// end_node is null, a new node will be created. Otherwise, the passed node
    /// will be used. Note that this visits every node, to make sure it is
    /// attached to all connected components. Note that if a graph has, say,
    /// heads but no tails, the start node will be attached buut the end node
    /// will be free-floating.
    void add_start_end_markers(int length,
                               char start_char, char end_char,
                               Node*& start_node, Node*& end_node,
                               id_t start_id = 0, id_t end_id = 0);

    /// Structure for managing parallel construction of a graph.
    // TODO: delete this since we don't use it anymore.
    struct Plan {
        VG* graph;
        map<long, vector<vcflib::VariantAllele> > alleles;
        // What alleles are visited by phasing paths? For each position and
        // allele index, stores a vector of flags, one per phase path.
        map<pair<long, int>, vector<bool>> phase_visits;
        // What alleles are visited by paths defining the alts of a variant? For
        // each position and allele index, stores a vector of variant ID, alt
        // number pairs.
        map<pair<long, int>, vector<pair<string, int>>> variant_alts;
        string seq;
        string name;
        // Make a new plan, moving the alleles map and phase visit vector map
        // into the plan.
        Plan(VG* graph,
             map<long, vector<vcflib::VariantAllele> >&& alleles,
             map<pair<long, int>, vector<bool>>&& phase_visits,
             map<pair<long, int>, vector<pair<string, int>>>&& variant_alts,
             string seq,
             string name)
            : graph(graph)
            , alleles(std::move(alleles))
            , phase_visits(std::move(phase_visits))
            , variant_alts(std::move(variant_alts))
            , seq(seq)
            , name(name) { };
    };


private:

    void init(void); ///< setup, ensures that gssw == NULL on startup
    /// Placeholder for functions that sometimes need to be passed an empty vector
    vector<id_t> empty_ids;
    /// Placeholder for functions that sometimes need to be passed an empty vector
    vector<pair<id_t, bool>> empty_edge_ends;

};

} // end namespace vg

#endif

#ifndef VG_H

#define VG_H

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
#include "gcsa.h"
#include "lcp.h"
#include "gssw_aligner.hpp"
#include "region.hpp"
#include "path.hpp"
#include "utility.hpp"
#include "alignment.hpp"

#include "vg.pb.h"
#include "hash_map.hpp"

#include "progress_bar.hpp"
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
#include "DetectSuperBubble.hpp"
#include "helperDefs.hpp"


// uncomment to enable verbose debugging to stderr
//#define debug

namespace vg {


// Represents a node traversed in a certain orientation. The default orientation
// is start to end, but if `backward` is set, represents the node being
// traversed end to start. A list of these can serve as an edit-free version of
// a path, especially if supplemented with a length and an initial node offset.
// A path node has a left and a right side, which are the start and end of the
// node if it is forward, or the end and start of the node if it is backward.
class NodeTraversal {
public:
    Node* node;
    bool backward;

    inline NodeTraversal(Node* node, bool backward = false): node(node), backward(backward) {
        // Nothing to do
    }

    inline NodeTraversal(): NodeTraversal(nullptr) {
        // Nothing to do
    }

    inline bool operator==(const NodeTraversal& other) const {
        return node == other.node && backward == other.backward;
    }

    inline bool operator!=(const NodeTraversal& other) const {
        return node != other.node || backward != other.backward;
    }

    inline bool operator<(const NodeTraversal& other) const {
        return node < other.node || (node == other.node && backward < other.backward);
    }

    // reverse complement the node traversal
    inline NodeTraversal reverse(void) const {
        return NodeTraversal(node, !backward);
    }

};

inline ostream& operator<<(ostream& out, const NodeTraversal& nodetraversal) {
    return out << nodetraversal.node->id() << " " << (nodetraversal.backward ? "rev" : "fwd");
}

// Represents one side of a Node, identified by ID, for the purposes of
// indexing edges.
class NodeSide {
public:
    id_t node;
    bool is_end;

    // We need this to be a converting constructor so we can represent the empty and deleted item keys in a pair_hash_map.
    inline NodeSide(id_t node, bool is_end = false): node(node), is_end(is_end) {
        // Nothing to do
    }

    inline NodeSide(): NodeSide(0, false) {
        // Nothing to do
    }

    inline bool operator==(const NodeSide& other) const {
        return node == other.node && is_end == other.is_end;
    }

    inline bool operator!=(const NodeSide& other) const {
        return node != other.node || is_end != other.is_end;
    }

    inline bool operator<(const NodeSide& other) const {
        return node < other.node || (node == other.node && is_end < other.is_end);
    }

    // Make an edge into a canonically ordered pair of NodeSides
    static inline pair<NodeSide, NodeSide> pair_from_edge(Edge* e) {
        return minmax(NodeSide(e->from(), !e->from_start()), NodeSide(e->to(), e->to_end()));
    }

    // Make an edge into a canonically ordered pair of NodeSides
    static inline pair<NodeSide, NodeSide> pair_from_edge(const Edge& e) {
        return minmax(NodeSide(e.from(), !e.from_start()), NodeSide(e.to(), e.to_end()));
    }

    // Make a canonically ordered pair of NodeSides from an edge off of the
    // start of a node, to another node in the given relative orientation.
    static inline pair<NodeSide, NodeSide> pair_from_start_edge(id_t start_id, const pair<id_t, bool>& oriented_other) {
        // If it's in the same relative orientation, we go to its end.
        return minmax(NodeSide(start_id, false), NodeSide(oriented_other.first, !oriented_other.second));
    }

    // Make a canonically ordered pair of NodeSides from an edge off of the
    // end of a node, to another node in the given relative orientation.
    static inline pair<NodeSide, NodeSide> pair_from_end_edge(id_t end_id, const pair<id_t, bool>& oriented_other) {
        // If it's in the same relative orientation, we go to its start.
        return minmax(NodeSide(end_id, true), NodeSide(oriented_other.first, oriented_other.second));
    }

    // reverse complement the node side
    inline NodeSide flip(void) const {
        return NodeSide(node, !is_end);
    }

};

// helpers to be more clear
NodeSide node_start(id_t id);
NodeSide node_end(id_t id);

// We create a struct that represents each kmer record we want to send to gcsa2
struct KmerPosition {
    string kmer;
    string pos;
    set<char> prev_chars;
    set<char> next_chars;
    set<string> next_positions;
};

struct SB_Input{
    int num_vertices;
    vector<pair<id_t, id_t> > edges;
};

inline ostream& operator<<(ostream& out, const NodeSide& nodeside) {
    return out << nodeside.node << " " << (nodeside.is_end ? "end" : "start");
}

}

// Go up to std namespace for a moment
namespace std {
// We need to implement a hash function for these if we want to be able to use them in keys.
template <> struct hash<vg::NodeSide>
{
    // Produce a hash of a NodeSide
    size_t operator()(const vg::NodeSide& item) const
    {
        // Hash it just as we would a pair.
        return hash<pair<id_t, bool>>()(make_pair(item.node, item.is_end));
    }
};
}

namespace vg {

// Represents a sequence graph. Graphs consist of nodes, connected by edges.
// Graphs are bidirected and may be cyclic. Nodes carry forward-oriented
// sequences. Edges are directed, with a "from" and to" node, and are generally
// used to connect the end of the "from" node to the start of the "to" node.
// However, edges can connect to either the start or end of either node, in
// general, as long as they do not allow the same node to be visited twice along
// a path. Graphs have "head" and "tail" nodes, which are overall at the
// left/right of the graph, with nothing before/after them. Because otherwise
// identifying these nodes (i.e. classifying a terminal node as a head or a
// tail) would require a topological sort, we require that all head and tail
// nodes be in the same relative orientation. Head nodes must have edges only to
// their right sides, and tail nodes must have edges only to their left sides.
// There must be no possible path in the graph containing two head nodes or two
// tail nodes.
class VG {

public:

    // protobuf-based representation
    // NB: we can't subclass this safely, so it's best as a member
    Graph graph;

    // manages paths of the graph
    // initialized by setting paths._paths = graph.paths
    Paths paths;

    // name
    string name;

    // current id
    id_t current_id;
    // todo
    //id_t min_id;
    //id_t max_id;

    // nodes by id
    hash_map<id_t, Node*> node_by_id;

    // edges by sides of nodes they connect
    // Since duplicate edges are not permitted, two edges cannot connect the same pair of node sides.
    // Each edge is indexed here with the smaller NodeSide first. The actual node order is recorded in the Edge object.
    pair_hash_map<pair<NodeSide, NodeSide>, Edge*> edge_by_sides;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    hash_map<Node*, int> node_index;

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    hash_map<Edge*, int> edge_index;

    // edges indexed by nodes they connect
    // Stores the destinations and backward flags for edges attached to the starts of nodes (whether that node is "from" or "to").
    hash_map<id_t, vector<pair<id_t, bool>>> edges_on_start;
    // Stores the destinations and backward flags for edges attached to the ends of nodes (whether that node is "from" or "to").
    hash_map<id_t, vector<pair<id_t, bool>>> edges_on_end;

    // Set the edge indexes through this function. Picks up the sides being
    // connected by the edge automatically, and silently drops the edge if they
    // are already connected.
    void set_edge(Edge*);
    void print_edges(void);

    // access the edge indexes through these functions
    // Get nodes and backward flags following edges that attach to this node's start
    vector<pair<id_t, bool>>& edges_start(Node* node);
    vector<pair<id_t, bool>>& edges_start(id_t id);
    // Get nodes and backward flags following edges that attach to this node's end
    vector<pair<id_t, bool>>& edges_end(Node* node);
    vector<pair<id_t, bool>>& edges_end(id_t id);
    // properties of the graph
    size_t size(void); // number of nodes
    size_t length(void);

    // clear everything
    //void clear(void);

    // constructors

    // default
    VG(void);

    // construct from protobufs
    VG(istream& in, bool showp = false);

    // construct from an arbitrary source of Graph protobuf messages (which
    // populates the given Graph and returns a flag for whether it's valid).
    VG(function<bool(Graph&)>& get_next_graph, bool showp = false);

    // construct from sets of nodes and edges (e.g. subgraph of another graph)
    VG(set<Node*>& nodes, set<Edge*>& edges);

    // construct from VCF
    VG(vcflib::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target,
       bool target_is_chrom,
       int vars_per_region,
       int max_node_size = 0,
       bool flat_input_vcf = false,
       bool load_phasing_paths = false,
       bool load_variant_alt_paths = false,
       bool showprog = false);
       
    // Build the graph from a bunch of alleles, organized by position.
    void from_alleles(const map<long, vector<vcflib::VariantAllele> >& altp,
                      const map<pair<long, int>, vector<bool>>& visits,
                      size_t num_phasings,
                      const map<pair<long, int>, vector<pair<string, int>>>& variant_alts,
                      string& seq,
                      string& chrom);
    void vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, vector<vcflib::VariantAllele> >& altp,
                                map<pair<long, int>, vector<bool>>* phase_visits,
                                map<pair<long, int>, vector<pair<string, int>>>* alt_allele_visits,
                                bool flat_input_vcf = false);
    void slice_alleles(map<long, vector<vcflib::VariantAllele> >& altp,
                       int start_pos,
                       int stop_pos,
                       int max_node_size);
    // chops up the nodes
    void dice_nodes(int max_node_size);
    // does the reverse --- combines nodes by removing edges where doing so has no effect on the graph labels
    void unchop(void);
    // the set of components that could be merged into single nodes without
    // changing the path space of the graph
    set<list<Node*>> simple_components(int min_size = 1);
    set<list<Node*>> simple_multinode_components(void);
    // strongly connected components of the graph
    set<set<id_t> > strongly_connected_components(void);
    // only multi-node components
    set<set<id_t> > multinode_strongly_connected_components(void);
    // true if the graph does not contain cycles
    bool is_acyclic(void);
    // removes all elements which are not in a strongly connected component
    void keep_multinode_strongly_connected_components(void);
    // does the specified node traversal self-loop?
    bool is_self_looping(NodeTraversal trav);
    // simple cycles following Johnson's elementary cycles algorithm
    set<list<NodeTraversal> > elementary_cycles(void);
    // concatenates the nodes into a new node with the same external linkage as the provided component
    Node* concat_nodes(const list<Node*>& nodes);
    // merge the nodes into a single node, preserving external linkages
    // use the sequence of the first node as the basis
    Node* merge_nodes(const list<Node*>& nodes);
    // uses unchop and sibling merging to simplify the graph into a normalized form
    void normalize(void);
    // turn the graph into a dag by copying strongly connected components expand_scc_steps times
    // and translating the edges in the component to flow through the copies in one direction
    VG dagify(uint32_t expand_scc_steps,
              map<id_t, pair<id_t, bool> >& node_translation,
              size_t target_min_walk_length = 0,
              size_t component_length_max = 0);
    // generate a new graph that unrolls the current one using backtracking (caution: exponential in branching)
    VG backtracking_unroll(uint32_t max_length, uint32_t max_depth,
                           map<id_t, pair<id_t, bool> >& node_translation);
    // represents the whole graph up to max_length across an inversion on the forward strand
    VG unfold(uint32_t max_length,
              map<id_t, pair<id_t, bool> >& node_translation);
    // assume two node translations, the over is based on the under; merge them
    map<id_t, pair<id_t, bool> > overlay_node_translations(const map<id_t, pair<id_t, bool> >& over,
                                                           const map<id_t, pair<id_t, bool> >& under);
    // use our topological sort to quickly break cycles in the graph, return the edges which are removed
    // very non-optimal, but fast
    vector<Edge> break_cycles(void);
    // removes pieces of the graph which are not part of any path
    void remove_non_path(void);
    // converts edges that are both from_start and to_end to "regular" ones from end to start
    void flip_doubly_reversed_edges(void);

    void from_gfa(istream& in, bool showp = false);
    void from_turtle(string filename, string baseuri, bool showp = false);

    // default constructor, destructor
    ~VG(void);

    // copy constructor
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

    // move constructor
    VG(VG&& other) noexcept {
        init();
        graph = other.graph;
        paths = other.paths;
        other.graph.Clear();
        rebuild_indexes();
        // should copy over indexes
    }

    // copy assignment operator
    VG& operator=(const VG& other) {
        VG tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    // move assignment operator
    VG& operator=(VG&& other) noexcept {
        std::swap(graph, other.graph);
        rebuild_indexes();
        return *this;
    }

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

    // literally merge protobufs
    void merge(Graph& g);
    void merge(VG& g);

    // clear the paths object (which indexes the graph.paths) and the graph paths themselves
    void clear_paths(void);
    // synchronizes in-memory indexes and protobuf graph
    void sync_paths(void);

    // merge protobufs after removing overlaps
    // good when there aren't many overlaps
    void merge_union(VG& g);
    // helper to merge_union
    void remove_duplicated_in(VG& g);

    // limit the local complexity of the graph, connecting pruned components to a head and tail node
    // depending on the direction which we come into the node when the edge_max is passed
    void prune_complex_paths(int length, int edge_max, Node* head_node, Node* tail_node);
    void prune_short_subgraphs(size_t min_size);

    // write to a stream in chunked graphs
    void serialize_to_ostream(ostream& out, id_t chunk_size = 1000);
    void serialize_to_file(const string& file_name, id_t chunk_size = 1000);

    // can we handle this with merge?
    //void concatenate(VG& g);

    id_t max_node_id(void);
    id_t min_node_id(void);
    // Squish the node IDs down into as small a space as possible. Fixes up paths itself.
    void compact_ids(void);
    // Add the given value to all node IDs. Preserves the paths.
    void increment_node_ids(id_t increment);
    // Subtract the given value from all the node IDs. Must not create a node with 0 or negative IDs. Invalidates the paths.
    void decrement_node_ids(id_t decrement);
    // Change the ID of the node with the first id to the second, new ID not
    // used by any node. Invalidates any paths containing the node, since they
    // are not updated.
    void swap_node_id(id_t node_id, id_t new_id);
    // Change the ID of the given node to the second, new ID not used by any
    // node. Invalidates the paths. Invalidates any paths containing the node,
    // since they are not updated.
    void swap_node_id(Node* node, id_t new_id);

    // Iteratively add when nodes and edges are novel. Good when there are very
    // many overlaps. TODO: If you are using this with warn on duplicates on,
    // and you know there shouldn't be any duplicates, maybe you should use
    // merge instead.
    
    // This version sorts paths on rank after adding in the path mappings from
    // the other graph.
    void extend(VG& g, bool warn_on_duplicates = false);
    // This version does not sort path mappings by rank. In order to preserve
    // paths, call paths.sort_by_mapping_rank() and paths.rebuild_mapping_aux()
    // after you are done adding in graphs to this graph.
    void extend(Graph& graph, bool warn_on_duplicates = false);

    // modify ids of the second graph to ensure we don't have conflicts
    // then attach tails of this graph to the heads of the other, and extend(g)
    void append(VG& g);

    // don't append or join the nodes in the graphs
    // just ensure that ids are unique, then apply extend
    void combine(VG& g);

    // edit the graph to include the path
    void include(const Path& path);
    // Edit the graph to include all the sequence and edges added by the given
    // paths. Can handle paths that visit nodes in any orientation.
    void edit(const vector<Path>& paths);

    // Find all the points at which a Path enters or leaves nodes in the graph. Adds
    // them to the given map by node ID of sets of bases in the node that will need
    // to become the starts of new nodes.
    void find_breakpoints(const Path& path, map<id_t, set<pos_t>>& breakpoints);

    // Take a map from node ID to a set of offsets at which new nodes should
    // start (which may include 0 and 1-past-the-end, which should be ignored),
    // break the specified nodes at those positions. Returns a map from old node
    // ID to a map from old node start position to new node pointer in the
    // graph. Note that the caller will have to crear and rebuild path rank
    // data.
    map<pos_t, Node*> ensure_breakpoints(const map<id_t, set<pos_t>>& breakpoints);

    // flips the breakpoints onto the forward strand

   // Given a path on nodes that may or may not exist, and a map from node ID
    map<id_t, set<pos_t>> forwardize_breakpoints(const map<id_t, set<pos_t>>& breakpoints);

    // Given a path on nodes that may or may not exist, and a map from node ID
    // in the path's node ID space to a table of offset and actual node, add in
    // all the new sequence and edges required by the path. The given path must
    // not contain adjacent perfect match edits in the same mapping (the removal
    // of which can be accomplished with the simplify() function).
    void add_nodes_and_edges(const Path& path, const map<pos_t, Node*>& node_translation);

    // Add in the given node, by value
    void add_node(const Node& node);
    void add_nodes(const vector<Node>& nodes);
    void add_edge(const Edge& edge);
    void add_edges(const vector<Edge>& edges);
    void add_nodes(const set<Node*>& nodes);
    void add_edges(const set<Edge*>& edges);

    id_t node_count(void);
    id_t edge_count(void);
    id_t total_length_of_nodes(void);
    // returns the rank of the node in the protobuf array that backs the graph
    int node_rank(Node* node);
    int node_rank(id_t id);
    // Number of edges attached to the start of a node
    int start_degree(Node* node);
    // Number of edges attached to the end of a node
    int end_degree(Node* node);
    // Number of edges attached to the left side of a NodeTraversal
    int left_degree(NodeTraversal node);
    // Number of edges attached to the right side of a NodeTraversal
    int right_degree(NodeTraversal node);
    // Get the edges of the specified node, and add them to the given vector.
    // Guaranteed to add each edge only once per call.
    void edges_of_node(Node* node, vector<Edge*>& edges);
    vector<Edge*> edges_of(Node* node);
    vector<Edge*> edges_from(Node* node);
    vector<Edge*> edges_to(Node* node);
    // Get the edges of the specified set of nodes, and add them to the given set of edge pointers.
    void edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges);

    // Sides on the other side of edges to this side of the node
    set<NodeSide> sides_to(NodeSide side);
    // Sides on the other side of edges from this side of the node
    set<NodeSide> sides_from(NodeSide side);
    // Sides from both sides of the node
    set<NodeSide> sides_from(id_t id);
    // Sides to both sides of the node
    set<NodeSide> sides_to(id_t id);
    // union of sides_to and sides_from
    set<NodeSide> sides_of(NodeSide side);
    // All sides connecting to this node
    set<pair<NodeSide, bool>> sides_context(id_t node_id);
    // Use sides_from an sides_to to determine if both nodes have the same context
    bool same_context(id_t id1, id_t id2);
    // determine if the node is an ancestor of this one by trying to find it in a given number of steps
    bool is_ancestor_prev(id_t node_id, id_t candidate_id);
    bool is_ancestor_prev(id_t node_id, id_t candidate_id, set<id_t>& seen, size_t steps = 64);
    // the same but in the other direction
    bool is_ancestor_next(id_t node_id, id_t candidate_id);
    bool is_ancestor_next(id_t node_id, id_t candidate_id, set<id_t>& seen, size_t steps = 64);
    // try to find a common ancestor by walking back up to steps from the first node
    id_t common_ancestor_prev(id_t id1, id_t id2, size_t steps = 64);
    // try to find a common ancestor by walking forward up to steps from the first node
    id_t common_ancestor_next(id_t id1, id_t id2, size_t steps = 64);
    // to-siblings are nodes which also have edges to them from the same nodes as this one
    set<NodeTraversal> siblings_to(const NodeTraversal& traversal);
    // from-siblings are nodes which also have edges to them from the same nodes as this one
    set<NodeTraversal> siblings_from(const NodeTraversal& traversal);
    // full to-siblings are nodes traversals which share exactly the same upstream NodeSides
    set<NodeTraversal> full_siblings_to(const NodeTraversal& trav);
    // full from-siblings are nodes traversals which share exactly the same downstream NodeSides
    set<NodeTraversal> full_siblings_from(const NodeTraversal& trav);
    // general siblings of
    set<Node*> siblings_of(Node* node);
    // removes easily-resolvable redundancy in the graph
    void simplify_siblings(void);
    // does so for all provided to-sibling sets
    void simplify_to_siblings(const set<set<NodeTraversal>>& to_sibs);
    // does so for all provided from-sibling sets
    void simplify_from_siblings(const set<set<NodeTraversal>>& from_sibs);
    // removes intransitive sibling sets, such as where (A, B, C) = S1 but C ∊ S2
    set<set<NodeTraversal>> transitive_sibling_sets(const set<set<NodeTraversal>>& sibs);
    // removes sibling sets which don't have identical orientation
    set<set<NodeTraversal>> identically_oriented_sibling_sets(const set<set<NodeTraversal>>& sibs);
    // determines of pos1 occurs directly before pos2
    bool adjacent(const Position& pos1, const Position& pos2);

    // use the VG class to generate ids
    Node* create_node(const string& seq, id_t id = 0);
    // find a particular node
    Node* get_node(id_t id);
    // Get the subgraph of a node and all the edges it is responsible for (i.e.
    // where it has the minimal ID) and add it into the given VG.
    void nonoverlapping_node_context_without_paths(Node* node, VG& g);
    void expand_context(VG& g, size_t steps, bool add_paths = true);

    // destroy the node at the given pointer. This pointer must point to a Node owned by the graph.
    void destroy_node(Node* node);
    // destroy the node with the given ID.
    void destroy_node(id_t id);
    bool has_node(id_t id);
    bool has_node(Node* node);
    bool has_node(const Node& node);
    Node* find_node_by_name_or_add_new(string name);
    void for_each_node(function<void(Node*)> lambda);
    void for_each_node_parallel(function<void(Node*)> lambda);
    // Go through all the nodes in the same connected component as the given node. Ignores relative orientation.
    void for_each_connected_node(Node* node, function<void(Node*)> lambda);
    void dfs(
        // called when node is first encountered
        const function<void(Node*)>& node_begin_fn,
        // called when node goes out of scope
        const function<void(Node*)>& node_end_fn,
        // called when an edge is encountered
        const function<void(Edge*)>& edge_fn,
        // called when an edge forms part of the DFS spanning tree
        const function<void(Edge*)>& tree_fn,
        // called when we meet an edge in the current tree component
        const function<void(Edge*)>& edge_curr_fn,
        // called when we meet an edge in an already-traversed tree component
        const function<void(Edge*)>& edge_cross_fn);
    // specialization of dfs for only handling nodes
    void dfs(const function<void(Node*)>& node_begin_fn,
             const function<void(Node*)>& node_end_fn);

    // is the graph empty?
    bool empty(void);

    // generates a digest of the serialized graph
    const string hash(void);

    // remove nodes with no sequence
    // these are created in some cases during the process of graph construction
    void remove_null_nodes(void);
    // remove a node but connect all of its predecessor and successor nodes with new edges
    void remove_node_forwarding_edges(Node* node);
    // remove null nodes but connect predecessors and successors, preserving structure
    void remove_null_nodes_forwarding_edges(void);
    // remove edges for which one of the nodes is not present
    void remove_orphan_edges(void);
    // removes edges representing an inversion and edges on the reverse complement
    void remove_inverting_edges(void);
    // true if the graph has inversions, false otherwise
    bool has_inverting_edges(void);

    // Keep paths in the given set of path names. Populates kept_names with the names of the paths it actually found to keep.
    // The paths specified may not overlap. Removes all nodes and edges not used by one of the specified paths.
    void keep_paths(set<string>& path_names, set<string>& kept_names);
    void keep_path(string& path_name);

    // path stats
    // starting from offset in the first node, how many edges do we cross?
    // path must be nonempty and longer than the given length. offset is
    // interpreted as relative to the first node in its on-path
    // orientation, and is inclusive.
    int path_edge_count(list<NodeTraversal>& path, int32_t offset, int path_length);
    // At what offset in its last node does the path starting at this offset in its first node end?
    // path must be nonempty and longer than the given length. offset is
    // interpreted as relative to the first node in its on-path
    // orientation, and is inclusive. Returned offset is remaining unused length
    // in the last node touched.
    int path_end_node_offset(list<NodeTraversal>& path, int32_t offset, int path_length);
    // converts the stored paths in this graph to alignments
    const vector<Alignment> paths_as_alignments(void);
    const string path_sequence(const Path& path);

    SB_Input vg_to_sb_input();
    vector<pair<id_t, id_t> > get_superbubbles(SB_Input sbi);
    vector<pair<id_t, id_t> > get_superbubbles();
    // edges
    // If the given edge cannot be created, returns null.
    // If the given edge already exists, returns the existing edge.
    Edge* create_edge(Node* from, Node* to, bool from_start = false, bool to_end = false);
    Edge* create_edge(id_t from, id_t to, bool from_start = false, bool to_end = false);
    // Makes a left-to-right edge from the left NodeTraversal to the right one, respecting orientations.
    Edge* create_edge(NodeTraversal left, NodeTraversal right);
    // Makes an edge connecting the given sides of nodes.
    Edge* create_edge(NodeSide side1, NodeSide side2);

    // This can take sides in any order
    Edge* get_edge(const NodeSide& side1, const NodeSide& side2);
    // This can take sides in any order
    Edge* get_edge(const pair<NodeSide, NodeSide>& sides);
    // This gets the edge connecting the given oriented nodes in the given order.
    Edge* get_edge(const NodeTraversal& left, const NodeTraversal& right);
    // destroy the edge at the given pointer. This pointer must point to an edge owned by the graph.
    void destroy_edge(Edge* edge);
    // destroy the edge between the given sides of nodes. These can be in either order.
    void destroy_edge(const NodeSide& side1, const NodeSide& side2);
    // This can take sides in any order
    void destroy_edge(const pair<NodeSide, NodeSide>& sides);
    // remove an edge from the node side indexes, so it doesn't show up when you
    // ask for the edges connected to the side of a node. Makes the edge
    // untraversable until the indexes are rebuilt.
    void unindex_edge_by_node_sides(const NodeSide& side1, const NodeSide& side2);
    void unindex_edge_by_node_sides(Edge* edge);
    // add an edge to the node side indexes. Doesn't touch the index of edges by
    // node pairs or the graph; those must be updated seperately.
    void index_edge_by_node_sides(Edge* edge);
    // Get the edge between the given node sides, which can be in either order.
    bool has_edge(const NodeSide& side1, const NodeSide& side2);
    // This can take sides in any order
    bool has_edge(const pair<NodeSide, NodeSide>& sides);
    bool has_edge(Edge* edge);
    bool has_edge(const Edge& edge);
    bool has_inverting_edge(Node* n);
    bool has_inverting_edge_from(Node* n);
    bool has_inverting_edge_to(Node* n);
    void for_each_edge(function<void(Edge*)> lambda);
    void for_each_edge_parallel(function<void(Edge*)> lambda);

    // connect node -> nodes
    // Connects from the right side of the first to the left side of the second
    void connect_node_to_nodes(NodeTraversal node, vector<NodeTraversal>& nodes);
    // You can optionally use the start of the first node instead of the end
    void connect_node_to_nodes(Node* node, vector<Node*>& nodes, bool from_start = false);
    // connect nodes -> node
    // Connects from the right side of the first to the left side of the second
    void connect_nodes_to_node(vector<NodeTraversal>& nodes, NodeTraversal node);
    // You can optionally use the end of the second node instead of the start
    void connect_nodes_to_node(vector<Node*>& nodes, Node* node, bool to_end = false);

    // utilities
    // These only work on forward nodes.

    // Divide a node at a given internal position. Inserts the new nodes in the
    // correct paths, but can't update the ranks, so they need to be cleared and
    // re-calculated by the caller.
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    // This version works on a collection of internal positions, in linear time.
    void divide_node(Node* node, vector<int> positions, vector<Node*>& parts);
    // Divide a path at a position. Also invalidates stored rank information.
    void divide_path(map<long, id_t>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    void to_dot(ostream& out,
                vector<Alignment> alignments = {},
                bool show_paths = false,
                bool walk_paths = false,
                bool annotate_paths = false,
                bool show_mappings = false,
                bool simple_mode = false,
                bool invert_edge_ports = false,
                bool color_variants = false,
                int random_seed = 0);


    void to_dot(ostream& out, vector<Alignment> alignments = {}, bool show_paths = false, bool walk_paths = false,
                            bool annotate_paths = false, bool show_mappings = false, bool invert_edge_ports = false, int random_seed = 0, bool color_variants = false);
   void to_gfa(ostream& out);
    void to_turtle(ostream& out, const string& rdf_base_uri, bool precompress);
    bool is_valid(bool check_nodes = true,
                  bool check_edges = true,
                  bool check_paths = true,
                  bool check_orphans = true);

    // topologically orders nodes
    // Makes sure that Nodes appear in the Protobuf Graph object in their topological sort order.
    void sort(void);
    // helper function, not really meant for external use
    void topological_sort(deque<NodeTraversal>& l);
    void swap_nodes(Node* a, Node* b);

    // Use a topological sort to order and orient the nodes, and then flip some
    // nodes around so that they are oriented the way they are in the sort.
    // Populates nodes_flipped with the ids of the nodes that have had their
    // orientations changed. TODO: update the paths that touch nodes that
    // flipped around
    void orient_nodes_forward(set<id_t>& nodes_flipped);

    // for each path assigns edits that describe a total match of the mapping to the node
    void force_path_match(void);
    // for each path, if a mapping has no edits then make it a perfect match against a node
    // (the same as force_path_match, but only for empty mappings)
    void fill_empty_path_mappings(void);

    // Align to the graph. The graph must be acyclic and contain only end-to-start edges.
    // Will modify the graph by re-ordering the nodes.
    Alignment align(const Alignment& alignment,
                    int32_t match = 2,
                    int32_t mismatch = 2,
                    int32_t gap_open = 3,
                    int32_t gap_extension = 1,
                    size_t max_query_graph_ratio = 0);
    Alignment align(const string& sequence,
                    int32_t match = 2,
                    int32_t mismatch = 2,
                    int32_t gap_open = 3,
                    int32_t gap_extension = 1,
                    size_t max_query_graph_ratio = 0);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

    // returns all node-crossing paths with up to length across node boundaries
    // considers each node in forward orientation to produce the kpaths around it
    void for_each_kpath(int k, bool path_only, int edge_max,
                        function<void(NodeTraversal)> handle_prev_maxed,
                        function<void(NodeTraversal)> handle_next_maxed,
                        function<void(list<NodeTraversal>::iterator, list<NodeTraversal>&)> lambda);
    void for_each_kpath_parallel(int k, bool path_only, int edge_max,
                                 function<void(NodeTraversal)> handle_prev_maxed,
                                 function<void(NodeTraversal)> handle_next_maxed,
                                 function<void(list<NodeTraversal>::iterator, list<NodeTraversal>&)> lambda);
    void for_each_kpath(int k, bool path_only, int edge_max,
                        function<void(NodeTraversal)> handle_prev_maxed,
                        function<void(NodeTraversal)> handle_next_maxed,
                        function<void(size_t,Path&)> lambda);
    void for_each_kpath_parallel(int k, bool path_only, int edge_max,
                                 function<void(NodeTraversal)> handle_prev_maxed,
                                 function<void(NodeTraversal)> handle_next_maxed,
                                 function<void(size_t,Path&)> lambda);
    void for_each_kpath_of_node(Node* node, int k, bool path_only, int edge_max,
                                function<void(NodeTraversal)> handle_prev_maxed,
                                function<void(NodeTraversal)> handle_next_maxed,
                                function<void(list<NodeTraversal>::iterator, list<NodeTraversal>&)> lambda);
    void for_each_kpath_of_node(Node* n, int k, bool path_only, int edge_max,
                                function<void(NodeTraversal)> handle_prev_maxed,
                                function<void(NodeTraversal)> handle_next_maxed,
                                function<void(size_t,Path&)> lambda);

    void kpaths(set<list<NodeTraversal> >& paths, int length, bool path_only, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    void kpaths(vector<Path>& paths, int length, bool path_only, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);

    void kpaths_of_node(Node* node, set<list<NodeTraversal> >& paths,
                        int length, bool path_only, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    void kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, bool path_only, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    void kpaths_of_node(id_t node_id, vector<Path>& paths, int length, bool path_only, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    // Given an oriented start node, a length in bp, a maximum number of edges
    // to cross, and a stack of nodes visited so far, fill in the set of paths
    // with all the paths starting at the oriented start node and going left off
    // its end no longer than the specified length, calling maxed_nodes on nodes
    // which can't be visited due to the edge-crossing limit. Produces paths
    // ending with the specified node. TODO: postfix should not be (potentially)
    // copied on every call.
    void prev_kpaths_from_node(NodeTraversal node, int length, bool path_only, int edge_max, bool edge_bounding,
                               list<NodeTraversal> postfix, set<list<NodeTraversal> >& walked_paths,
                               const vector<string>& followed_paths,
                               function<void(NodeTraversal)>& maxed_nodes);
    // Do the same as prec_kpaths_from_node, except going right, producing a path starting with the specified node.
    void next_kpaths_from_node(NodeTraversal node, int length, bool path_only, int edge_max, bool edge_bounding,
                               list<NodeTraversal> prefix, set<list<NodeTraversal> >& walked_paths,
                               const vector<string>& followed_paths,
                               function<void(NodeTraversal)>& maxed_nodes);

    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(id_t from, id_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    // traversal
    // Get the nodes attached to the left side of the given NodeTraversal, in their proper orientations.
    void nodes_prev(NodeTraversal n, vector<NodeTraversal>& nodes);
    vector<NodeTraversal> nodes_prev(NodeTraversal n);
    // (same but using set) traversals before this node on the same strand
    set<NodeTraversal> travs_to(NodeTraversal node);

    // Get the nodes attached to the right side of the given NodeTraversal, in their proper orientations.
    void nodes_next(NodeTraversal n, vector<NodeTraversal>& nodes);
    vector<NodeTraversal> nodes_next(NodeTraversal n);
    // (uses set) traversals after this node on the same strand
    set<NodeTraversal> travs_from(NodeTraversal node);

    // traversals before this node on the same strand
    set<NodeTraversal> travs_of(NodeTraversal node);

    // Count the nodes attached to the left side of the given NodeTraversal
    int node_count_prev(NodeTraversal n);
    // Count the nodes attached to the right side of the given NodeTraversal
    int node_count_next(NodeTraversal n);

    // paths
    Path create_path(const list<NodeTraversal>& nodes);
    Path create_path(const vector<NodeTraversal>& nodes);
    string path_string(const list<NodeTraversal>& nodes);
    // Assumes the path covers the entirety of any nodes visited. Handles backward nodes.
    string path_string(const Path& path);
    void expand_path(const list<NodeTraversal>& path, vector<NodeTraversal>& expanded);
    // Fill in the node_start map with the first index along the path at which each node appears.
    // Caller is responsible for dealing with orientations.
    void node_starts_in_path(const list<NodeTraversal>& path,
                             map<Node*, int>& node_start);
    // true if nodes share all paths and the mappings they share in these paths
    // are adjacent
    bool nodes_are_perfect_path_neighbors(id_t id1, id_t id2);
    // true if the mapping completely covers the node it maps to and is a perfect match
    bool mapping_is_total_match(const Mapping& m);
    // concatenate the mappings for a pair of nodes; handles multiple mappings per path
    map<string, vector<Mapping>> concat_mappings_for_node_pair(id_t id1, id_t id2);
    // for a list of nodes that we want to concatenate
    map<string, vector<Mapping>> concat_mappings_for_nodes(const list<Node*>& nodes);
    // helper function
    map<string, map<int, Mapping>>
        concat_mapping_groups(map<string, map<int, Mapping>>& r1,
                              map<string, map<int, Mapping>>& r2);

    // These versions handle paths in which nodes can be traversed multiple
    // times. Unfortunately since we're throwing non-const iterators around, we
    // can't take the input path as const.
    void expand_path(list<NodeTraversal>& path, vector<list<NodeTraversal>::iterator>& expanded);
    // To get the starts out of the map this produces, you need to dereference
    // the iterator and then get the address of the NodeTraversal (stored in the
    // list) that you are talking about.
    void node_starts_in_path(list<NodeTraversal>& path,
                             map<NodeTraversal*, int>& node_start);

    // kmers
    void for_each_kmer_parallel(int kmer_size,
                                bool path_only,
                                int edge_max,
                                function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                                int stride = 1,
                                bool allow_dups = false,
                                bool allow_negatives = false);
    void for_each_kmer(int kmer_size,
                       bool path_only,
                       int edge_max,
                       function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                       int stride = 1,
                       bool allow_dups = false,
                       bool allow_negatives = false);
    void for_each_kmer_of_node(Node* node,
                               int kmer_size,
                               bool path_only,
                               int edge_max,
                               function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                               int stride = 1,
                               bool allow_dups = false,
                               bool allow_negatives = false);

    // for gcsa2. For the given kmer of the given length starting at the given
    // offset into the given Node along the given path, fill in end_node and
    // end_offset with where the end of the kmer falls (counting from the right
    // side of the NodeTraversal), prev_chars with the characters that preceed
    // it, next_chars with the characters that follow it, prev_ and
    // next_positions with the ((node ID, orientation), offset) pairs of the
    // places you can come from/go next (from the right end of the kmer).
    // Refuses to follow more than edge_max edges. Offsets are in the path
    // orientation.
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

    // Do the GCSA2 kmers for a node. head_node and tail_node must both be non-
    // null, but only one of those nodes actually needs to be in the graph. They
    // will be examined directly to get their representative characters. They
    // also don't need to be actually owned by the graph; they can be copies.
    void gcsa_handle_node_in_graph(Node* node, int kmer_size, bool path_only,
                                   int edge_max, int stride,
                                   bool forward_only,
                                   Node* head_node, Node* tail_node,
                                   function<void(KmerPosition&)> lambda);

    // GCSA kmers are the kmers in the graph with each node
    // existing in both its forward and reverse-complement orientation. Node IDs
    // in the GCSA graph are 2 * original node ID, +1 if the GCSA node
    // represents the reverse complement, and +0 if it does not. Non-reversing
    // edges link the forward copy of the from node to the forward copy of the
    // to node, and similarly for the reverse complement copies, while reversing
    // edges link the forward copy of the from node to the *reverse complement*
    // copy of the to node, and visa versa. This allows us to index both the
    // forward and reverse strands of every node, and to deal with GCSA's lack
    // of support for reversing edges, with the same trick. Note that
    // start_tail_id, if zero, will be replaced with the ID actually used for the
    // start/end node before lambda is ever called.
    void for_each_gcsa_kmer_position_parallel(int kmer_size, bool path_only,
                                              int edge_max, int stride,
                                              bool forward_only,
                                              id_t& head_id, id_t& tail_id,
                                              function<void(KmerPosition&)> lambda);

    void get_gcsa_kmers(int kmer_size, bool path_only,
                        int edge_max, int stride,
                        bool forward_only,
                        const function<void(vector<gcsa::KMer>&, bool)>& handle_kmers,
                        id_t& head_id, id_t& tail_id);

    void write_gcsa_kmers(int kmer_size, bool path_only,
                          int edge_max, int stride,
                          bool forward_only,
                          ostream& out,
                          id_t& head_id, id_t& tail_id);

    // write the kmers to a tmp file with the given base, return the name of the file
    string write_gcsa_kmers_to_tmpfile(int kmer_size,
                                       bool paths_only,
                                       bool forward_only,
                                       id_t& head_id, id_t& tail_id,
                                       size_t doubling_steps = 2,
                                       size_t size_limit = 200,
                                       const string& base_file_name = ".vg-kmers-tmp-");

    // construct the GCSA index for this graph
    void build_gcsa_lcp(gcsa::GCSA*& gcsa,
                        gcsa::LCPArray*& lcp,
                        int kmer_size,
                        bool paths_only,
                        bool forward_only,
                        size_t doubling_steps = 2,
                        size_t size_limit = 200,
                        const string& base_file_name = ".vg-kmers-tmp-");

    // for pruning graph prior to indexing with gcsa2
    // takes all nodes that would introduce paths of > edge_max edge crossings, removes them, and links their neighbors to
    // head_node or tail_node depending on which direction the path extension was stopped
    void prune_complex(int path_length, int edge_max, Node* head_node, Node* tail_node);
    // wraps the graph with heads and tails before doing the prune
    // utility function for preparing for indexing
    void prune_complex_with_head_tail(int path_length, int edge_max);

private:
    // Call the given function on each kmer. If parallel is specified, goes
    // through nodes one per thread. If node is not null, looks only at kmers of
    // that specific node.
    void _for_each_kmer(int kmer_size,
                        bool path_only,
                        int edge_max,
                        function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)> lambda,
                        bool parallel,
                        int stride,
                        bool allow_dups,
                        bool allow_negatives,
                        Node* node = nullptr);


public:

    // reads
    // note that even if either_strand is false, having backward nodes in the
    // graph will result in some reads from the global reverse strand.
    Alignment random_read(size_t read_len, mt19937& rng, id_t min_id, id_t max_id, bool either_strand);

    // subgraphs
    void disjoint_subgraphs(list<VG>& subgraphs);
    // Get the head nodes (nodes with edges only to their right sides). These are required to be oriented forward.
    void head_nodes(vector<Node*>& nodes);
    vector<Node*> head_nodes(void);
    bool is_head_node(id_t id);
    bool is_head_node(Node* node);
    // distance from head of node to beginning of graph, or -1 if limit exceeded
    int32_t distance_to_head(NodeTraversal node, int32_t limit = 1000);
    int32_t distance_to_head(NodeTraversal node, int32_t limit,
                             int32_t dist, set<NodeTraversal>& seen);
    // Get the tail nodes (nodes with edges only to their left sides). These are required to be oriented forward.
    vector<Node*> tail_nodes(void);
    void tail_nodes(vector<Node*>& nodes);
    bool is_tail_node(id_t id);
    bool is_tail_node(Node* node);
    // distance from tail of node to end of graph, or -1 if limit exceeded
    int32_t distance_to_tail(NodeTraversal node, int32_t limit = 1000);
    int32_t distance_to_tail(NodeTraversal node, int32_t limit,
                             int32_t dist, set<NodeTraversal>& seen);
    int32_t distance_to_tail(id_t id, int32_t limit = 1000);
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    // join head nodes of graph to common null node, creating a new single head.
    Node* join_heads(void);
    // or heads and tails to common new single head or tail (optionally from the start/to the end).
    void join_heads(Node* node, bool from_start = false);
    void join_tails(Node* node, bool to_end = false);

    // add singular head and tail null nodes to graph
    void wrap_with_null_nodes(void);

    // Add a start node and an end node, where all existing heads in the graph
    // are connected to the start node, and all existing tails in the graph are
    // connected to the end node. Any connected components in the graph which do
    // not have either are connected to the start at an arbitrary point, and the
    // end node from nodes going to that arbitrary point. If start_node or
    // end_node is null, a new node will be created. Otherwise, the passed node
    // will be used. Note that this visits every node, to make sure it is
    // attached to all connected components. Note that if a graph has, say,
    // heads but no tails, the start node will be attached buut the end node
    // will be free-floating.
    void add_start_end_markers(int length,
                               char start_char, char end_char,
                               Node*& start_node, Node*& end_node,
                               id_t start_id = 0, id_t end_id = 0);

    bool show_progress;
    string progress_message;
    long progress_count;
    long last_progress;
    ProgressBar* progress;
    void create_progress(const string& message, long count);
    void create_progress(long count);
    void update_progress(long i);
    void destroy_progress(void);

    // for managing parallel construction
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

    void init(void); // setup, ensures that gssw == NULL on startup
    // placeholders for empty
    vector<id_t> empty_ids;
    vector<pair<id_t, bool>> empty_edge_ends;

};

} // end namespace vg

#endif

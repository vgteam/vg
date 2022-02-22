#ifndef VG_VG_HPP_INCLUDED

#define VG_VG_HPP_INCLUDED

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
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include "aligner.hpp"
#include "ssw_aligner.hpp"
#include "region.hpp"
#include "path.hpp"
#include "utility.hpp"
#include "alignment.hpp"
#include "mem.hpp"

#include <vg/vg.pb.h>
#include <vg/io/stream.hpp>
#include <vg/io/protobuf_emitter.hpp>
#include "hash_map.hpp"

#include "progressive.hpp"
#include "lru_cache.h"

#include "Variant.h"
#include "Fasta.h"

#include "swap_remove.hpp"

#include "pictographs.hpp"
#include "colors.hpp"

#include "types.hpp"

#include "nodetraversal.hpp"
#include "nodeside.hpp"

#include "handle.hpp"

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
class VG : public Progressive, public MutablePathDeletableHandleGraph {

public:
    
    ////////////////////////////////////////////////////////////////////////////
    // Custom, concatenateable serialization
    ////////////////////////////////////////////////////////////////////////////
    
    /// Write the contents of this graph to an ostream.
    virtual void serialize(ostream& out) const;
    
    /// Sets the contents of this graph to the contents of a serialized graph from
    /// an istream. The serialized graph must be from the same implementation of the
    /// HandleGraph interface as is calling deserialize(). Can only be called by an
    /// empty graph.
    virtual void deserialize(istream& in);

    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    virtual nid_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    
    /// Return the total number of edges in the graph.
    virtual size_t get_edge_count() const;
    
    // We use the default, linear get_total_length 
    
    /// Get the minimum node ID used in the graph, if any are used
    virtual nid_t min_node_id() const;
    /// Get the maximum node ID used in the graph, if any are used
    virtual nid_t max_node_id() const;
    
    /// Efficiently get the number of edges attached to one side of a handle.
    /// Uses the VG graph's internal degree index.
    virtual size_t get_degree(const handle_t& handle, bool go_left) const;
    
    /// Efficiently check for the existence of an edge using VG graph's internal
    /// index of node sides.
    virtual bool has_edge(const handle_t& left, const handle_t& right) const;
    
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    virtual char get_base(const handle_t& handle, size_t index) const;
    
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    virtual string get_subsequence(const handle_t& handle, size_t index, size_t size) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Returns the number of paths stored in the graph
    virtual size_t get_path_count() const;
    
    /// Determine if a path name exists and is legal to get a path handle for.
    virtual bool has_path(const string& path_name) const;
    
    /// Look up the path handle for the given path name
    virtual path_handle_t get_path_handle(const string& path_name) const;
    
    /// Look up the name of a path from a handle to it
    virtual string get_path_name(const path_handle_t& path_handle) const;
    
    /// Look up whether a path is circular
    virtual bool get_is_circular(const path_handle_t& path_handle) const;
    
    /// Returns the number of node steps in the path
    virtual size_t get_step_count(const path_handle_t& path_handle) const;
    
    /// Get a node handle (node ID and orientation) from a handle to an step on a path
    virtual handle_t get_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the path that an step is on
    virtual path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Get a handle to the first step, or in a circular path to an arbitrary step
    /// considered "first". If the path is empty, returns the past-the-last step
    /// returned by path_end.
    virtual step_handle_t path_begin(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position past the end of a path. This position is
    /// return by get_next_step for the final step in a path in a non-circular path.
    /// Note that get_next_step will *NEVER* return this value for a circular path.
    virtual step_handle_t path_end(const path_handle_t& path_handle) const;
    
    /// Get a handle to the last step, which will be an arbitrary step in a circular path that
    /// we consider "last" based on our construction of the path. If the path is empty
    /// then the implementation must return the same value as path_front_end().
    virtual step_handle_t path_back(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position before the beginning of a path. This position is
    /// return by get_previous_step for the first step in a path in a non-circular path.
    /// Note: get_previous_step will *NEVER* return this value for a circular path.
    virtual step_handle_t path_front_end(const path_handle_t& path_handle) const;
    
    /// Returns a handle to the next step on the path. If the given step is the final step
    /// of a non-circular path, returns the past-the-last step that is also returned by
    /// path_end. In a circular path, the "last" step will loop around to the "first" (i.e.
    /// the one returned by path_begin).
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    virtual step_handle_t get_next_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the previous step on the path. If the given step is the first
    /// step of a non-circular path, this method has undefined behavior. In a circular path,
    /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
    /// the "last" step.
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    virtual step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
    /// Returns true if the step is not the last step in a non-circular path.
    virtual bool has_next_step(const step_handle_t& step_handle) const;
    
    /// Returns true if the step is not the first step in a non-circular path.
    virtual bool has_previous_step(const step_handle_t& step_handle) const;
    
    /// Execute a function on each path in the graph
    virtual bool for_each_path_handle_impl(const function<bool(const path_handle_t&)>& iteratee) const;
    
    /// Loop over the steps of a handle in paths.
    virtual bool for_each_step_on_handle_impl(const handle_t& handle,
                                                    const function<bool(const step_handle_t&)>& iteratee) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Mutable handle-based interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Create a new node with the given sequence and return the handle.
    /// The sequence may not be empty.
    virtual handle_t create_handle(const string& sequence);

    /// Create a new node with the given id and sequence, then return the handle.
    /// The sequence may not be empty.
    /// The ID must be strictly greater than 0.
    virtual handle_t create_handle(const string& sequence, const nid_t& id);
    
    /// Remove the node belonging to the given handle and all of its edges.
    /// Destroys any paths in which the node participates.
    /// Invalidates the destroyed handle.
    /// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
    /// May **NOT** be called during parallel for_each_handle iteration.
    /// May **NOT** be called on the node from which edges are being followed during follow_edges.
    /// May **NOT** be called during iteration over paths, if it would destroy a path.
    /// May **NOT** be called during iteration along a path, if it would destroy that path.
    virtual void destroy_handle(const handle_t& handle);
    
    /// Create an edge connecting the given handles in the given order and orientations.
    virtual void create_edge(const handle_t& left, const handle_t& right);
    
    /// Remove the edge connecting the given handles in the given order and orientations.
    virtual void destroy_edge(const handle_t& left, const handle_t& right);
    
    /// Remove all nodes and edges. Does not update any stored paths.
    virtual void clear();
    
    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated.
    virtual void swap_handles(const handle_t& a, const handle_t& b);
    
    /// Alter the node that the given handle corresponds to so the orientation
    /// indicated by the handle becomes the node's local forward orientation.
    /// Rewrites all edges pointing to the node and the node's sequence to
    /// reflect this. Invalidates all handles to the node (including the one
    /// passed). Returns a new, valid handle to the node in its new forward
    /// orientation. Note that it is possible for the node's ID to change.
    virtual handle_t apply_orientation(const handle_t& handle);
    
    /// Split a handle's underlying node at the given offsets in the handle's
    /// orientation. Returns all of the handles to the parts. Other handles to
    /// the node being split may be invalidated. The split pieces stay in the
    /// same local forward orientation as the original node, but the returned
    /// handles come in the order and orientation appropriate for the handle
    /// passed in.
    virtual vector<handle_t> divide_handle(const handle_t& handle, const vector<size_t>& offsets);
    
    /// Adjust the representation of the graph in memory to improve performance.
    /// Optionally, allow the node IDs to be reassigned to further improve
    /// performance.
    /// Note: Ideally, this method is called one time once there is expected to be
    /// few graph modifications in the future.
    virtual void optimize(bool allow_id_reassignment = true);
    
    /// Reorder the graph's internal structure to match that given.
    /// This sets the order that is used for iteration in functions like for_each_handle.
    /// If compact_ids is true, may (but will not necessarily) compact the id space of the graph to match the ordering, from 1->|ordering|.
    /// In other cases, node IDs will be preserved.
    /// This may be a no-op in the case of graph implementations that do not have any mechanism to maintain an ordering.
    /// This may invalidate outstanding handles.
    /// Returns true if node IDs actually were adjusted to match the given order, and false if they remain unchanged.
    virtual bool apply_ordering(const std::vector<handle_t>& order, bool compact_ids = false);
    
    /// No-op function (required by MutableHandleGraph interface)
    virtual void set_id_increment(const nid_t& min_id);
    
    /// Add the given value to all node IDs. Preserves the paths.
    virtual void increment_node_ids(nid_t increment);
    
    /// Reassign all node IDs as specified by the old->new mapping function.
    virtual void reassign_node_ids(const std::function<nid_t(const nid_t&)>& get_new_id);

    ////////////////////////////////////////////////////////////////////////////
    // Mutable path handle interface
    ////////////////////////////////////////////////////////////////////////////

    /// Destroy the given path. Invalidates handles to the path and its node steps.
    virtual void destroy_path(const path_handle_t& path);

    /// Create a path with the given name.
    virtual path_handle_t create_path_handle(const string& name,
                                             bool is_circular = false);
    
    /// Append a visit to a node to the given path
    virtual step_handle_t append_step(const path_handle_t& path, const handle_t& to_append);
    
    /// Append a visit to a node to the given path
    virtual step_handle_t prepend_step(const path_handle_t& path, const handle_t& to_prepend);
    
    
    /// Delete a segment of a path and rewrite it as some other sequence of steps. Returns a pair
    ///  of step_handle_t's that indicate the range of the new segment in the path.
    virtual pair<step_handle_t, step_handle_t> rewrite_segment(const step_handle_t& segment_begin,
                                                               const step_handle_t& segment_end,
                                                               const vector<handle_t>& new_segment);
    
    /// Make a path circular or non-circular. If the path is becoming circular, the
    /// last step is joined to the first step. If the path is becoming linear, the
    /// step considered "last" is unjoined from the step considered "first" according
    /// to the method path_begin.
    virtual void set_circularity(const path_handle_t& path, bool circular);
    
public:
    
    ////////////////////////////////////////////////////////////////////////////
    // Actual backing storage
    ////////////////////////////////////////////////////////////////////////////

    /// Protobuf-based representation.
    // NB: we can't subclass this safely, so it's best as a member
    Graph graph;

    /// Manages paths of the graph.
    /// Initialized by setting paths._paths = graph.paths.
    Paths paths;

    /// Name of the graph.
    string name;

    /// Current id for Node to be added next.
    nid_t current_id;
    // todo
    //nid_t min_id;
    //nid_t max_id;

    /// `Node`s by id.
    hash_map<nid_t, Node*> node_by_id;

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
    hash_map<nid_t, vector<pair<nid_t, bool>>> edges_on_start;
    /// Stores the destinations and backward flags for edges attached to the ends of nodes (whether that node is "from" or "to").
    hash_map<nid_t, vector<pair<nid_t, bool>>> edges_on_end;

    /// Set the edge indexes through this function. Picks up the sides being
    /// connected by the edge automatically, and silently drops the edge if they
    /// are already connected.
    void set_edge(Edge*);
    void print_edges(void);

    // access the edge indexes through these functions
    
    /// Get nodes and backward flags following edges that attach to this node's start.
    vector<pair<nid_t, bool>>& edges_start(Node* node);
    /// Get nodes and backward flags following edges that attach to this node's start.
    vector<pair<nid_t, bool>>& edges_start(nid_t id);
    /// Get nodes and backward flags following edges that attach to this node's end.
    vector<pair<nid_t, bool>>& edges_end(Node* node);
    /// Get nodes and backward flags following edges that attach to this node's end.
    vector<pair<nid_t, bool>>& edges_end(nid_t id);
    
    // properties of the graph
    // TODO: redundant with handle graph methods!
    size_t size(void); ///< Number of nodes
    size_t length(void); ///< Total sequence length

    // Clear everything
    //void clear(void);

    // constructors

    /// Default constructor.
    VG(void);

    /// Construct from Graph objects serialized in a tagged group stream.
    VG(istream& in, bool showp = false, bool warn_on_duplicates = true);
    void from_istream(istream& in, bool showp = false, bool warn_on_duplicates = true);

    /// Construct from an arbitrary source of Graph protobuf messages, with the message source in control of execution.
    /// Takes a function that takes a callback to call with each Graph object to incorporate.
    /// The Graph object may be moved by the VG.
    VG(const function<void(const function<void(Graph&)>&)>& send_graphs, bool showp = false, bool warn_on_duplicates = true);

    /// Construct from a single Protobuf graph. The same as making an empty VG and using extend().
    VG(const Graph& from, bool showp = false, bool warn_on_duplicates = true);

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
    map<nid_t, vcflib::Variant> get_node_nid_to_variant(vcflib::VariantCallFile vfile);

    // This index stores the same info as alt_paths, but allows us to annotate nodes
    // that flank the variant nodes and variant paths containing no nodes (e.g. deletion edges).
    // The SnarlTraversals are named identically to alt_paths (_alt_[a-z0-9]*_[0-9]*).
    map<string, SnarlTraversal> variant_to_traversal;
    
private:
    /// Does the specified node have any self-loops?
    bool is_self_looping(Node* node);
    /// Use the orientation of the first node as the basis.
    Node* merge_nodes(const list<Node*>& nodes);
public:
    
    /// Turn the graph into a dag by copying strongly connected components expand_scc_steps times
    /// and translating the edges in the component to flow through the copies in one direction.
    /// Assumes that all nodes in the graph are articulated on one consistent strand.
    /// Tolerates doubly-reversing edges in the input graph.
    VG dagify(uint32_t expand_scc_steps,
              unordered_map<nid_t, pair<nid_t, bool> >& node_translation,
              size_t target_min_walk_length = 0,
              size_t component_length_max = 0);
    /// Ensure that all traversals up to max_length are represented as a path on one strand or the other
    /// without taking an inverting edge. All inverting edges are converted to non-inverting edges to
    /// reverse complement nodes. If no inverting edges are present, the strandedness of all nodes is
    /// the same as the input graph. If inverting edges are present, node strandedness is arbitrary.
    VG unfold(uint32_t max_length,
              unordered_map<nid_t, pair<nid_t, bool> >& node_translation);
    /// Create the reverse complemented graph with topology preserved. Record translation in provided map.
    VG reverse_complement_graph(unordered_map<nid_t, pair<nid_t, bool>>& node_translation);
    /// Record the translation of this graph into itself in the provided map.
    void identity_translation(unordered_map<nid_t, pair<nid_t, bool>>& node_translation);
    
    /// Assume two node translations, the over is based on the under; merge them.
    unordered_map<nid_t, pair<nid_t, bool> > overlay_node_translations(const unordered_map<nid_t, pair<nid_t, bool> >& over,
                                                                     const unordered_map<nid_t, pair<nid_t, bool> >& under);
    /// Use our topological sort to quickly break cycles in the graph, return the edges which are removed.
    /// Very non-optimal, but fast.
    vector<Edge> break_cycles(void);
    /// Remove pieces of the graph which are not part of any path.
    void remove_non_path(void);
    /// Remove pieces of the graph which are part of some path.
    void remove_path(void);
    /// Get all of the edges that are on any path.
    set<Edge*> get_path_edges(void);
    
    /// Convert edges that are both from_start and to_end to "regular" ones from end to start.
    void flip_doubly_reversed_edges(void);

    /// Build a graph from a Turtle stream.
    void from_turtle(string filename, string baseuri, bool showp = false);

    /// Destructor.
    ~VG(void);

    /// Copy constructor.
    VG(const VG& other);

    /// Move constructor.
    VG(VG&& other) noexcept;

    /// Copy assignment operator.
    VG& operator=(const VG& other);

    /// Move assignment operator.
    VG& operator=(VG&& other) noexcept;

    // TODO: document all these

    void build_indexes(void);
    void build_node_indexes(void);
    void build_edge_indexes(void);
    void build_indexes_no_init_size(void);
    void build_node_indexes_no_init_size(void);
    void build_edge_indexes_no_init_size(void);
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

    /// Send chunked graphs to a function that will write them somewhere.
    /// Used to internally implement saving to many destinations.
    /// Graph will be serialized in internal storage order.
    /// This is NOT const because it synchronizes path ranks before saving!
    /// We allow the emitted Graph to move.
    void serialize_to_function(const function<void(Graph&)>& emit, nid_t chunk_size = 1000);
    /// Write chunked graphs to a ProtobufEmitter that will write them to a stream.
    /// Use when combining multiple VG objects in a stream.
    /// Graph will be serialized in internal storage order.
    void serialize_to_emitter(vg::io::ProtobufEmitter<Graph>& emitter, nid_t chunk_size = 1000);
    /// Write to a stream in chunked graphs. Adds an EOF marker.
    /// Use when this VG will be the only thing in the stream.
    /// Graph will be serialized in internal storage order.
    void serialize_to_ostream(ostream& out, nid_t chunk_size = 1000);
    /// Write the graph to a file, with an EOF marker.
    /// Graph will be serialized in internal storage order.
    void serialize_to_file(const string& file_name, nid_t chunk_size = 1000);

    // can we handle this with merge?
    //void concatenate(VG& g);

    /// Squish the node IDs down into as small a space as possible. Fixes up paths itself.
    void compact_ids(void);
    /// Squish the node IDs down into as small a space as possible. Fixes up paths itself.
    /// Record translation in provided map.
    void compact_ids(hash_map<nid_t, nid_t> & new_id);
    /// Subtract the given value from all the node IDs. Must not create a node with 0 or negative IDs. Invalidates the paths.
    void decrement_node_ids(nid_t decrement);
    /// Change the ID of the node with the first id to the second, new ID not
    /// used by any node. Invalidates any paths containing the node, since they
    /// are not updated.
    void swap_node_id(nid_t node_id, nid_t new_id);
    /// Change the ID of the given node to the second, new ID not used by any
    /// node. Invalidates the paths. Invalidates any paths containing the node,
    /// since they are not updated.
    void swap_node_id(Node* node, nid_t new_id);
    
    /// Topologically sort the graph, and then apply that sort to re- order the nodes in the backing
    /// data structure. The sort is guaranteed to be stable. This sort is well-defined  on graphs that
    /// are not DAGs, but instead of finding a topological sort it does a heuristic sort to minimize
    /// a feedback arc set.
    void sort();
    
    /// Order the backing graph data structure by node ID
    void id_sort();
    
    /// Iteratively add when nodes and edges are novel. Good when there are very
    /// many overlaps. TODO: If you are using this with warn on duplicates on,
    /// and you know there shouldn't be any duplicates, maybe you should use
    /// merge instead.
    /// This version sorts paths on rank after adding in the path mappings from
    /// the other graph.
    void extend(const VG& g, bool warn_on_duplicates = false);
    /// This version does not sort path mappings by rank. In order to preserve
    /// paths, call Paths::sort_by_mapping_rank() and
    /// Paths::rebuild_mapping_aux() after you are done adding in graphs to this
    /// graph.
    void extend(const Graph& graph, bool warn_on_duplicates = false);
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
    /// paths. Can handle paths that visit nodes in any orientation.
    /// If out_translations is given, a vector of Translations will be written to it,
    /// one per node existing after the edit, describing
    /// how each new or conserved node is embedded in the old graph.
    /// Note that this method sorts the graph and rebuilds the path index, so it should
    /// not be called in a loop.
    ///
    /// If update_paths is true, the paths will be modified to reflect their
    /// embedding in the modified graph. If save_paths is true, the paths as
    /// embedded in the graph will be added to the graph's set of paths. If
    /// break_at_ends is true (or save_paths is true), nodes will be broken at
    /// the ends of paths that start/end woth perfect matches, so the paths can
    /// be added to the vg graph's paths object.
    void edit(vector<Path>& paths_to_add,
              vector<Translation>* out_translations = nullptr,
              bool save_paths = false,
              bool update_paths = false,
              bool break_at_ends = false);

    /// Streaming version of above.  Instead of reading a list of paths into memory
    /// all at once, a file stream is opened from the given path and used to go one-by-one.
     /// Instead of an option to
    /// updtate the in-memory list, an optional output path for the paths is used
    ///
    /// todo: duplicate less code between the two versions. 
    void edit(const string& paths_to_add_path,
              vector<Translation>* out_translations = nullptr,
              bool save_paths = false,
              const string& out_gam_path = "",
              bool break_at_ends = false,
              bool remove_softclips = false);
    
    /// %Edit the graph to include all the sequences and edges added by the
    /// given path. Returns a vector of Translations, one per original-node
    /// fragment. Completely novel nodes are not mentioned, and nodes with no
    /// Translations are assumed to be carried through unchanged. Invalidates
    /// the rank-based Paths index. Does not sort the graph. Suitable for
    /// calling in a loop.
    ///
    /// Can attach newly created nodes on the left of the path to the given set
    /// of dangling NodeSides, and populates the set at the end with the
    /// NodeSide corresponding to the end of the path. This mechanism allows
    /// edits that hit the end of a node to be attached to what comes
    /// before/after the node by the caller, as this function doesn't handle
    /// that.
    vector<Translation> edit_fast(const Path& path, set<NodeSide>& dangling, size_t max_node_size = 1024);

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

    
    /// Return the number of nodes in the graph
    size_t node_count(void) const;
    /// Count the number of edges in the graph.
    size_t edge_count(void) const;
    /// Get the rank of the node in the protobuf array that backs the graph.
    int node_rank(Node* node);
    /// Get the rank of the node in the protobuf array that backs the graph.
    int node_rank(nid_t id);
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
    set<NodeSide> sides_from(nid_t id);
    /// Get the sides to both sides of the node.
    set<NodeSide> sides_to(nid_t id);
    /// Union of sides_to and sides_from.
    set<NodeSide> sides_of(NodeSide side);
    /// Get all sides connecting to this node.
    set<pair<NodeSide, bool>> sides_context(nid_t node_id);
    /// Use sides_from an sides_to to determine if both nodes have the same context.
    bool same_context(nid_t id1, nid_t id2);
    /// Determine if the node is a prev ancestor of this one.
    bool is_ancestor_prev(nid_t node_id, nid_t candidate_id);
    /// Determine if the node is a prev ancestor of this one by trying to find it in a given number of steps.
    bool is_ancestor_prev(nid_t node_id, nid_t candidate_id, set<nid_t>& seen, size_t steps = 64);
    /// Determine if the node is a next ancestor of this one.
    bool is_ancestor_next(nid_t node_id, nid_t candidate_id);
    /// Determine if the node is a next ancestor of this one by trying to find it in a given number of steps.
    bool is_ancestor_next(nid_t node_id, nid_t candidate_id, set<nid_t>& seen, size_t steps = 64);
    /// Try to find a common ancestor by walking back up to steps from the first node.
    nid_t common_ancestor_prev(nid_t id1, nid_t id2, size_t steps = 64);
    /// Try to find a common ancestor by walking forward up to steps from the first node
    nid_t common_ancestor_next(nid_t id1, nid_t id2, size_t steps = 64);
    
    /// Determine if pos1 occurs directly before pos2.
    bool adjacent(const Position& pos1, const Position& pos2);

    /// Create a node. Use the VG class to generate ids.
    Node* create_node(const string& seq);
    /// Create a node. Use a specified, nonzero node ID.
    Node* create_node(const string& seq, nid_t id);
    /// Find a particular node.
    Node* get_node(nid_t id);
    const Node* get_node(nid_t id) const;
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
    void destroy_node(nid_t id);
    /// Determine if the graph has a node with the given ID.
    bool has_node(nid_t id) const;
    /// Determine if the graph contains the given node.
    bool has_node(const Node* node) const;
    /// Determine if the graph contains the given node.
    bool has_node(const Node& node) const;
    /// Find a node with the given name, or create a new one if none is found.
    Node* find_node_by_name_or_add_new(string name);
    /// Run the given function on every node.
    void for_each_node(function<void(Node*)> lambda);
    void for_each_node(function<void(const Node*)> lambda) const;
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
        const unordered_set<NodeTraversal>* sinks);

    /// Specialization of dfs for only handling nodes.
    void dfs(const function<void(NodeTraversal)>& node_begin_fn,
             const function<void(NodeTraversal)>& node_end_fn,
             const vector<NodeTraversal>* sources = NULL,
             const unordered_set<NodeTraversal>* sinks = NULL);

    /// Specialization of dfs for only handling nodes + break function.
    void dfs(const function<void(NodeTraversal)>& node_begin_fn,
             const function<void(NodeTraversal)>& node_end_fn,
             const function<bool(void)>& break_fn);

    /// Is the graph empty?
    bool empty(void) const;

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
    void keep_paths(const set<string>& path_names, set<string>& kept_names);
    void keep_path(const string& path_name);

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

    /// Takes in a pathname and the nucleotide position (like from a vcf) and
    /// returns the node id which contains that position.
    nid_t get_node_at_nucleotide(string pathname, int nuc);

    // edges
    /// Create an edge.
    /// If the given edge cannot be created, returns null.
    /// If the given edge already exists, returns the existing edge.
    Edge* create_edge(Node* from, Node* to, bool from_start = false, bool to_end = false);
    /// Create an edge.
    /// If the given edge cannot be created, returns null.
    /// If the given edge already exists, returns the existing edge.
    Edge* create_edge(nid_t from, nid_t to, bool from_start = false, bool to_end = false);
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
    /// Destroy the edge at the given pointer. This pointer must point to an edge owned by the graph.
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
    bool has_edge(const NodeSide& side1, const NodeSide& side2) const;
    /// Determine if the graph has an edge. This can take sides in any order.
    bool has_edge(const pair<NodeSide, NodeSide>& sides) const;
    /// Determine if the graph has an edge. This can take sides in any order.
    bool has_edge(Edge* edge) const;
    /// Determine if the graph has an edge. This can take sides in any order.
    bool has_edge(const Edge& edge) const;
    /// Determine if the graph has an inverting edge on the given node.
    bool has_inverting_edge(Node* n);
    /// Determine if the graph has an inverting edge from the given node.
    bool has_inverting_edge_from(Node* n);
    /// Determine if the graph has an inverting edge to the given node.
    bool has_inverting_edge_to(Node* n);
    /// Run the given function for each edge.
    void for_each_edge(function<void(Edge*)> lambda);
    void for_each_edge(function<void(const Edge*)> lambda) const;
    /// Run the given function for each edge, in parallel.
    void for_each_edge_parallel(function<void(Edge*)> lambda);


    /// Circularize a subgraph / path using the head / tail nodes.
    void circularize(nid_t head, nid_t tail);
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
    void divide_node(Node* node, vector<int>& positions, vector<Node*>& parts);
    /// Divide a path at a position. Also invalidates stored rank information.
    void divide_path(map<long, nid_t>& path, long pos, Node*& left, Node*& right);
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
                bool noseq_mode = false,
                bool invert_edge_ports = false,
                bool color_variants = false,
                bool ultrabubble_labeling = false,
                bool skip_missing_nodes = false,
                bool ascii_labels = false,
                int random_seed = 0);

    /// Convert the graph to Turtle format.
    void to_turtle(ostream& out, const string& rdf_base_uri, bool precompress);
    /// Determine if the graph is valid or not, according to the specified criteria.
    bool is_valid(bool check_nodes = true,
                  bool check_edges = true,
                  bool check_paths = true,
                  bool check_orphans = true);

    /// Swap the given nodes. TODO: what does that mean?
    void swap_nodes(Node* a, Node* b);
    
    /// Align without base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const string& sequence,
                    const Aligner* aligner,
                    bool traceback = true,
                    bool acyclic_and_sorted = false,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    bool banded_global = false,
                    size_t band_padding_override = 0,
                    size_t max_span = 0,
                    size_t unroll_length = 0,
                    int xdrop_alignment = 0);                // 1 for forward, >1 for reverse, see constructor for the X-drop threshold
    
    /// Align without base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const Alignment& alignment,
                    const Aligner* aligner,
                    const vector<MaximalExactMatch>& mems,
                    bool traceback = true,
                    bool acyclic_and_sorted = false,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    bool banded_global = false,
                    size_t band_padding_override = 0,
                    size_t max_span = 0,
                    size_t unroll_length = 0,
                    int xdrop_alignment = 0);
    Alignment align(const Alignment& alignment,
                    const Aligner* aligner,
                    bool traceback = true,
                    bool acyclic_and_sorted = false,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    bool banded_global = false,
                    size_t band_padding_override = 0,
                    size_t max_span = 0,
                    size_t unroll_length = 0,
                    int xdrop_alignment = 0);
    
    /// Align with default Aligner.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const Alignment& alignment,
                    bool traceback = true,
                    bool acyclic_and_sorted = false,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    bool banded_global = false,
                    size_t band_padding_override = 0,
                    size_t max_span = 0,
                    size_t unroll_length = 0,
                    int xdrop_alignment = 0);                // 1 for forward, >1 for reverse
    /// Align with default Aligner.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align(const string& sequence,
                    bool traceback = true,
                    bool acyclic_and_sorted = false,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    bool banded_global = false,
                    size_t band_padding_override = 0,
                    size_t max_span = 0,
                    size_t unroll_length = 0,
                    int xdrop_alignment = 0);                // 1 for forward, >1 for reverse
    
    /// Align with base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align_qual_adjusted(const Alignment& alignment,
                                  const QualAdjAligner* qual_adj_aligner,
                                  const vector<MaximalExactMatch>& mems,
                                  bool traceback = true,
                                  bool acyclic_and_sorted = false,
                                  size_t max_query_graph_ratio = 0,
                                  bool pinned_alignment = false,
                                  bool pin_left = false,
                                  bool banded_global = false,
                                  size_t band_padding_override = 0,
                                  size_t max_span = 0,
                                  size_t unroll_length = 0,
                                  int xdrop_alignment = 0);              // 1 for forward, >1 for reverse
    Alignment align_qual_adjusted(const Alignment& alignment,
                                  const QualAdjAligner* qual_adj_aligner,
                                  bool traceback = true,
                                  bool acyclic_and_sorted = false,
                                  size_t max_query_graph_ratio = 0,
                                  bool pinned_alignment = false,
                                  bool pin_left = false,
                                  bool banded_global = false,
                                  size_t band_padding_override = 0,
                                  size_t max_span = 0,
                                  size_t unroll_length = 0,
                                  int xdrop_alignment = 0);              // 1 for forward, >1 for reverse
    /// Align with base quality adjusted scores.
    /// Align to the graph.
    /// May modify the graph by re-ordering the nodes.
    /// May add nodes to the graph, but cleans them up afterward.
    Alignment align_qual_adjusted(const string& sequence,
                                  const QualAdjAligner* qual_adj_aligner,
                                  bool traceback = true,
                                  bool acyclic_and_sorted = false,
                                  size_t max_query_graph_ratio = 0,
                                  bool pinned_alignment = false,
                                  bool pin_left = false,
                                  bool banded_global = false,
                                  size_t band_padding_override = 0,
                                  size_t max_span = 0,
                                  size_t unroll_length = 0,
                                  int xdrop_alignment = 0);              // 1 for forward, >1 for reverse
    
    
    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(nid_t from, nid_t to, vector<Path>& paths);
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
    /// Expand a path. TODO: what does that mean?
    void expand_path(const list<NodeTraversal>& path, vector<NodeTraversal>& expanded);
    /// Fill in the node_start map with the first index along the path at which each node appears.
    /// Caller is responsible for dealing with orientations.
    void node_starts_in_path(const list<NodeTraversal>& path,
                             map<Node*, int>& node_start);
    /// Return true if the mapping completely covers the node it maps to and is a perfect match.
    bool mapping_is_total_match(const Mapping& m);
    /// Concatenate the mappings for a pair of nodes; handles multiple mappings per path.
    map<string, vector<mapping_t>> concat_mappings_for_node_pair(nid_t id1, nid_t id2);

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
    /// band_padding_override gives the band padding to use for banded global
    /// alignment. In banded global mode, if the band padding override is
    /// nonzero, permissive banding is not used, and instead the given band
    /// padding is provided. If the band padding override is not provided, the
    /// max span is used as the band padding and permissive banding is enabled.
    Alignment align(const Alignment& alignment,
                    const Aligner* aligner,
                    const QualAdjAligner* qual_adj_aligner,
                    const vector<MaximalExactMatch>& mems,
                    bool traceback = true,
                    bool acyclic_and_sorted = false,
                    size_t max_query_graph_ratio = 0,
                    bool pinned_alignment = false,
                    bool pin_left = false,
                    bool banded_global = false,
                    size_t band_padding_override = 0,
                    size_t max_span = 0,
                    size_t unroll_length = 0,
                    int xdrop_alignment = 0);                // 1 for forward, >1 for reverse


public:

    /// Generate random reads.
    /// Note that even if either_strand is false, having backward nodes in the
    /// graph will result in some reads from the global reverse strand.
    Alignment random_read(size_t read_len, mt19937& rng, nid_t min_id, nid_t max_id, bool either_strand);

    /// Get the head nodes (nodes with edges only to their right sides). These are required to be oriented forward.
    void head_nodes(vector<Node*>& nodes);
    /// Get the head nodes (nodes with edges only to their right sides). These are required to be oriented forward.
    vector<Node*> head_nodes(void);
    /// Determine if a node is a head node.
    bool is_head_node(nid_t id);
    /// Determine if a node is a head node.
    bool is_head_node(Node* node);
    /// Get the tail nodes (nodes with edges only to their left sides). These are required to be oriented forward.
    vector<Node*> tail_nodes(void);
    /// Get the tail nodes (nodes with edges only to their left sides). These are required to be oriented forward.
    void tail_nodes(vector<Node*>& nodes);
    /// Determine if a node is a tail node.
    bool is_tail_node(nid_t id);
    /// Determine if a node is a tail node.
    bool is_tail_node(Node* node);
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
                               nid_t& start_id, nid_t& end_id);

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
    vector<nid_t> empty_ids;
    /// Placeholder for functions that sometimes need to be passed an empty vector
    vector<pair<nid_t, bool>> empty_edge_ends;

    bool warned_about_rewrites = false;
};



} // end namespace vg

#endif

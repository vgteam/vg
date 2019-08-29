#ifndef VG_XG_HPP_INCLUDED
#define VG_XG_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <map>
#include <queue>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>

#include <vg/vg.pb.h>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "position.hpp"
#include "graph.hpp"
#include "path.hpp"
#include "handle.hpp"
#include "utility.hpp"

// We can have DYNAMIC or SDSL-based gPBWTs
#define MODE_DYNAMIC 1
#define MODE_SDSL 2

#define GPBWT_MODE MODE_SDSL

#if GPBWT_MODE == MODE_DYNAMIC
#include "dynamic.hpp"
#endif

namespace xg {

using namespace std;
using namespace sdsl;
using namespace vg;

class XGPath;
//typedef pair<int64_t, bool> Side;
typedef int64_t id_t; // generic id type
// node sides
typedef int64_t side_t;
id_t side_id(const side_t& side);
bool side_is_end(const side_t& side);
side_t make_side(id_t id, bool is_end);
// node traversals
typedef pair<int64_t, int32_t> trav_t; // meant to encode pos+side or pos+strand
id_t trav_id(const trav_t& trav);
bool trav_is_rev(const trav_t& trav);
int32_t trav_rank(const trav_t& trav);
trav_t make_trav(id_t id, bool is_end, int32_t rank);

/**
 * Thrown when attempting to interpret invalid data as an XG index.
 */
class XGFormatError : public runtime_error {
    // Use the runtime_error constructor
    using runtime_error::runtime_error;
};

/**
 * Provides succinct storage for a graph, its positional paths, and a set of
 * embedded threads.
 */
class XG : public PathPositionHandleGraph, public SerializableHandleGraph, public VectorizableHandleGraph {
public:
    
    ////////////////////////////////////////////////////////////////////////////
    // Here are the ways we can construct XG objects (from graph data or files)
    ////////////////////////////////////////////////////////////////////////////
    
    XG(void) = default;
    ~XG(void);
    
    // Construct an XG index by loading from a stream. Throw an XGFormatError if
    // the stream does not produce a valid XG file.
    XG(istream& in);
    XG(Graph& graph);
    XG(function<void(function<void(Graph&)>)> get_chunks);
    
    // We cannot move, assign, or copy until we add code to point SDSL suppots
    // at the new addresses for their vectors.
    XG(const XG& other) = delete;
    XG(XG&& other) = delete;
    XG& operator=(const XG& other) = delete;
    XG& operator=(XG&& other) = delete;
    
    void from_stream(istream& in, bool validate_graph = false,
        bool print_graph = false);

    void from_graph(Graph& graph, bool validate_graph = false,
        bool print_graph = false);
    // Load the graph by calling a function that calls us back with graph chunks.
    // The function passed in here is responsible for looping.
    void from_callback(function<void(function<void(Graph&)>)> get_chunks,
        bool validate_graph = false, bool print_graph = false);

    /// build the graph from another path handle graph
    void from_path_handle_graph(const PathHandleGraph& graph, bool validate_graph = false,
        bool print_graph = false);
    /// handlegraphified version of from_callback().  same logic but ranks aren't processed
    void from_path_handle_graph_callback(function<void(function<void(const PathHandleGraph&)>)> get_chunks,
        bool validate_graph = false, bool print_graph = false);

    /// Actually build the graph
    /// Note that path_nodes is a map to make the output deterministic in path order.
    void build(vector<pair<id_t, string> >& node_label,
               unordered_map<side_t, vector<side_t> >& from_to,
               unordered_map<side_t, vector<side_t> >& to_from,
               map<string, vector<trav_t> >& path_nodes,
               unordered_set<string>& circular_paths,
               bool validate_graph,
               bool print_graph);
               
    // What's the maximum XG version number we can read with this code?
    const static uint32_t MAX_INPUT_VERSION = 12;
    // What's the version we serialize?
    const static uint32_t OUTPUT_VERSION = 12;
               
    // Load this XG index from a stream. Throw an XGFormatError if the stream
    // does not produce a valid XG file.
    void load(istream& in);
    
    // Alias for load() to match the SerializableHandleGraph interface
    void deserialize(istream& in);
    
    void serialize(std::ostream& out) const;
    // Save this XG index to a stream.
    size_t serialize_and_measure(std::ostream& out,
                                 sdsl::structure_tree_node* v = NULL,
                                 std::string name = "") const;
                     
    
    ////////////////////////////////////////////////////////////////////////////
    // Basic API
    ////////////////////////////////////////////////////////////////////////////
    
    // General public statisitcs
    size_t seq_length = 0;
    size_t node_count = 0;
    size_t edge_count = 0;
    size_t path_count = 0;
    
    const uint64_t* sequence_data(void) const;
    const size_t sequence_bit_size(void) const;
    size_t id_to_rank(int64_t id) const;
    int64_t rank_to_id(size_t rank) const;
    size_t max_node_rank(void) const;
    bool has_node(int64_t id) const;
    /// Get the node ID at the given sequence position. Works in 1-based coordinates.
    nid_t node_at_vector_offset(const size_t& pos) const;
    size_t node_vector_offset(const nid_t& id) const;
    Node node(int64_t id) const; // gets node sequence
    string node_sequence(int64_t id) const;
    size_t node_length(int64_t id) const;
    char pos_char(int64_t id, bool is_rev, size_t off) const; // character at position
    string pos_substr(int64_t id, bool is_rev, size_t off, size_t len = 0) const; // substring in range
    // these provide a way to get an index for each node and edge in the g_iv structure and are used by gPBWT
    size_t node_graph_idx(int64_t id) const;
    size_t edge_graph_idx(const Edge& edge) const;
    size_t edge_index(const edge_t& edge) const;

    size_t get_g_iv_size() const;

    ////////////////////////////////////////////////////////////////////////////
    // Here is the old low-level API that needs to be restated in terms of the 
    // locally traversable graph API and then removed.
    ////////////////////////////////////////////////////////////////////////////

    /// Returns true if the given edge is present in the given orientation, and false otherwise.
    bool has_edge(int64_t id1, bool is_start, int64_t id2, bool is_end) const;
    /// Returns true if the given edge is present in either orientation, and false otherwise.
    bool has_edge(const Edge& edge) const;
    
    vector<Edge> edges_of(int64_t id) const;
    vector<Edge> edges_to(int64_t id) const;
    vector<Edge> edges_from(int64_t id) const;
    vector<Edge> edges_on_start(int64_t id) const;
    vector<Edge> edges_on_end(int64_t id) const;
    int indegree(int64_t id) const;
    int outdegree(int64_t id) const;
    
    /// Get the rank of the edge, or numeric_limits<size_t>.max() if no such edge exists.
    // Given an edge which is in the graph in some orientation, return the edge
    // oriented as it actually appears.
    Edge canonicalize(const Edge& edge) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Here is the new locally traversable graph storage API
    ////////////////////////////////////////////////////////////////////////////
    
    /// use the unified graph storage to return a node subgraph
    Graph node_subgraph_id(int64_t id) const;
    Graph node_subgraph_g(int64_t g) const;
    
    /// provide the graph context up to a given length from the current position; step by nodes
    Graph graph_context_id(const pos_t& pos, int64_t length) const;
    Graph graph_context_g(const pos_t& pos, int64_t length) const;
    
    /// return an edge from the three-part encoding used in the graph vector
    /// Edge type encoding:
    /// 1: end to start
    /// 2: end to end
    /// 3: start to start
    /// 4: start to end
    Edge edge_from_encoding(int64_t from, int64_t to, int type) const;
    void idify_graph(Graph& graph) const;
    
    /// a numerical code for the edge type (based on the two reversal flags)
    int edge_type(bool from_start, bool to_end) const;
    int edge_type(const Edge& edge) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Here is the handle graph API
    ////////////////////////////////////////////////////////////////////////////
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const;
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
    /// continue.
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    /// Get the minimum node ID used in the graph, if any are used
    virtual id_t min_node_id() const;
    /// Get the maximum node ID used in the graph, if any are used
    virtual id_t max_node_id() const;
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    virtual char get_base(const handle_t& handle, size_t index) const;
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    virtual string get_subsequence(const handle_t& handle, size_t index, size_t size) const;
    
    // TODO: There's currently no really good efficient way to implement
    // get_degree; we have to decode each edge to work out what node side it is
    // on. So we use the default implementation.
    
    ////////////////////////
    // Path handle graph API
    ////////////////////////
   
    /// Returns the number of paths stored in the graph
    size_t get_path_count() const;
    /// Determine if a path with a given name exists
    bool has_path(const string& path_name) const;
    /// Look up the path handle for the given path name
    path_handle_t get_path_handle(const string& path_name) const;
    /// Look up the name of a path from a handle to it
    string get_path_name(const path_handle_t& path_handle) const;
    /// Look up whether a path is circular
    bool get_is_circular(const path_handle_t& path_handle) const;
    /// Returns the number of node steps in the path
    size_t get_step_count(const path_handle_t& path_handle) const;
    /// Get a node handle (node ID and orientation) from a handle to a step on a path
    handle_t get_handle_of_step(const step_handle_t& step_handle) const;
    /// Returns a handle to the path that an step is on
    path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    /// Get a handle to the first step, or in a circular path to an arbitrary step
    /// considered "first". If the path is empty, returns the past-the-last step
    /// returned by path_end.
    step_handle_t path_begin(const path_handle_t& path_handle) const;
    /// Get a handle to a fictitious position past the end of a path. This position is
    /// return by get_next_step for the final step in a path in a non-circular path.
    /// Note that get_next_step will *NEVER* return this value for a circular path.
    step_handle_t path_end(const path_handle_t& path_handle) const;
    /// Get a handle to the last step, which will be an arbitrary step in a circular path that
    /// we consider "last" based on our construction of the path. If the path is empty
    /// then the implementation must return the same value as path_front_end().
    step_handle_t path_back(const path_handle_t& path_handle) const;
    /// Get a handle to a fictitious position before the beginning of a path. This position is
    /// return by get_previous_step for the first step in a path in a non-circular path.
    /// Note: get_previous_step will *NEVER* return this value for a circular path.
    step_handle_t path_front_end(const path_handle_t& path_handle) const;
    /// Returns true if the step is not the last step in a non-circular path.
    bool has_next_step(const step_handle_t& step_handle) const;
    /// Returns true if the step is not the first step in a non-circular path.
    bool has_previous_step(const step_handle_t& step_handle) const;
    /// Returns a handle to the next step on the path. If the given step is the final step
    /// of a non-circular path, returns the past-the-last step that is also returned by
    /// path_end. In a circular path, the "last" step will loop around to the "first" (i.e.
    /// the one returned by path_begin).
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    step_handle_t get_next_step(const step_handle_t& step_handle) const;
    /// Returns a handle to the previous step on the path. If the given step is the first
    /// step of a non-circular path, this method has undefined behavior. In a circular path,
    /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
    /// the "last" step.
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
    /// Execute a function on each path in the graph
    bool for_each_path_handle_impl(const function<bool(const path_handle_t&)>& iteratee) const;
    /// Executes a function on each step of a handle in any path.
    bool for_each_step_on_handle_impl(const handle_t& handle, const function<bool(const step_handle_t&)>& iteratee) const;
    
    /// Returns the total length of sequence in the path
    size_t get_path_length(const path_handle_t& path_handle) const;
    
    /// Returns the position along the path of the beginning of this step measured in
    /// bases of sequence. In a circular path, positions start at the step returned by
    /// path_begin().
    size_t get_position_of_step(const step_handle_t& step) const;
    
    /// Returns the step at this position, measured in bases of sequence starting at
    /// the step returned by path_begin(). If the position is past the end of the
    /// path, returns path_end().
    step_handle_t get_step_at_position(const path_handle_t& path, const size_t& position) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Higher-level graph API
    ////////////////////////////////////////////////////////////////////////////
    
    // use_steps flag toggles whether dist refers to steps or length in base pairs
    void neighborhood(int64_t id, size_t dist, Graph& g, bool use_steps = true) const;
    void for_path_range(const string& name, int64_t start, int64_t stop, function<void(int64_t node_id, bool rev)> lambda, bool is_rev = false) const;
    void get_path_range(const string& name, int64_t start, int64_t stop, Graph& g, bool is_rev = false) const;
    // basic method to query regions of the graph
    // add_paths flag allows turning off the (potentially costly, and thread-locking) addition of paths
    // when these are not necessary
    // use_steps flag toggles whether dist refers to steps or length in base pairs
    void expand_context(Graph& g, size_t dist, bool add_paths = true, bool use_steps = true,
                        bool expand_forward = true, bool expand_backward = true,
                        int64_t until_node = 0) const;

    // expand by steps (original and default)
    void expand_context_by_steps(Graph& g, size_t steps, bool add_paths = true,
                                 bool expand_forward = true, bool expand_backward = true,
                                 int64_t until_node = 0) const;
    // expand by length
    void expand_context_by_length(Graph& g, size_t length, bool add_paths = true,
                                  bool expand_forward = true, bool expand_backward = true,
                                  int64_t until_node = 0) const;
    // get the nodes one step from the graph
    void get_connected_nodes(Graph& g) const;
    void get_id_range(int64_t id1, int64_t id2, Graph& g) const;
    // walk forward in id space, collecting nodes, until at least length bases covered
    // (or end of graph reached).  if forward is false, go backwards...
    void get_id_range_by_length(int64_t id1, int64_t length, Graph& g, bool forward) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Here is the paths API
    ////////////////////////////////////////////////////////////////////////////

    /// Pull out the path with the given name.
    Path path(const string& name) const;
    /// Get the path string
    string path_string(const Path& path);
    /// Get the path as an alignment
    Alignment path_as_alignment(const string& name);
    /// Convert the Path to an alignment
    Alignment path_as_alignment(const Path& path);
    /// Get all the paths as alignments
    vector<Alignment> paths_as_alignments(void);
    /// Get the path object by name
    const XGPath& get_path(const string& name) const;
    /// Returns the rank of the path with the given name, or 0 if no such path exists.
    size_t path_rank(const string& name) const;
    /// Returns the ranks of the paths prefixed by the given string
    vector<size_t> path_ranks_by_prefix(const string& prefix) const;
    /// Returns the names of the paths prefixed by the given string
    vector<string> path_names_by_prefix(const string& prefix) const;
    /// Returns the paths with the given prefix
    vector<Path> paths_by_prefix(const string& prefix) const;
    /// Returns the maxiumum rank of any existing path. A path does exist at this rank.
    size_t max_path_rank(void) const;
    /// Get the name of the path at the given rank. Ranks begin at 1.
    string path_name(size_t rank) const;
    vector<size_t> paths_of_node(int64_t id) const;
    map<string, vector<Mapping>> node_mappings(int64_t id) const;
    bool path_contains_node(const string& name, int64_t id) const;
    void add_paths_to_graph(map<int64_t, Node*>& nodes, Graph& g) const;
    size_t node_occs_in_path(int64_t id, const string& name) const;
    size_t node_occs_in_path(int64_t id, size_t rank) const;
    vector<size_t> node_ranks_in_path(int64_t id, const string& name) const;
    vector<size_t> node_ranks_in_path(int64_t id, size_t rank) const;
    /// Get the positions (but not the orientations) of the given node on the given path.
    // See also: oriented_occurrences_on_path, which gives rank and orientation
    vector<size_t> position_in_path(int64_t id, const string& name) const;
    vector<size_t> position_in_path(int64_t id, size_t rank) const;
    map<string, vector<size_t> > position_in_paths(int64_t id, bool is_rev = false, size_t offset = 0) const;
    
    map<string, vector<size_t> > distance_in_paths(int64_t id1, bool is_rev1, size_t offset1,
                                                   int64_t id2, bool is_rev2, size_t offset2) const;
    int64_t min_distance_in_paths(int64_t id1, bool is_rev1, size_t offset1,
                                  int64_t id2, bool is_rev2, size_t offset2) const;
    /// Get the ID of the node that covers the given 0-based position along the path.
    int64_t node_at_path_position(const string& name, size_t pos) const;
    /// Get the 0-based start position in the path that covers the given 0-based position along the path.
    size_t node_start_at_path_position(const string& name, size_t pos) const;
    /// Get the graph position at the given 0-based path position
    pos_t graph_pos_at_path_position(const string& name, size_t pos) const;
    size_t path_length(const string& name) const;
    size_t path_length(size_t rank) const;
    /// Return true if the path with the given name is circular, and false if it is not. The path must exist.
    bool path_is_circular(const string& name) const;
    /// Return true if the path with the given rank is circular, and false if it is not. The path must exist.
    bool path_is_circular(size_t rank) const;
    
    /// returns true if the paths are on the same connected component of the graph (constant time)
    bool paths_on_same_component(size_t path_rank_1, size_t path_rank_2) const;
    
    /// returns the ranks (NOT base positions) and orientations of a given node on a path
    vector<pair<size_t, bool>> oriented_occurrences_on_path(int64_t id, size_t path) const;
    
    /// returns the ranks (NOT base positions) and orientations of a given node on a set of paths
    vector<pair<size_t, vector<pair<size_t, bool>>>> oriented_occurrences_on_paths(int64_t id, vector<size_t>& paths) const;
    
    /// returns all of the paths that a node traversal occurs on, the rank of these occurrences on the path
    /// and the orientation of the occurrences. false indicates that the traversal occurs in the same
    /// orientation as in the path, true indicates.
    vector<pair<size_t, vector<pair<size_t, bool>>>> oriented_paths_of_node(int64_t id) const;
    
private:

    ////////////////////////////////////////////////////////////////////////////
    // Here is the New Way (locally traversable graph storage)
    // Everything should be rewritten in terms of these members
    ////////////////////////////////////////////////////////////////////////////

    /// locally traversable graph storage
    /// 
    /// Encoding designed for efficient compression, cache locality, and relativistic traversal of the graph.
    ///
    /// node := { header, edges_to, edges_from }
    /// header := { node_id, node_start, node_length, edges_to_count, edges_from_count }
    /// node_id := integer
    /// node_start := integer (offset in s_iv)
    /// node_length := integer
    /// edges_to_count := integer
    /// edges_from_count := integer
    /// edges_to := { edge_to, ... }
    /// edges_from := { edge_from, ... }
    /// edge_to := { offset_to_previous_node, edge_type }
    /// edge_from := { offset_to_next_node, edge_type }
    int_vector<> g_iv;
    /// delimit node records to allow lookup of nodes in g_civ by rank
    bit_vector g_bv;
    rank_support_v<1> g_bv_rank;
    bit_vector::select_1_type g_bv_select;
    
    // Let's define some offset ints
    const static int G_NODE_ID_OFFSET = 0;
    const static int G_NODE_SEQ_START_OFFSET = 1;
    const static int G_NODE_LENGTH_OFFSET = 2;
    const static int G_NODE_TO_COUNT_OFFSET = 3;
    const static int G_NODE_FROM_COUNT_OFFSET = 4;
    const static int G_NODE_HEADER_LENGTH = 5;
    
    const static int G_EDGE_OFFSET_OFFSET = 0;
    const static int G_EDGE_TYPE_OFFSET = 1;
    const static int G_EDGE_LENGTH = 2;
    
    /// This is a utility function for the edge exploration. It says whether we
    /// want to visit an edge depending on its type, whether we're the to or
    /// from node, whether we want to look left or right, and whether we're
    /// forward or reverse on the node.
    bool edge_filter(int type, bool is_to, bool want_left, bool is_reverse) const;
    
    // This loops over the given number of edge records for the given g node,
    // starting at the given start g vector position. For all the edges that are
    // wanted by edge_filter given the is_to, want_left, and is_reverse flags,
    // the iteratee is called. Returns true if the iteratee never returns false,
    // or false (and stops iteration) as soon as the iteratee returns false.
    bool do_edges(const size_t& g, const size_t& start, const size_t& count,
        bool is_to, bool want_left, bool is_reverse, const function<bool(const handle_t&)>& iteratee) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Here are the bits we need to keep around to talk about the sequence
    ////////////////////////////////////////////////////////////////////////////
    
    // sequence/integer vector
    int_vector<> s_iv;
    // node starts in sequence, provides id schema
    // rank_1(i) = id
    // select_1(id) = i
    bit_vector s_bv; // node positions in siv
    rank_support_v<1> s_bv_rank;
    bit_vector::select_1_type s_bv_select;
    
    ////////////////////////////////////////////////////////////////////////////
    // And here are the bits for tracking actual node IDs
    ////////////////////////////////////////////////////////////////////////////

    // maintain old ids from input, ranked as in s_iv and s_bv
    int64_t min_id = 0; // id ranges don't have to start at 0
    int64_t max_id = 0;
    int_vector<> r_iv; // ids-id_min is the rank

    ////////////////////////////////////////////////////////////////////////////
    // Here is path storage
    ////////////////////////////////////////////////////////////////////////////

    // paths: serialized as bitvectors over nodes and edges
    int_vector<> pn_iv; // path names
    csa_wt<> pn_csa; // path name compressed suffix array
    bit_vector pn_bv;  // path name starts in uncompressed version of csa
    rank_support_v<1> pn_bv_rank;
    bit_vector::select_1_type pn_bv_select;
    int_vector<> pi_iv; // path ids by rank in the path names

    // probably these should get compressed, for when we have whole genomes with many chromosomes
    // the growth in required memory is quadratic but the stored matrix is sparse
    vector<XGPath*> paths; // path entity membership

    // node->path membership
    int_vector<> np_iv;
    bit_vector np_bv; // entity delimiters in ep_iv
    rank_support_v<1> np_bv_rank;
    bit_vector::select_1_type np_bv_select;
    
    char start_marker = '#';
    char end_marker = '$';
};

class XGPath {
public:
    XGPath(void) = default;
    ~XGPath(void) = default;
    // Path name is required here only for complaining intelligently when
    // something goes wrong. We can also spit out the total unique members,
    // because in here is the most efficient place to count them.
    XGPath(const string& path_name,
           const vector<trav_t>& path,
           bool is_circular,
           size_t node_count,
           XG& graph,
           size_t* unique_member_count_out = nullptr);
    // Path names are stored in the XG object, in a compressed fashion, and are
    // not duplicated here.
    
    // These contain rank and select supports and so cannot move or be copied
    // without code to update them.
    XGPath(const XGPath& other) = delete;
    XGPath(XGPath&& other) = delete;
    XGPath& operator=(const XGPath& other) = delete;
    XGPath& operator=(XGPath&& other) = delete;
    int64_t min_node_id = 0;
    wt_gmr<> ids;
    sd_vector<> directions; // forward or backward through nodes
    int_vector<> positions;
    int_vector<> ranks;
    bit_vector offsets;
    rank_support_v<1> offsets_rank;
    bit_vector::select_1_type offsets_select;
    bool is_circular = false;
    void load(istream& in, uint32_t file_version, const function<int64_t(size_t)>& rank_to_id);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* v = NULL,
                     std::string name = "") const;
    // Get a mapping. Note that the mapping will not have its lengths filled in.
    Mapping mapping(size_t offset, const function<int64_t(id_t)>& node_length) const;

    // Get the node orientation at a 0-based offset.
    id_t node(size_t offset) const;
    bool is_reverse(size_t offset) const;
    id_t local_id(id_t id) const;
    id_t external_id(id_t id) const;
    id_t node_at_position(size_t pos) const;
    size_t offset_at_position(size_t pos) const;
};


Mapping new_mapping(const string& name, int64_t id, size_t rank, bool is_reverse);
void to_text(ostream& out, Graph& graph);

// Determine if two edges are equivalent (the same or one is the reverse of the other)
bool edges_equivalent(const Edge& e1, const Edge& e2);

// Given two equivalent edges, return false if they run in the same direction,
// and true if they are articulated in opposite directions.
bool relative_orientation(const Edge& e1, const Edge& e2);

// Return true if we can only arrive at the start of the given oriented node by
// traversing the given edge in reverse, and false if we can do it by traversing
// the edge forwards. (For single-side self loops, this always beans false.) The
// edge must actually attach to the start of the given oriented node.
bool arrive_by_reverse(const Edge& e, int64_t node_id, bool node_is_reverse);

// Like above, but returns true if we can only ever depart via the edge from the
// node by going in reverse over the edge. Also always false for reversing self
// loops.
bool depart_by_reverse(const Edge& e, int64_t node_id, bool node_is_reverse);

// Make an edge from its fields (generally for comparison)
Edge make_edge(int64_t from, bool from_start, int64_t to, bool to_end);

// Position parsing helpers for CLI
void extract_pos(const string& pos_str, int64_t& id, bool& is_rev, size_t& off);
void extract_pos_substr(const string& pos_str, int64_t& id, bool& is_rev, size_t& off, size_t& len);

}

#endif

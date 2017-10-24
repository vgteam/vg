#ifndef VG_XG_HPP_INCLUDED
#define VG_XG_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <map>
#include <queue>
#include <omp.h>
#include <unordered_map>
#include "cpp/vg.pb.h"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "hash_map_set.hpp"
#include "position.hpp"
#include "graph.hpp"
#include "path.hpp"
#include "handle.hpp"

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
class XG : public HandleGraph {
public:
    
    ////////////////////////////////////////////////////////////////////////////
    // Here are the ways we can construct XG objects (from graph data or files)
    ////////////////////////////////////////////////////////////////////////////
    
    XG(void) : seq_length(0),
               node_count(0),
               edge_count(0),
               path_count(0),
               start_marker('#'),
               end_marker('$') { }
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
        bool print_graph = false, bool store_threads = false,
        bool is_sorted_dag = false);
    void from_graph(Graph& graph, bool validate_graph = false,
        bool print_graph = false, bool store_threads = false,
        bool is_sorted_dag = false);
    // Load the graph by calling a function that calls us back with graph chunks.
    // The function passed in here is responsible for looping.
    // If is_sorted_dag is true and store_threads is true, we store the threads
    // with an algorithm that only works on topologically sorted DAGs, but which
    // is faster.
    void from_callback(function<void(function<void(Graph&)>)> get_chunks,
        bool validate_graph = false, bool print_graph = false,
        bool store_threads = false, bool is_sorted_dag = false); 
    void build(vector<pair<id_t, string> >& node_label,
               unordered_map<side_t, vector<side_t> >& from_to,
               unordered_map<side_t, vector<side_t> >& to_from,
               map<string, vector<trav_t> >& path_nodes,
               bool validate_graph,
               bool print_graph,
               bool store_threads,
               bool is_sorted_dag);
               
    // What's the maximum XG version number we can read with this code?
    const static uint32_t MAX_INPUT_VERSION = 5;
    // What's the version we serialize?
    const static uint32_t OUTPUT_VERSION = 5;
               
    // Load this XG index from a stream. Throw an XGFormatError if the stream
    // does not produce a valid XG file.
    void load(istream& in);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* v = NULL,
                     std::string name = "");
                     
    
    ////////////////////////////////////////////////////////////////////////////
    // Basic API
    ////////////////////////////////////////////////////////////////////////////
    
    // General public statisitcs
    size_t seq_length;
    size_t node_count;
    size_t edge_count;
    size_t path_count;
    
    const uint64_t* sequence_data(void) const;
    const size_t sequence_bit_size(void) const;
    size_t id_to_rank(int64_t id) const;
    int64_t rank_to_id(size_t rank) const;
    size_t max_node_rank(void) const;
    bool has_node(int64_t id) const;
    /// Get the node ID at the given sequence position. Works in 1-based coordinates.
    int64_t node_at_seq_pos(size_t pos) const;
    size_t node_start(int64_t id) const;
    Node node(int64_t id) const; // gets node sequence
    string node_sequence(int64_t id) const;
    size_t node_length(int64_t id) const;
    char pos_char(int64_t id, bool is_rev, size_t off) const; // character at position
    string pos_substr(int64_t id, bool is_rev, size_t off, size_t len = 0) const; // substring in range
    size_t node_graph_idx(int64_t id) const;
    size_t edge_graph_idx(const Edge& edge) const;
    
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
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse) const;
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
    virtual void follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual void for_each_handle(const function<bool(const handle_t&)>& iteratee) const;
    /// Return the number of nodes in the graph
    virtual size_t node_size() const;

    ////////////////////////////////////////////////////////////////////////////
    // Higher-level graph API
    ////////////////////////////////////////////////////////////////////////////
    
    // use_steps flag toggles whether dist refers to steps or length in base pairs
    void neighborhood(int64_t id, size_t dist, Graph& g, bool use_steps = true) const;
    void for_path_range(const string& name, int64_t start, int64_t stop, function<void(int64_t node_id)> lambda, bool is_rev = false) const;
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

    // Pull out the path with the given name.
    Path path(const string& name) const;
    // Returns the rank of the path with the given name, or 0 if no such path
    // exists.
    size_t path_rank(const string& name) const;
    // Returns the maxiumum rank of any existing path. A path does exist at this
    // rank.
    size_t max_path_rank(void) const;
    // Get the name of the path at the given rank. Ranks begin at 1.
    string path_name(size_t rank) const;
    vector<size_t> paths_of_node(int64_t id) const;
    map<string, vector<Mapping>> node_mappings(int64_t id) const;
    bool path_contains_node(const string& name, int64_t id) const;
    void add_paths_to_graph(map<int64_t, Node*>& nodes, Graph& g) const;
    size_t node_occs_in_path(int64_t id, const string& name) const;
    size_t node_occs_in_path(int64_t id, size_t rank) const;
    vector<size_t> node_ranks_in_path(int64_t id, const string& name) const;
    vector<size_t> node_ranks_in_path(int64_t id, size_t rank) const;
    vector<size_t> position_in_path(int64_t id, const string& name) const;
    vector<size_t> position_in_path(int64_t id, size_t rank) const;
    map<string, vector<size_t> > position_in_paths(int64_t id, bool is_rev = false, size_t offset = 0) const;
    map<string, vector<size_t> > distance_in_paths(int64_t id1, bool is_rev1, size_t offset1,
                                                   int64_t id2, bool is_rev2, size_t offset2) const;
    int64_t min_distance_in_paths(int64_t id1, bool is_rev1, size_t offset1,
                                  int64_t id2, bool is_rev2, size_t offset2) const;
    /// Get the ID of the node that covers the given 0-based position along the path.
    int64_t node_at_path_position(const string& name, size_t pos) const;
    /// Get the Mapping that covers the given 0-based position along the path.
    Mapping mapping_at_path_position(const string& name, size_t pos) const;
    /// Get the 0-based start position in the path that covers the given 0-based position along the path.
    size_t node_start_at_path_position(const string& name, size_t pos) const;
    Alignment target_alignment(const string& name, size_t pos1, size_t pos2, const string& feature) const;
    size_t path_length(const string& name) const;
    size_t path_length(size_t rank) const;
    // nearest node (in steps) that is in a path, and the paths
    pair<int64_t, vector<size_t> > nearest_path_node(int64_t id, int max_steps = 16) const;
    int64_t min_approx_path_distance(int64_t id1, int64_t id2) const;
    
    /// returns all of the paths that a node traversal occurs on, the rank of these occurrences on the path
    /// and the orientation of the occurrences. false indicates that the traversal occurs in the same
    /// orientation as in the path, true indicates.
    vector<pair<size_t, vector<pair<size_t, bool>>>> oriented_paths_of_node(int64_t id) const;
    
    /// sets a pointer to a memoized result from oriented_paths_of_node. if no memo is provided, the result will be
    /// queried and stored in local_paths_var and the pointer will be set to point to this variable so that no
    /// additional code paths are necessary
    void memoized_oriented_paths_of_node(int64_t id, vector<pair<size_t, vector<pair<size_t, bool>>>>& local_paths_var,
                                         vector<pair<size_t, vector<pair<size_t, bool>>>>*& paths_of_node_ptr_out,
                                         unordered_map<id_t, vector<pair<size_t, vector<pair<size_t, bool>>>>>* paths_of_node_memo = nullptr) const;
    
    /// returns a the memoized result from get_handle if a memo is provided that the result has been queried previously
    /// otherwise, returns the result of get_handle directly and stores it in the memo if one is provided
    handle_t memoized_get_handle(int64_t id, bool rev, unordered_map<pair<int64_t, bool>, handle_t>* handle_memo = nullptr) const;
    
    /// the oriented distance (positive if pos2 is further along the path than pos1, otherwise negative)
    /// estimated by the distance along the nearest shared path to the two positions plus the distance
    /// to traverse to that path. returns numeric_limits<int64_t>::max() if no pair of nodes that occur same
    /// on the strand of a path are reachable within the max search distance (measured in sequence length,
    /// not node length).
    int64_t closest_shared_path_oriented_distance(int64_t id1, size_t offset1, bool rev1,
                                                  int64_t id2, size_t offset2, bool rev2,
                                                  size_t max_search_dist = 100,
                                                  unordered_map<id_t, vector<pair<size_t, vector<pair<size_t, bool>>>>>* paths_of_node_memo = nullptr,
                                                  unordered_map<pair<int64_t, bool>, handle_t>* handle_memo = nullptr) const;
    
    /// returns a vector of (node id, is reverse, offset) tuples that are found by jumping a fixed oriented distance
    /// along path(s) from the given position. if the position is not on a path, searches from the position to a path
    /// and adds/subtracts the search distance to the jump depending on the search direction. returns an empty vector
    /// if there is no path within the max search distance or if the jump distance goes past the end of the path
    vector<tuple<int64_t, bool, size_t>> jump_along_closest_path(int64_t id, bool is_rev, size_t offset, int64_t jump_dist, size_t max_search_dist,
                                                                 unordered_map<id_t, vector<pair<size_t, vector<pair<size_t, bool>>>>>* paths_of_node_memo = nullptr,
                                                                 unordered_map<pair<int64_t, bool>, handle_t>* handle_memo = nullptr) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // gPBWT API
    ////////////////////////////////////////////////////////////////////////////
    
#if GPBWT_MODE == MODE_SDSL
    // We keep our strings in instances of this cool run-length-compressed wavelet tree.
    using rank_select_int_vector = sdsl::wt_rlmn<sdsl::sd_vector<>>;
#elif GPBWT_MODE == MODE_DYNAMIC
    using rank_select_int_vector = dyn::rle_str;
#endif

    // We need the w function, which we call the "where_to" function. It tells
    // you, from a given visit at a given side, what visit offset if you go to
    // another side.
    int64_t where_to(int64_t current_side, int64_t visit_offset, int64_t new_side) const;
    
    // This is another version of the where_to function which requires that you
    // supply two vectors of edges.
    // edges_into_new -> edges going into new_side
    // edges_out_of_old -> edges coming out of current_side
    // this is to save the overhead of re-extracting these edge-vectors in cases
    // where you're calling where_to between the same two sides (but with 
    // different offsets) many times. Otherwise use version above
    int64_t where_to(int64_t current_side, int64_t visit_offset, int64_t new_side,
      vector<Edge>& edges_into_new, vector<Edge>& edges_out_of_old) const;
    
    // We define a thread visit that's much smaller than a Protobuf Mapping.
    struct ThreadMapping {
        int64_t node_id;
        bool is_reverse;
        /// We need comparison for deduplication in sets and canonically orienting threads
        bool operator<(const ThreadMapping& other) const {
            return tie(node_id, is_reverse) < tie(other.node_id, other.is_reverse);
        }
    };
    
    // we have a public function for querying the contents of the h_iv vector
    int64_t node_height(ThreadMapping node) const;
    
    // We define a thread as just a vector of these things, instead of a bulky
    // Path.
    using thread_t = vector<ThreadMapping>;
    
    // Count matches to a subthread among embedded threads

    /// Insert a thread. Name must be unique or empty.
    /// bs_bake() and tn_bake() need to be called before queries.
    void insert_thread(const thread_t& t, const string& name);
    /// Insert a whole group of threads. Names should be unique or empty (though
    /// they aren't used yet). The indexed graph must be a DAG, at least in the
    /// subset traversed by the threads. (Reversing edges are fine, but the
    /// threads in a node must all run in the same direction.) This uses a
    /// special efficient batch insert algorithm for DAGs that lets us just scan
    /// the graph and generate nodes' B_s arrays independently. This must be
    /// called only once, and no threads can have been inserted previously.
    /// Otherwise the gPBWT data structures will be left in an inconsistent
    /// state.
    void insert_threads_into_dag(const vector<thread_t>& t, const vector<string>& names);
    /// Read all the threads embedded in the graph.
    map<string, list<thread_t> > extract_threads(bool extract_reverse) const;
    /// Extract a particular thread by name. Name may not be empty.
    thread_t extract_thread(const string& name) const;
    /// Extract a set of threads matching a pattern.
    map<string, list<thread_t> > extract_threads_matching(const string& pattern, bool reverse) const;
    /// Extract a particular thread, referring to it by its offset at node; step
    /// it out to a maximum of max_length
    thread_t extract_thread(xg::XG::ThreadMapping node, int64_t offset, int64_t max_length);
    /// Count matches to a subthread among embedded threads
    size_t count_matches(const thread_t& t) const;
    size_t count_matches(const Path& t) const;
    
    /**
     * Represents the search state for the graph PBWT, so that you can continue
     * a search with more of a thread, or backtrack.
     *
     * By default, represents an un-started search (with no first visited side)
     * that can be extended to the whole collection of visits to a side.
     */
    struct ThreadSearchState {
        // What side have we just arrived at in the search?
        int64_t current_side = 0;
        // What is the first visit at that side that is selected?
        int64_t range_start = 0;
        // And what is the past-the-last visit that is selected?
        int64_t range_end = numeric_limits<int64_t>::max();
        
        // How many visits are selected?
        inline int64_t count() {
            return range_end - range_start;
        }
        
        // Return true if the range has nothing selected.
        inline bool is_empty() {
            return range_end <= range_start;
        }
    };
    
    // Extend a search with the given section of a thread.
    void extend_search(ThreadSearchState& state, const thread_t& t) const;

    /// Extend a search with the given single ThreadMapping.
    void extend_search(ThreadSearchState& state, const ThreadMapping& t) const;
    
    /// Select only the threads (if any) starting with a particular
    /// ThreadMapping, and not those continuing through it.
    ThreadSearchState select_starting(const ThreadMapping& start) const;

    /// Take a node id and side and return the side id
    int64_t id_rev_to_side(int64_t id, bool is_rev) const;

    /// Take a side and give a node id / rev pair
    pair<int64_t, bool> side_to_id_rev(int64_t side) const;

    /// The number of threads starting at this side
    int64_t threads_starting_on_side(int64_t side) const;
    
    /// Given a side and offset, return the id of the thread starting there (or 0 if none)
    int64_t thread_starting_at(int64_t side, int64_t offset) const;

    /// Given a thread id and the reverse state get the starting side and offset
    pair<int64_t, int64_t> thread_start(int64_t thread_id) const;

    /// Given a thread id, return its name
    string thread_name(int64_t thread_id) const;

    /// Gives the thread start for the given thread
    pair<int64_t, int64_t> thread_start(int64_t thread_id, bool is_rev) const;

    /// Gives the thread ids of those whose names start with this pattern
    vector<int64_t> threads_named_starting(const string& pattern) const;
    
    /// Select only the threads (if any) continuing through a particular
    /// ThreadMapping, and not those starting there.
    ThreadSearchState select_continuing(const ThreadMapping& start) const;

    // Dump the whole B_s array to the given output stream as a report.
    void bs_dump(ostream& out) const;
    
    char start_marker;
    char end_marker;
    
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
    /// edge_to := { offset_to_next_node, edge_type }
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
    
    // And some masks
    const static size_t HIGH_BIT = (size_t)1 << 63;
    const static size_t LOW_BITS = 0x7FFFFFFFFFFFFFFF;
    
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
    int_vector<> i_iv;
    int64_t min_id; // id ranges don't have to start at 0
    int64_t max_id;
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

    ////////////////////////////////////////////////////////////////////////////
    // Succinct thread storage (the gPBWT)
    ////////////////////////////////////////////////////////////////////////////
    
    // Threads are haplotype paths in the graph with no edits allowed, starting
    // and stopping at node boundaries.
    
    // TODO: Explain the whole graph PBWT extension here
    
    // Basically we keep usage counts for every element in the graph, and and
    // array of next-node-start sides for each side in the graph. We number
    // sides as 2 * xg internal node ID, +1 if it's a right side. This leaves us
    // 0 and 1 free for representing the null destination side and to use as a
    // per-side array run separator, respectively.
    
    // This holds, for each node and edge, in each direction (with indexes as in
    // the entity vector g_iv, *2, and +1 for reverse), the usage count (i.e.
    // the number of times it is visited by encoded threads). This doesn't have
    // to be dynamic since the length will never change. Remember that entity
    // ranks are 1-based, so if you have an entity rank you have to subtract 1
    // to get its position here. We have to track separately for both directions
    // because, even though when everything is inserted the usage counts in both
    // directions are the same, while we're inserting a thread in one direction
    // and not (yet) the other, the usage counts in both directions will be
    // different.
    int_vector<> h_iv;  // only used in construction
    vlc_vector<> h_civ;

    // This (as an extension to the algorithm described in the paper) holds the
    // number of threads beginning at each node. This isn't any extra
    // information relative to what's in the usage count array, but it's cheaper
    // (probably) to maintain this rather than to scan through all the edges on
    // a side every time.
    // ts stands for "thread start"
    int_vector<> ts_iv;  // only used in construction
    vlc_vector<> ts_civ;
    
#if GPBWT_MODE == MODE_SDSL
    // We use this for creating the sub-parts of the uncompressed B_s arrays.
    // We don't really support rank and select on this.
    vector<string> bs_arrays;
#endif
    
    // This holds the concatenated Benedict arrays, with BS_SEPARATOR separating
    // them, and BS_NULL noting the null side (i.e. the thread ends at this
    // node). Instead of holding destination sides, we actually hold the index
    // of the edge that gets taken to the destination side, out of all edges we
    // could take leaving the node. We offset all the values up by 2, to make
    // room for the null sentinel and the separator. Currently the separator
    // isn't used; we just place these by side.
    rank_select_int_vector bs_single_array;

    // thread name storage
    // CSA that lets us look up names efficiently, build from ordered null-delimited names
    // thread ids are taken to be the rank in the source text for tn_csa
    csa_bitcompressed<> tn_csa;
    // allows us to go from positions in the CSA to thread ids
    // rank(i) gives us our thread index for a position in tn_csa's source
    // select(i) gives us the thread name start for a given thread id
    sd_vector<> tn_cbv;
    sd_vector<>::rank_1_type tn_cbv_rank;
    sd_vector<>::select_1_type tn_cbv_select;
    // allows us to go from thread ids to thread start positions, enabling named queries of the graph
    vlc_vector<> tin_civ; // from thread id to side id / reverse thread id to side id
    vlc_vector<> tio_civ; // from thread id to offset / reverse thread id to offset
    // thread starts ordered by their identifiers so we can map from sides into thread ids
    wt_int<> side_thread_wt;
    
    // Holds the names of threads while they are being inserted, before the
    // succinct name representation is built.
    string names_str;

    // A "destination" is either a local edge number + 2, BS_NULL for stopping,
    // or possibly BS_SEPARATOR for cramming multiple Benedict arrays into one.
    using destination_t = size_t;
    
    // Constants used as sentinels in bs_iv above.
    const static destination_t BS_SEPARATOR;
    const static destination_t BS_NULL;
    
    // We access this only through these wrapper methods, because we're going to
    // swap out functionality.
    // Sides are from 1-based node ranks, so start at 2.
    // Get the item in a B_s array for a side at an offset.
    destination_t bs_get(int64_t side, int64_t offset) const;
    // Get the rank of a position among positions pointing to a certain
    // destination from a side.
    size_t bs_rank(int64_t side, int64_t offset, destination_t value) const;
    // Set the whole B_s array for a size. May throw an error if B_s for that
    // side has already been set (as overwrite is not necessarily possible).
    void bs_set(int64_t side, vector<destination_t> new_array);
    // Insert into the B_s array for a side
    void bs_insert(int64_t side, int64_t offset, destination_t value);
    
    // Prepare the B_s array data structures for query. After you call this, you
    // shouldn't call bset or bs_insert.
    void bs_bake();
    
    // Prepare the succinct thread name representation for queries
    void tn_bake();
};

class XGPath {
public:
    XGPath(void) { }
    ~XGPath(void) { }
    // Path name is required here only for complaining intelligently when
    // something goes wrong. We can also spit out the total unique members,
    // because in here is the most efficient place to count them.
    XGPath(const string& path_name,
           const vector<trav_t>& path,
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
    
    rrr_vector<> nodes;
    rrr_vector<>::rank_1_type nodes_rank;
    rrr_vector<>::select_1_type nodes_select;
    wt_gmr<> ids;
    sd_vector<> directions; // forward or backward through nodes
    int_vector<> positions;
    int_vector<> ranks;
    bit_vector offsets;
    rank_support_v<1> offsets_rank;
    bit_vector::select_1_type offsets_select;
    void load(istream& in);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* v = NULL,
                     std::string name = "") const;
    // Get a mapping. Note that the mapping will not have its lengths filled in.
    Mapping mapping(size_t offset) const; // 0-based
};


Mapping new_mapping(const string& name, int64_t id, size_t rank, bool is_reverse);
void parse_region(const string& target, string& name, int64_t& start, int64_t& end);
void to_text(ostream& out, Graph& graph);

// Serialize a rank_select_int_vector in an SDSL serialization compatible way. Returns the number of bytes written.
size_t serialize(XG::rank_select_int_vector& to_serialize, ostream& out,
    sdsl::structure_tree_node* parent, const std::string name);

// Deserialize a rank_select_int_vector in an SDSL serialization compatible way.
void deserialize(XG::rank_select_int_vector& target, istream& in);

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

// Helpers for when we're picking up parts of the graph without returning full Node objects
char reverse_complement(const char& c);
string reverse_complement(const string& seq);

// Position parsing helpers for CLI
void extract_pos(const string& pos_str, int64_t& id, bool& is_rev, size_t& off);
void extract_pos_substr(const string& pos_str, int64_t& id, bool& is_rev, size_t& off, size_t& len);

}

#endif

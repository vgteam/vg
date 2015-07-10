#ifndef VG_H
#define VG_H

#include <vector>
#include <set>
#include <string>
#include <deque>
#include <list>
#include <omp.h>
#include <unistd.h>
#include <limits.h>
#include <algorithm>
#include <random>

#include "gssw.h"
#include "gssw_aligner.hpp"
#include "region.hpp"
#include "path.hpp"
#include "utility.hpp"

#include "vg.pb.h"
#include "hash_map.hpp"

#include "progress_bar.hpp"
#include "lru_cache.h"

#include "Variant.h"
#include "Fasta.h"

#include "swap_remove.hpp"

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
    
    inline bool operator<(const NodeTraversal& other) const {
        return node < other.node || (node == other.node && backward < other.backward);
    }
};

// Represents a sequence graph. Graphs consist of nodes, connected by edges.
// Cycles are not currently permitted. Nodes carry forward-oriented sequences.
// Edges are directed, with a "from" and to" node, and are generally used to
// connect the end of the "from" node to the start of the "to" node. However,
// edges can connect to either the start or end of either node, in general, as
// long as they do not allow the same node to be visited twice along a path.
// Graphs have "head" and "tail" nodes, which are overall at the left/right of
// the graph, with nothing before/after them. Because otherwise identifying
// these nodes (i.e. classifying a terminal node as a head or a tail) would
// require a topological sort, we require that all head and tail nodes be in the
// same relative orientation. Head nodes must have edges only to their right
// sides, and tail nodes must have edges only to their left sides. There must be
// no possible path in the graph containing two head nodes or two tail nodes.
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
    int64_t current_id;
    // todo
    //int64_t min_id;
    //int64_t max_id;

    // nodes by id
    hash_map<int64_t, Node*> node_by_id;

    // edges by nodes they connect
    // Since cycles and duplicate edges are not permitted, two edges cannot connect the same pair of nodes.
    // Each edge is indexed here with the smaller ID first. The actual node order is recorded in the Edge object.
    pair_hash_map<pair<int64_t, int64_t>, Edge*> edge_by_id;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    hash_map<Node*, int> node_index;

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    hash_map<Edge*, int> edge_index;

    // edges indexed by nodes they connect
    // Stores the destinations and backward flags for edges attached to the starts of nodes (whether that node is "from" or "to").
    hash_map<int64_t, vector<pair<int64_t, bool>>> edges_on_start;
    // Stores the destinations and backward flags for edges attached to the ends of nodes (whether that node is "from" or "to").
    hash_map<int64_t, vector<pair<int64_t, bool>>> edges_on_end;

    // set the edge indexes through this function
    void set_edge(int64_t from, int64_t to, Edge*);
    void print_edges(void);

    // access the edge indexes through these functions
    // Get nodes and backward flags following edges that attach to this node's start
    vector<pair<int64_t, bool>>& edges_start(Node* node);
    vector<pair<int64_t, bool>>& edges_start(int64_t id);
    // Get nodes and backward flags following edges that attach to this node's end
    vector<pair<int64_t, bool>>& edges_end(Node* node);
    vector<pair<int64_t, bool>>& edges_end(int64_t id);
    
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
       int vars_per_region,
       int max_node_size = 0,
       bool showprog = false);
    void from_alleles(const map<long, set<vcflib::VariantAllele> >& altp,
                      string& seq,
                      string& chrom);
    void vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, set<vcflib::VariantAllele> >& altp,
                                int start_pos,
                                int stop_pos,
                                int max_node_size = 0);
    void slice_alleles(map<long, set<vcflib::VariantAllele> >& altp,
                       int start_pos,
                       int stop_pos,
                       int max_node_size);
    void dice_nodes(int max_node_size);

    void from_gfa(istream& in, bool showp = false);


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

    // todo
    // set vg up to not build indexes
    // providing very light runtime
    // this would disable a ton of functions
    // and change the deserialization semantics
    //
    //void no_indexes(void);
    //void yes_indexes(void);

    void build_indexes(void);
    void index_paths(void);
    void clear_indexes(void);
    void clear_indexes_no_resize(void);
    void resize_indexes(void);
    void rebuild_indexes(void);

    // literally merge protobufs
    void merge(Graph& g);
    void merge(VG& g);

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
    void serialize_to_ostream(ostream& out, int64_t chunk_size = 1000);
    void serialize_to_file(const string& file_name, int64_t chunk_size = 1000);

    // can we handle this with merge?
    //void concatenate(VG& g);

    int64_t max_node_id(void);
    int64_t min_node_id(void);
    void compact_ids(void);
    void increment_node_ids(int64_t increment);
    void decrement_node_ids(int64_t decrement);
    void swap_node_id(int64_t node_id, int64_t new_id);
    void swap_node_id(Node* node, int64_t new_id);

    // iteratively add when nodes and edges are novel
    // good when there are very many overlaps
    void extend(VG& g);
    void extend(Graph& graph);

    // modify ids of the second graph to ensure we don't have conflicts
    // then attach tails of this graph to the heads of the other, and extend(g)
    void append(VG& g);

    // don't append or join the nodes in the graphs
    // just ensure that ids are unique, then apply extend
    void combine(VG& g);

    // edit the graph to include the path
    void include(const Path& path);
    // or a set of mappings against one node
    void edit_node(int64_t node_id,
                   const vector<Mapping>& mappings,
                   map<pair<size_t, int64_t>, pair<set<Node*>, set<Node*>>>& cut_trans);
    // for each node, modify it with the associated mappings
    void edit(const map<int64_t, vector<Mapping> >& mappings,
              map<pair<int64_t, size_t>, pair<int64_t, size_t> >& del_f,
              map<pair<int64_t, size_t>, pair<int64_t, size_t> >& del_t);
    void edit(const vector<Path>& paths);
    
    // Add in the given node, by value
    void add_node(Node& node);
    void add_nodes(vector<Node>& nodes);
    void add_edge(Edge& edge);
    void add_edges(vector<Edge>& edges);
    void add_nodes(set<Node*>& nodes);
    void add_edges(set<Edge*>& edges);

    int64_t node_count(void);
    int64_t edge_count(void);
    int64_t total_length_of_nodes(void);
    // Number of edges attached to the start of a node
    int start_degree(Node* node);
    // Number of edges attached to the end of a node
    int end_degree(Node* node);
    // Number of edges attached to the left side of a NodeTraversal
    int left_degree(NodeTraversal node);
    // Number of edges attached to the right side of a NodeTraversal
    int right_degree(NodeTraversal node);
    void edges_of_node(Node* node, vector<Edge*>& edges);
    void edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges);

    // use the VG class to generate ids
    Node* create_node(string seq, int64_t id = 0);
    Node* get_node(int64_t id);
    void node_context(Node* node, VG& g);
    // destroy the node at the given pointer. This pointer must point to a Node owned by the graph.
    void destroy_node(Node* node);
    // destroy the node with the given ID.
    void destroy_node(int64_t id);
    bool has_node(int64_t id);
    bool has_node(Node* node);
    bool has_node(Node& node);
    void for_each_node(function<void(Node*)> lambda);
    void for_each_node_parallel(function<void(Node*)> lambda);

    // is the graph empty?
    bool empty(void);

    // remove nodes with no sequence
    // these are created in some cases during the process of graph construction
    void remove_null_nodes(void);
    // remove a node but connect all of its predecessor and successor nodes with new edges 
    void remove_node_forwarding_edges(Node* node);
    // remove null nodes but connect predecessors and successors, preserving structure
    void remove_null_nodes_forwarding_edges(void);

    // remove edges for which one of the nodes is not present
    void remove_orphan_edges(void);
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

    // edges
    // If the given edge cannot be created (because it would create a cycle), returns null.
    // If the given edge already exists, returns the existing edge.
    Edge* create_edge(Node* from, Node* to, bool from_start = false, bool to_end = false);
    Edge* create_edge(int64_t from, int64_t to, bool from_start = false, bool to_end = false);
    // Makes a left-to-right edge from the left NodeTraversal to the right one, respecting orientations.
    Edge* create_edge(NodeTraversal left, NodeTraversal right);
    // TODO: version that takes ID, bool pairs?
    
    // This can take nodes in any order
    Edge* get_edge(int64_t node1, int64_t node2);
    // destroy the edge at the given pointer. This pointer must point to an edge owned by the graph.
    void destroy_edge(Edge* edge);
    // destroy the edge between the nodes with the given IDs. These IDs can be in either order.
    void destroy_edge(int64_t node1, int64_t node2);
    // remove an edge from the node side indexes, so it doesn't show up when you
    // ask for the edges connected to the side of a node. Makes the edge
    // untraversable until the indexes are rebuilt.
    void unindex_edge_by_node_sides(int64_t node1, int64_t node2);
    void unindex_edge_by_node_sides(Edge* edge);
    // Get the edge between the given nodes. Sicne only one edge can connect a
    // pair of nodes (due to acyclicity), these nodes can be in either order.
    bool has_edge(int64_t node1, int64_t node2);
    bool has_edge(Edge* edge);
    bool has_edge(Edge& edge);
    void for_each_edge(function<void(Edge*)> lambda);
    void for_each_edge_parallel(function<void(Edge*)> lambda);

    // connect node -> nodes
    // Connects from the right side of the first to the left side of the second
    void connect_node_to_nodes(NodeTraversal node, vector<NodeTraversal>& nodes);
    void connect_node_to_nodes(Node* node, vector<Node*>& nodes);
    // connect nodes -> node
    // Connects from the right side of the first to the left side of the second
    void connect_nodes_to_node(vector<NodeTraversal>& nodes, NodeTraversal node);
    void connect_nodes_to_node(vector<Node*>& nodes, Node* node);

    // utilities
    // These only work on forward nodes.
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_path(map<long, int64_t>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    void to_dot(ostream& out, vector<Alignment> alignments = {});
    void to_gfa(ostream& out);
    bool is_valid(void);

    // topologically orders nodes
    void sort(void);
    // helper function, not really meant for external use
    void topological_sort(deque<NodeTraversal>& l);
    void swap_nodes(Node* a, Node* b);

    Alignment& align(Alignment& alignment);
    Alignment align(string& sequence);
    //Alignment& align(Alignment& alignment);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

    // returns all node-crossing paths with up to length across node boundaries
    // considers each node in forward orientation to produce the kpaths around it
    void for_each_kpath(int k, int edge_max,
                        function<void(NodeTraversal)> handle_prev_maxed,
                        function<void(NodeTraversal)> handle_next_maxed,
                        function<void(NodeTraversal,list<NodeTraversal>&)> lambda);
    void for_each_kpath_parallel(int k, int edge_max,
                                 function<void(NodeTraversal)> handle_prev_maxed,
                                 function<void(NodeTraversal)> handle_next_maxed,
                                 function<void(NodeTraversal,list<NodeTraversal>&)> lambda);
    void for_each_kpath(int k, int edge_max,
                        function<void(NodeTraversal)> handle_prev_maxed,
                        function<void(NodeTraversal)> handle_next_maxed,
                        function<void(Node*,Path&)> lambda);
    void for_each_kpath_parallel(int k, int edge_max,
                                 function<void(NodeTraversal)> handle_prev_maxed,
                                 function<void(NodeTraversal)> handle_next_maxed,
                                 function<void(Node*,Path&)> lambda);
    void for_each_kpath_of_node(Node* node, int k, int edge_max,
                                function<void(NodeTraversal)> handle_prev_maxed,
                                function<void(NodeTraversal)> handle_next_maxed,
                                function<void(NodeTraversal,list<NodeTraversal>&)> lambda);
    void for_each_kpath_of_node(Node* n, int k, int edge_max,
                                function<void(NodeTraversal)> handle_prev_maxed,
                                function<void(NodeTraversal)> handle_next_maxed,
                                function<void(Node*,Path&)> lambda);

    void kpaths(set<list<NodeTraversal> >& paths, int length, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    void kpaths(vector<Path>& paths, int length, int edge_max,
                function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);

    void kpaths_of_node(Node* node, set<list<NodeTraversal> >& paths,
                        int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    void kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    void kpaths_of_node(int64_t node_id, vector<Path>& paths, int length, int edge_max,
                        function<void(NodeTraversal)> prev_maxed, function<void(NodeTraversal)> next_maxed);
    // Given an oriented start node, a length in bp, a maximum number of edges
    // to cross, and a stack of nodes visited so far, fill in the set of paths
    // with all the paths starting at the oriented start node and going left no
    // longer than the specified length, calling maxed_nodes on nodes which
    // can't be visited due to the edge-crossing limit. Produces paths ending
    // with the specified node.
    // TODO: postfix should not be (potentially) copied on every call.
    void prev_kpaths_from_node(NodeTraversal node, int length, int edge_max,
                               list<NodeTraversal> postfix, set<list<NodeTraversal> >& paths,
                               function<void(NodeTraversal)>& maxed_nodes);
    // Do the same as prec_kpaths_from_node, except going right, producing a path starting with the specified node.
    void next_kpaths_from_node(NodeTraversal node, int length, int edge_max,
                               list<NodeTraversal> prefix, set<list<NodeTraversal> >& paths,
                               function<void(NodeTraversal)>& maxed_nodes);

    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(int64_t from, int64_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    string path_sequence(const Path& path);

    // traversal
    // Get the nodes attached to the left side of the given NodeTraversal, in their proper orientations.
    void nodes_prev(NodeTraversal n, vector<NodeTraversal>& nodes);
    // Get the nodes attached to the right side of the given NodeTraversal, in their proper orientations.
    void nodes_next(NodeTraversal n, vector<NodeTraversal>& nodes);
    // Count the nodes attached to the left side of the given NodeTraversal
    int node_count_prev(NodeTraversal n);
    // Count the nodes attached to the right side of the given NodeTraversal
    int node_count_next(NodeTraversal n);

    // paths
    Path create_path(const list<NodeTraversal>& nodes);
    Path create_path(const vector<NodeTraversal>& nodes);
    string path_string(const list<NodeTraversal>& nodes);
    // Assumes the path covers the entirety of any nodes visited. Handles backward nodes.
    string path_string(Path& path);
    void expand_path(const list<NodeTraversal>& path, vector<NodeTraversal>& expanded);
    // Fill in the node_start map with the first index along the path at which each node appears.
    // Caller is responsible for dealing with orientations.
    void node_starts_in_path(const list<NodeTraversal>& path,
                             map<Node*, int>& node_start);

    // kmers
    void for_each_kmer_parallel(int kmer_size,
                                int edge_max,
                                function<void(string&, NodeTraversal, int, list<NodeTraversal>&, VG&)> lambda,
                                int stride = 1,
                                bool allow_dups = false,
                                bool allow_negatives = false);
    void for_each_kmer(int kmer_size,
                       int edge_max,
                       function<void(string&, NodeTraversal, int, list<NodeTraversal>&, VG&)> lambda,
                       int stride = 1,
                       bool allow_dups = false,
                       bool allow_negatives = false);

    // for gcsa2. For the given kmer of the given length starting at the given
    // offset into the given Node along the given path, fill in prev_chars with
    // the characters that preceed it, next_chars with the characters that
    // follow it, and next_positions with the ((node ID, orientation), offset)
    // pairs of the places you can go next (from the right end of the kmer).
    // Refuses to follow more than edge_max edges. Offsets are in the path
    // orientation.
    void kmer_context(string& kmer,
                      int kmer_size,
                      int edge_max,
                      list<NodeTraversal>& path,
                      NodeTraversal node,
                      int32_t offset,
                      set<char>& prev_chars,
                      set<char>& next_chars,
                      set<pair<pair<int64_t, bool>, int32_t> >& next_positions);
    // for pruning graph prior to indexing with gcsa2
    // takes all nodes that would introduce paths of > edge_max edge crossings, removes them, and links their neighbors to
    // head_node or tail_node depending on which direction the path extension was stopped
    void prune_complex(int path_length, int edge_max, Node* head_node, Node* tail_node);

private:
    void _for_each_kmer(int kmer_size,
                        int edge_max,
                        function<void(string&, NodeTraversal, int, list<NodeTraversal>&, VG&)> lambda,
                        bool parallel,
                        int stride,
                        bool allow_dups,
                        bool allow_negatives);


public:

    // reads
    // note that even if either_strand is false, having backward nodes in the
    // graph will result in some reads from the global reverse strand.
    string random_read(int length, mt19937& rng, int64_t min_id, int64_t max_id, bool either_strand);

    // subgraphs
    void disjoint_subgraphs(list<VG>& subgraphs);
    // Get the head nodes (nodes with edges only to their right sides). These are required to be oriented forward.
    void head_nodes(vector<Node*>& nodes);
    vector<Node*> head_nodes(void);
    // Get the tail nodes (nodes with edges only to their left sides). These are required to be oriented forward.
    vector<Node*> tail_nodes(void);
    void tail_nodes(vector<Node*>& nodes);
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    // join head nodes of graph to common null node, creating a new single head.
    Node* join_heads(void);
    // or heads and tails to common new single head or tail.
    void join_heads(Node* node);
    void join_tails(Node* node);

    // add singular head and tail null nodes to graph
    void wrap_with_null_nodes(void);
    
    // Connect all existing head nodes to the given head node, and all existing
    // tail nodes to the given tail node. If either is null, it is created with
    // the specified length from the appropriate specified start/stop character,
    // and the corresponding pointer is updated. Used to prepare for indexing
    // with GCSA; length must be at least the GCSA kmer length to be used. Nodes
    // created here are owned by the graph, and will be deleted when the VG
    // object is deleted.
    void add_start_and_end_markers(int length, char start_char, char end_char,
                                   Node*& head_node, Node*& tail_node,
                                   int64_t head_id = 0, int64_t tail_id = 0);

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
        map<long, set<vcflib::VariantAllele> >* alleles;
        string seq;
        string name;
        Plan(VG* g,
             map<long, set<vcflib::VariantAllele> >* a,
             string s,
             string n)
            : graph(g)
            , alleles(a)
            , seq(s)
            , name(n) { };
        ~Plan(void) { delete alleles; }
    };


private:

    void init(void); // setup, ensures that gssw == NULL on startup
    // placeholders for empty
    vector<int64_t> empty_ids;
    vector<pair<int64_t, bool>> empty_edge_ends;
};

// utility functions

bool allATGC(string& s);
void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
string cigar_string(vector<pair<int, char> >& cigar);
string mapping_string(const string& source, const Mapping& mapping);
void divide_invariant_mapping(Mapping& orig, Mapping& left, Mapping& right, int offset, Node* nl, Node* nr);

} // end namespace vg

#endif
